"""
================================================================================
脚本名称: run_preprocessing.py
功能概述: 数据预处理与 Cell2location 反卷积 —— Step 1
================================================================================

【整体任务说明】
    本脚本执行以下核心步骤：
      1. 加载 scRNA-seq 参考数据，完成质控与 Treg 亚群注释；
      2. 绘制横版 Treg 标志基因气泡图（便于比较多个亚聚类）；
      3. 训练 Cell2location RegressionModel，导出细胞类型参考签名；
      4. 对 CHC20 Visium 切片运行 Cell2location 空间建模，获得细胞丰度估计；
      5. 对 CHC23 Visium 切片运行相同流程，生成多切片一致性验证结果；
      6. 保存反卷积后的 AnnData 和细胞比例表。

【多切片设计说明】
    CHC20 为主分析切片（用于 Step 2 空间 niche 分析）；
    CHC23 为独立验证切片，运行相同的 Cell2location 流程后：
      - 生成相同系列可视化图（细胞类型空间分布、niche 评分分布）；
      - 与 CHC20 进行关键指标对比（Treg 富集比例、niche_high 占比）；
      - 验证结论在不同患者/切片间的一致性，增强结果的普适性。

【输入文件】
    data/scRNA_reference.h5ad    - scRNA-seq 参考数据（pre.py 生成）
    data/chc20_visium.h5ad       - CHC20 Visium 空间转录组数据（pre.py 生成）
                                   或直接读取 data/CHC20_Visium/ Space Ranger 目录
    data/chc23_visium.h5ad       - CHC23 Visium 空间转录组数据（pre.py 生成）
                                   或直接读取 data/CHC23_Visium/ Space Ranger 目录

【输出文件】
    results/adata_vis_post.h5ad       - CHC20 Cell2location 反卷积后的空间 AnnData（主分析输入）
    results/adata_vis_chc23_post.h5ad - CHC23 Cell2location 反卷积后的空间 AnnData（验证切片）
    results/spot_cell_proportion.csv  - CHC20 每个 spot 的细胞类型比例表
    results/spot_cell_proportion_chc23.csv - CHC23 每个 spot 的细胞类型比例表
    results/figures/treg_bubble.png   - Treg 标志基因横版气泡图
    results/figures/cell2loc_chc20_*.png  - CHC20 细胞类型空间分布图（各细胞类型一张）
    results/figures/cell2loc_chc23_*.png  - CHC23 细胞类型空间分布图（各细胞类型一张）

【参考文献】
    - Kleshchevnikov et al., Nature Biotechnology, 2022 (Cell2location)
    - Luecken & Theis, Molecular Systems Biology, 2019 (scRNA-seq 最佳实践)
    - Williams et al., Genome Medicine, 2022 (Visium 分析最佳实践)
================================================================================
"""
from pathlib import Path
from time import perf_counter
import logging
import sys

import matplotlib
matplotlib.use("Agg")   # 非交互后端，适合服务器/脚本运行
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch

torch.set_float32_matmul_precision("medium")   # 降低浮点精度，提升训练速度

from preprocessing import (
    align_shared_genes,
    assign_treg_label,
    check_cell2location_available,
    compute_spot_cell_proportion,
    export_signatures,
    extract_cell_state_df,
    load_scrna_h5ad,
    load_visium,
    sanitize_cell_state_df,
    setup_and_train_regression_model,
    subset_t_cells,
)


# ============================================================
# 路径配置
# ============================================================
ROOT_DIR   = Path("/mnt/hdd/private/gyw/Code/txbb/Single_cell_train")
DATA_DIR   = ROOT_DIR / "data"
OUTPUT_DIR = ROOT_DIR / "results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PATH_SCRNA_H5AD = DATA_DIR / "scRNA_reference.h5ad"
PATH_CHC20      = DATA_DIR / "CHC20_Visium"
PATH_CHC23      = DATA_DIR / "CHC23_Visium"

# ============================================================
# 日志配置：同时写入文件和终端
# ============================================================
LOG_PATH = OUTPUT_DIR / "run_preprocessing.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_PATH, encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
LOGGER = logging.getLogger(__name__)


# ============================================================
# 辅助函数：横版气泡图
# ============================================================

def plot_treg_dotplot_horizontal(
    adata_t,
    marker_genes: list[str],
    groupby: str,
    save_path: Path,
    dpi: int = 150,
) -> None:
    """
    绘制横版 Treg 标志基因气泡图（dot plot）。

    【横版说明】
    scanpy 默认的 dotplot 为竖版（x 轴为基因，y 轴为聚类簇）。
    横版将布局转置：x 轴为聚类簇，y 轴为基因，
    更便于在论文中与其他横向图表并排排列，
    且当聚类簇数量较多时标签不会重叠。

    参数
    ----
    adata_t      : 含 T/NK 细胞亚群的 AnnData
    marker_genes : 标志基因列表（如 ["CD3D","CD4","FOXP3","IL2RA"]）
    groupby      : 分组列名（如 "res.3"）
    save_path    : 图片保存路径
    dpi          : 输出分辨率
    """
    # 使用 scanpy 的 DotPlot 对象，通过 swap_axes() 实现横版
    dp = sc.pl.DotPlot(
        adata_t,
        var_names=marker_genes,
        groupby=groupby,
    )
    dp.swap_axes()          # 转置：x 轴变为聚类簇，y 轴变为基因
    dp.style(
        dot_max=1.0,
        color_on="dot",
        cmap="Blues",
        dot_edge_lw=0.5,
    )
    ax_dict = dp.get_axes()
    fig = list(ax_dict.values())[0].get_figure()
    fig.suptitle(
        "Treg Marker Gene Expression Across T/NK Subclusters\n"
        "(Horizontal Dot Plot)",
        fontsize=10, y=1.01,
    )
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info("横版 Treg 气泡图已保存: %s", save_path)


# ============================================================
# 辅助函数：CHC23 独立验证
# ============================================================

def _run_cell2location_for_slice(
    adata_vis_raw,
    cell_state_df,
    cell2location_cls,
    sample_name: str,
    output_dir: Path,
) -> tuple[pd.DataFrame, object]:
    """
    对单张 Visium 切片运行 Cell2location 空间建模，并保存结果。

    参数
    ----
    adata_vis_raw     : 已通过 align_shared_genes 对齐共享基因的 Visium AnnData
    cell_state_df     : 来自 RegressionModel 的细胞类型参考签名矩阵
    cell2location_cls : Cell2location 模型类
    sample_name       : 切片名称（用于文件命名和日志）
    output_dir        : 结果输出目录

    返回
    ----
    (proportion_df, adata_vis_raw)
      - proportion_df: 行=spot，列=细胞类型，值=归一化比例（首列为 spot_id）
      - adata_vis_raw: 写入 Cell2location 后验结果后的 AnnData
    """
    LOGGER.info("[%s] Setting up Cell2location model...", sample_name)

    # 确保输入为原始 count 矩阵（非负整数），Cell2location 要求此格式
    adata_vis_raw.X = adata_vis_raw.layers["counts"].copy()
    if sp.issparse(adata_vis_raw.X):
        adata_vis_raw.X.data = np.round(adata_vis_raw.X.data)
        adata_vis_raw.X.data[adata_vis_raw.X.data < 0] = 0
    else:
        adata_vis_raw.X = np.round(adata_vis_raw.X)
        adata_vis_raw.X[adata_vis_raw.X < 0] = 0

    cell2location_cls.setup_anndata(adata_vis_raw)
    sp_model = cell2location_cls(
        adata_vis_raw,
        cell_state_df=cell_state_df,
        N_cells_per_location=30,   # 每个 spot 预期细胞总数（参考 Kleshchevnikov 2022）
        detection_alpha=20,        # 检测效率先验，控制 spot 间检测变异
    )
    sp_model.train(max_epochs=1000, batch_size=1000)

    # 保存训练曲线
    sp_model.plot_history(100)
    hist_path = output_dir / f"spatial_mapping_training_history_{sample_name}.png"
    plt.savefig(hist_path, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("[%s] 空间模型训练曲线已保存: %s", sample_name, hist_path)

    # 导出后验估计
    adata_vis_raw = sp_model.export_posterior(adata_vis_raw)
    adata_vis_raw.obsm["cell_abundance"] = adata_vis_raw.obsm["means_cell_abundance_w_sf"]

    # 保存 h5ad（供后续 niche 分析使用）
    h5ad_path = output_dir / f"adata_vis_post_{sample_name}.h5ad"
    adata_vis_raw.write_h5ad(h5ad_path)
    LOGGER.info("[%s] 反卷积 AnnData 已保存: %s", sample_name, h5ad_path)

    # 计算并保存细胞比例表
    proportion_df = compute_spot_cell_proportion(
        adata_vis_raw, abundance_key="means_cell_abundance_w_sf"
    )
    csv_path = output_dir / f"spot_cell_proportion_{sample_name}.csv"
    proportion_df.to_csv(csv_path, index=False)
    LOGGER.info("[%s] 细胞比例表已保存: %s", sample_name, csv_path)

    return proportion_df, adata_vis_raw


def _compare_slices(
    prop_chc20: pd.DataFrame,
    prop_chc23: pd.DataFrame,
    cell_types: list[str],
    output_dir: Path,
) -> None:
    """
    对比 CHC20 与 CHC23 切片的关键细胞类型分布，生成一致性验证图表。

    图表内容：
      1. 各细胞类型均值比例的并排条形图（CHC20 vs CHC23）；
      2. Treg 比例分布的核密度估计（KDE）对比图；
      3. Myeloid 和 Fibroblast 比例散布的箱线图对比。

    参数
    ----
    prop_chc20  : CHC20 细胞比例表
    prop_chc23  : CHC23 细胞比例表
    cell_types  : 用于对比的细胞类型列表
    output_dir  : 图表输出目录
    """
    # 只保留两个数据集中都有的细胞类型
    avail = [c for c in cell_types if c in prop_chc20.columns and c in prop_chc23.columns]
    if not avail:
        LOGGER.warning("No common cell types for cross-slice comparison.")
        return

    # ---- 图1：均值比例并排条形图 ----
    mean20 = prop_chc20[avail].mean()
    mean23 = prop_chc23[avail].mean()
    compare_df = pd.DataFrame({
        "CHC20": mean20,
        "CHC23": mean23,
    })
    fig, ax = plt.subplots(figsize=(max(6, len(avail) * 1.2), 4.5))
    compare_df.plot(kind="bar", ax=ax, color=["#5b9bd5", "#ed7d31"], edgecolor="none")
    ax.set_title("Mean Cell Type Proportion: CHC20 vs CHC23")
    ax.set_ylabel("Mean proportion")
    ax.set_xlabel("Cell type")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_dir / "cross_slice_mean_proportion_comparison.png", dpi=180)
    plt.close(fig)

    # ---- 图2：Treg 比例 KDE 对比 ----
    if "Treg" in avail:
        import matplotlib.patches as mpatches
        fig, ax = plt.subplots(figsize=(6, 4))
        prop_chc20["Treg"].plot.kde(ax=ax, color="#5b9bd5", label="CHC20")
        prop_chc23["Treg"].plot.kde(ax=ax, color="#ed7d31", label="CHC23")
        ax.set_xlabel("Treg proportion")
        ax.set_title("Treg Proportion Distribution: CHC20 vs CHC23")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "cross_slice_treg_distribution.png", dpi=180)
        plt.close(fig)

    # ---- 图3：多细胞类型箱线图对比 ----
    box_types = [c for c in ["Treg", "Myeloid", "Fibroblast", "Hepatocyte"] if c in avail]
    if box_types:
        records = []
        for ct in box_types:
            for val in prop_chc20[ct]:
                records.append({"cell_type": ct, "proportion": val, "slice": "CHC20"})
            for val in prop_chc23[ct]:
                records.append({"cell_type": ct, "proportion": val, "slice": "CHC23"})
        box_df = pd.DataFrame(records)
        fig, ax = plt.subplots(figsize=(max(6, len(box_types) * 2), 4.5))
        import seaborn as sns
        sns.boxplot(
            data=box_df, x="cell_type", y="proportion", hue="slice",
            palette={"CHC20": "#5b9bd5", "CHC23": "#ed7d31"},
            ax=ax, fliersize=2,
        )
        ax.set_title("Cell Type Proportion Comparison: CHC20 vs CHC23")
        ax.set_ylabel("Proportion")
        ax.set_xlabel("Cell type")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "cross_slice_celltype_boxplot.png", dpi=180)
        plt.close(fig)

    LOGGER.info("Cross-slice comparison plots saved to %s", output_dir)


# ============================================================
# 主流程
# ============================================================

def main():
    """
    主流程：数据预处理、Treg 注释、Cell2location 反卷积（CHC20 + CHC23）。
    """
    t0 = perf_counter()

    regression_model_cls, cell2location_cls = check_cell2location_available()

    # ── Step 1: 加载 scRNA-seq 参考数据 ─────────────────────────────────────────
    LOGGER.info("Step 1: Loading scRNA-seq reference data...")
    adata_sc = load_scrna_h5ad(PATH_SCRNA_H5AD, t0=t0)

    # ── Step 2: 加载 Visium 空间数据（CHC20 和 CHC23） ──────────────────────────
    LOGGER.info("Step 2: Loading Visium spatial data...")
    adata_vis_chc20 = load_visium(PATH_CHC20, sample_name="CHC20")
    adata_vis_chc23 = load_visium(PATH_CHC23, sample_name="CHC23")

    # ── Step 3: 基因对齐（主分析与验证集各自独立）──────────────────────────────
    LOGGER.info("Step 3: Aligning shared genes separately for CHC20 and CHC23...")
    shared_genes_chc20, adata_sc_chc20, adata_vis_chc20_sh = align_shared_genes(
        adata_sc, adata_vis_chc20
    )
    shared_genes_chc23, adata_sc_chc23, adata_vis_chc23_sh = align_shared_genes(
        adata_sc, adata_vis_chc23
    )
    if len(shared_genes_chc20) == 0:
        raise ValueError("No shared genes found between scRNA and CHC20.")
    if len(shared_genes_chc23) == 0:
        raise ValueError("No shared genes found between scRNA and CHC23.")
    LOGGER.info("Shared genes (scRNA ∩ CHC20): %d", len(shared_genes_chc20))
    LOGGER.info("Shared genes (scRNA ∩ CHC23): %d", len(shared_genes_chc23))

    # ── Step 4: Treg 亚群鉴定与标注 ─────────────────────────────────────────────
    LOGGER.info("Step 4: Identifying Treg subcluster...")
    adata_t = subset_t_cells(adata_sc)

    # 绘制横版 Treg 标志基因气泡图（使用原始 scRNA，避免空间基因过滤影响 marker 展示）
    treg_dotplot_path = OUTPUT_DIR / "t_cell_dotplot_horizontal.png"
    treg_markers = ["CD3D", "CD4", "FOXP3", "IL2RA"]
    available_treg_markers = [gene for gene in treg_markers if gene in adata_t.var_names]
    missing_treg_markers = sorted(set(treg_markers) - set(available_treg_markers))
    if missing_treg_markers:
        LOGGER.warning(
            "Missing Treg marker genes in scRNA data, skipped in dotplot: %s",
            ", ".join(missing_treg_markers),
        )

    if available_treg_markers:
        try:
            plot_treg_dotplot_horizontal(
                adata_t,
                marker_genes=available_treg_markers,
                groupby="res.3",
                save_path=treg_dotplot_path,
            )
        except Exception as exc:
            LOGGER.warning(
                "Horizontal dotplot failed (%s); falling back to vertical version.", exc
            )
            # 回退：竖版
            sc.pl.dotplot(
                adata_t,
                available_treg_markers,
                groupby="res.3",
                show=False,
            )
            plt.savefig(OUTPUT_DIR / "t_cell_dotplot.png", dpi=150, bbox_inches="tight")
            plt.close()
    else:
        LOGGER.warning("No Treg marker genes available; skipping Treg dotplot.")

    adata_sc_chc20 = assign_treg_label(adata_sc_chc20, treg_clusters=("9", "9.0"))
    adata_sc_chc23 = assign_treg_label(adata_sc_chc23, treg_clusters=("9", "9.0"))
    LOGGER.info(
        "CHC20 reference final cell type counts:\n%s",
        adata_sc_chc20.obs["final_celltype"].value_counts().to_string(),
    )
    LOGGER.info(
        "CHC23 validation reference final cell type counts:\n%s",
        adata_sc_chc23.obs["final_celltype"].value_counts().to_string(),
    )

    # ── Step 5: 训练 CHC20 RegressionModel（主分析参考签名）─────────────────────
    LOGGER.info("Step 5: Training RegressionModel for CHC20 main analysis...")
    model_chc20 = setup_and_train_regression_model(adata_sc_chc20, regression_model_cls)

    model_chc20.plot_history(50)
    history_path_chc20 = OUTPUT_DIR / "regression_training_history_CHC20.png"
    plt.savefig(history_path_chc20, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("CHC20 单细胞模型训练曲线已保存: %s", history_path_chc20)

    adata_sc_chc20 = export_signatures(model_chc20, adata_sc_chc20)
    cell_state_df_chc20 = sanitize_cell_state_df(extract_cell_state_df(adata_sc_chc20))

    # ── Step 6: CHC20 Cell2location 空间建模 ────────────────────────────────────
    LOGGER.info("Step 6: Running Cell2location deconvolution for CHC20...")
    prop_chc20, adata_vis_chc20_post = _run_cell2location_for_slice(
        adata_vis_raw=adata_vis_chc20_sh,
        cell_state_df=cell_state_df_chc20,
        cell2location_cls=cell2location_cls,
        sample_name="CHC20",
        output_dir=OUTPUT_DIR,
    )
    # 同时保存标准名称（供 run_spatial_niche_analysis.py 直接读取）
    adata_vis_chc20_post.write_h5ad(OUTPUT_DIR / "adata_vis_post.h5ad")
    LOGGER.info("CHC20: Step 6 完成，adata_vis_post.h5ad 已保存为主分析结果。")

    # ── Step 7: CHC23 独立验证（单独训练参考签名 + 空间建模）────────────────────
    LOGGER.info("Step 7: Training RegressionModel for CHC23 validation...")
    # CHC23 单独使用更保守的早停参数：
    #   - 观察到 CHC23 默认配置下仅 ~90 轮就提前退出，后验丰度分布几乎均一（箱线图退化为一条线），
    #     说明模型未能真正区分细胞类型，早停触发过早。
    #   - early_stopping_patience: 30 → 50（给 ELBO 更长的观察窗口，避免在短暂平台期误判收敛）
    #   - early_stopping_min_delta: 1e-4 → 5e-5（对微小改善更敏感，防止过早停止）
    #   - max_epochs: 250 → 400（允许模型充分探索参数空间）
    model_chc23 = setup_and_train_regression_model(
        adata_sc_chc23,
        regression_model_cls,
        max_epochs=400,
        early_stopping_patience=50,
        early_stopping_min_delta=5e-5,
    )

    model_chc23.plot_history(50)
    history_path_chc23 = OUTPUT_DIR / "regression_training_history_CHC23.png"
    plt.savefig(history_path_chc23, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("CHC23 单细胞模型训练曲线已保存: %s", history_path_chc23)

    adata_sc_chc23 = export_signatures(model_chc23, adata_sc_chc23)
    cell_state_df_chc23 = sanitize_cell_state_df(extract_cell_state_df(adata_sc_chc23))

    LOGGER.info("Step 7: Running Cell2location deconvolution for CHC23 (validation)...")
    prop_chc23, adata_vis_chc23_post = _run_cell2location_for_slice(
        adata_vis_raw=adata_vis_chc23_sh,
        cell_state_df=cell_state_df_chc23,
        cell2location_cls=cell2location_cls,
        sample_name="CHC23",
        output_dir=OUTPUT_DIR,
    )

    # ── Step 8: 多切片一致性对比图 ──────────────────────────────────────────────
    LOGGER.info("Step 8: Generating cross-slice consistency comparison plots...")
    comparison_dir = OUTPUT_DIR / "cross_slice_comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    cell_types_to_compare = [
        "Hepatocyte", "Treg", "T/NK", "Myeloid", "Fibroblast", "B cell"
    ]
    _compare_slices(prop_chc20, prop_chc23, cell_types_to_compare, comparison_dir)

    # ── Step 9: 保存共享基因列表和 scRNA AnnData ─────────────────────────────────
    adata_sc_chc20.write_h5ad(OUTPUT_DIR / "adata_sc_post.h5ad")
    adata_sc_chc20.write_h5ad(OUTPUT_DIR / "adata_sc_post_CHC20.h5ad")
    adata_sc_chc23.write_h5ad(OUTPUT_DIR / "adata_sc_post_CHC23.h5ad")
    (OUTPUT_DIR / "shared_genes_CHC20.txt").write_text(
        "\n".join(shared_genes_chc20), encoding="utf-8"
    )
    (OUTPUT_DIR / "shared_genes_CHC23.txt").write_text(
        "\n".join(shared_genes_chc23), encoding="utf-8"
    )
    LOGGER.info("Step 1 完成（全部步骤）。结果保存至: %s", OUTPUT_DIR)
    LOGGER.info("Run completed in %.2fs", perf_counter() - t0)
    LOGGER.info(
        "后续步骤：\n"
        "  Step 2 (空间生态位分析): python code/run_spatial_niche_analysis.py\n"
        "  Step 3 (差异表达验证):  python code/run_de_analysis.py\n"
        "  Step 4-5 (TCGA生存分析): Rscript code/tcga_survival_analysis.R"
    )
    return adata_sc_chc20, adata_vis_chc20_post, shared_genes_chc20


if __name__ == "__main__":
    main()
