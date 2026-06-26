"""
================================================================================
脚本名称: run_preprocessing.py
功能概述: 数据预处理与 Cell2location 双切片反卷积 —— Step 1
================================================================================

【整体任务说明】
    本脚本是分析流程的 Step 1，串联 preprocessing.py 模块中的所有核心函数，
    对 CHC20（主分析）和 CHC23（独立验证）两张切片依次完成以下 9 个步骤：

      Step 1  加载 scRNA-seq 参考数据，执行基础质控（过滤低表达细胞/基因）；
      Step 2  加载 CHC20 / CHC23 Visium 空间数据，保留 is_tissue=1 的有效 spot；
      Step 3  取 scRNA 与两张 Visium 切片的共享基因子集，确保建模基因空间一致；
      Step 4  提取 T/NK 细胞亚群，将指定聚类簇（默认 "9" / "9.0"）重注释为 Treg；
              绘制横版 Treg 标志基因气泡图（dotplot），便于论文并排展示；
      Step 5  用 CHC20 scRNA-seq 训练 RegressionModel（第一阶段参考签名学习），
              保存 CHC20 训练曲线图；
      Step 6  用 CHC20 的参考签名运行 Cell2location 空间建模（第二阶段反卷积），
              获得每个 spot 的细胞丰度估计，保存主分析结果；
      Step 7  独立为 CHC23 训练 RegressionModel（采用更保守的早停参数：
              patience=50, min_delta=5e-5, max_epochs=400，防止过早收敛），
              并对 CHC23 运行空间建模，生成验证切片结果；
      Step 8  生成 CHC20 / CHC23 多切片一致性对比图（各细胞类型比例条形图），
              保存至 results/cross_slice_comparison/；
      Step 9  保存两张切片各自的 scRNA AnnData 和共享基因列表。

【多切片设计说明】
    CHC20 为主分析切片（用于 Step 2 空间 niche 分析）；
    CHC23 为独立验证切片（由 run_chc23_validation.py 消费），两者分别独立
    训练参考签名模型，以排除样本间批次效应对反卷积结果的干扰。

【输入文件】
    data/scRNA_reference.h5ad    - scRNA-seq 参考数据（pre.py 生成）
    data/chc20_visium.h5ad       - CHC20 Visium 空间转录组数据（pre.py 生成）
                                   或直接读取 data/CHC20_Visium/ Space Ranger 目录
    data/chc23_visium.h5ad       - CHC23 Visium 空间转录组数据（pre.py 生成）
                                   或直接读取 data/CHC23_Visium/ Space Ranger 目录

【输出文件】
    results/adata_vis_post.h5ad               - CHC20 反卷积后的空间 AnnData（主分析输入）
    results/adata_vis_post_CHC20.h5ad         - 同上，带 CHC20 后缀的副本
    results/adata_vis_post_CHC23.h5ad         - CHC23 反卷积后的空间 AnnData（验证切片）
    results/adata_sc_post.h5ad                - CHC20 scRNA-seq 参考数据（含后验签名）
    results/adata_sc_post_CHC20.h5ad          - 同上，带 CHC20 后缀的副本
    results/adata_sc_post_CHC23.h5ad          - CHC23 scRNA-seq 参考数据（含后验签名）
    results/spot_cell_proportion_CHC20.csv    - CHC20 每个 spot 的细胞类型比例表
    results/spot_cell_proportion_CHC23.csv    - CHC23 每个 spot 的细胞类型比例表
    results/shared_genes_CHC20.txt            - CHC20 scRNA × Visium 共享基因列表
    results/shared_genes_CHC23.txt            - CHC23 scRNA × Visium 共享基因列表
    results/regression_training_history_CHC20.png - CHC20 RegressionModel 训练曲线
    results/regression_training_history_CHC23.png - CHC23 RegressionModel 训练曲线
    results/treg_bubble.png                   - Treg 标志基因横版气泡图
    results/t_cell_dotplot.png                - T 细胞聚类气泡图（竖版回退版本）
    results/cross_slice_comparison/           - 多切片一致性对比图（细胞类型比例）
    results/run_preprocessing.log            - 全流程运行日志

【依赖关系】
    上游：pre.py（生成 .h5ad 数据文件）
    下游：run_spatial_niche_analysis.py（读取 adata_vis_post.h5ad 进行 Niche 分析）
          run_chc23_validation.py（读取 adata_vis_post_CHC23.h5ad 进行独立验证）
          colocation.py / run_de_analysis.py（读取细胞比例表和空间数据）

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
# 路径配置（动态推断，无硬编码路径）
# 直接运行此文件时使用默认 CHC20/CHC23；
# 通过 run_chc20.py / run_hcc4r.py / run_hcc6nr.py 调用时传入自定义路径。
# ============================================================
_REPO_ROOT = Path(__file__).resolve().parents[2]   # code/pipeline/ -> code/ -> 项目根
DATA_DIR   = _REPO_ROOT / "data"

# ── 默认值（直接运行时生效）
PATH_SCRNA_H5AD = DATA_DIR / "scRNA_reference.h5ad"
PATH_SAMPLE1    = DATA_DIR / "CHC20_Visium"    # 主分析切片
PATH_SAMPLE2    = DATA_DIR / "CHC23_Visium"    # 独立验证切片
SAMPLE1_NAME    = "CHC20"
SAMPLE2_NAME    = "CHC23"
OUTPUT_DIR      = _REPO_ROOT / "results" / SAMPLE1_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 日志配置：同时写入文件和终端（日志路径依赖 OUTPUT_DIR，在 main() 中重新初始化）
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


def _compare_slices_generic(
    prop_s1: pd.DataFrame,
    prop_s2: pd.DataFrame,
    name1: str,
    name2: str,
    cell_types: list[str],
    output_dir: Path,
) -> None:
    """
    对比两张 Visium 切片的关键细胞类型分布，生成一致性验证图表。
    支持任意切片名称（name1/name2），替代原来硬编码 CHC20/CHC23 的版本。

    图表内容：
      1. 各细胞类型均值比例的并排条形图；
      2. Treg 比例分布的核密度估计（KDE）对比图；
      3. Myeloid / Fibroblast / Hepatocyte / Treg 比例箱线图对比。

    参数
    ----
    prop_s1    : 切片 1 细胞比例表
    prop_s2    : 切片 2 细胞比例表
    name1      : 切片 1 标签，如 "CHC20" / "HCC4R"
    name2      : 切片 2 标签，如 "CHC23" / "HCC6NR"
    cell_types : 用于对比的细胞类型列表
    output_dir : 图表输出目录
    """
    avail = [c for c in cell_types if c in prop_s1.columns and c in prop_s2.columns]
    if not avail:
        LOGGER.warning("No common cell types for cross-slice comparison.")
        return

    colors = {"s1": "#5b9bd5", "s2": "#ed7d31"}

    # ---- 图1：均值比例并排条形图 ----
    compare_df = pd.DataFrame({
        name1: prop_s1[avail].mean(),
        name2: prop_s2[avail].mean(),
    })
    fig, ax = plt.subplots(figsize=(max(6, len(avail) * 1.2), 4.5))
    compare_df.plot(kind="bar", ax=ax, color=[colors["s1"], colors["s2"]], edgecolor="none")
    ax.set_title(f"Mean Cell Type Proportion: {name1} vs {name2}")
    ax.set_ylabel("Mean proportion")
    ax.set_xlabel("Cell type")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_dir / "cross_slice_mean_proportion_comparison.png", dpi=180)
    plt.close(fig)

    # ---- 图2：Treg 比例 KDE 对比 ----
    if "Treg" in avail:
        fig, ax = plt.subplots(figsize=(6, 4))
        prop_s1["Treg"].plot.kde(ax=ax, color=colors["s1"], label=name1)
        prop_s2["Treg"].plot.kde(ax=ax, color=colors["s2"], label=name2)
        ax.set_xlabel("Treg proportion")
        ax.set_title(f"Treg Proportion Distribution: {name1} vs {name2}")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "cross_slice_treg_distribution.png", dpi=180)
        plt.close(fig)

    # ---- 图3：多细胞类型箱线图对比 ----
    box_types = [c for c in ["Treg", "Myeloid", "Fibroblast", "Hepatocyte"] if c in avail]
    if box_types:
        import seaborn as sns
        records = []
        for ct in box_types:
            for val in prop_s1[ct]:
                records.append({"cell_type": ct, "proportion": val, "slice": name1})
            for val in prop_s2[ct]:
                records.append({"cell_type": ct, "proportion": val, "slice": name2})
        box_df = pd.DataFrame(records)
        fig, ax = plt.subplots(figsize=(max(6, len(box_types) * 2), 4.5))
        sns.boxplot(
            data=box_df, x="cell_type", y="proportion", hue="slice",
            palette={name1: colors["s1"], name2: colors["s2"]},
            ax=ax, fliersize=2,
        )
        ax.set_title(f"Cell Type Proportion Comparison: {name1} vs {name2}")
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

def main(
    path_scrna: Path = PATH_SCRNA_H5AD,
    path_sample1: Path = PATH_SAMPLE1,
    path_sample2: Path | None = PATH_SAMPLE2,
    sample1_name: str = SAMPLE1_NAME,
    sample2_name: str = SAMPLE2_NAME,
    output_dir: Path = OUTPUT_DIR,
):
    """
    主流程：数据预处理、Treg 注释、Cell2location 反卷积。

    参数
    ----
    path_scrna    : scRNA 参考 h5ad 路径
    path_sample1  : 主分析 Visium 目录（必填）
    path_sample2  : 验证 Visium 目录（可选；None 则跳过 Step 7/8）
    sample1_name  : 主分析切片名（用于文件命名），如 "CHC20" / "HCC4R"
    sample2_name  : 验证切片名，如 "CHC23" / "HCC6NR"
    output_dir    : 结果输出目录（不同数据集应传入不同路径）
    """
    global OUTPUT_DIR, LOGGER, LOG_PATH
    OUTPUT_DIR = output_dir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 重新初始化日志（输出到各数据集自己的目录）
    LOG_PATH = OUTPUT_DIR / "run_preprocessing.log"
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(LOG_PATH, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    LOGGER = logging.getLogger(__name__)
    LOGGER.info("=== run_preprocessing.py ===")
    LOGGER.info("  sample1      : %s -> %s", sample1_name, path_sample1)
    LOGGER.info("  sample2      : %s -> %s", sample2_name, path_sample2)
    LOGGER.info("  output_dir   : %s", OUTPUT_DIR)

    t0 = perf_counter()

    regression_model_cls, cell2location_cls = check_cell2location_available()

    # ── Step 1: 加载 scRNA-seq 参考数据 ─────────────────────────────────────────
    LOGGER.info("Step 1: Loading scRNA-seq reference data...")
    adata_sc = load_scrna_h5ad(path_scrna, t0=t0)

    # ── Step 2: 加载 Visium 空间数据 ────────────────────────────────────────────
    LOGGER.info("Step 2: Loading Visium spatial data...")
    adata_vis_s1 = load_visium(path_sample1, sample_name=sample1_name)
    adata_vis_s2 = load_visium(path_sample2, sample_name=sample2_name) if path_sample2 else None

    # ── Step 3: 基因对齐 ────────────────────────────────────────────────────────
    LOGGER.info("Step 3: Aligning shared genes...")
    shared_genes_s1, adata_sc_s1, adata_vis_s1_sh = align_shared_genes(adata_sc, adata_vis_s1)
    if len(shared_genes_s1) == 0:
        raise ValueError(f"No shared genes found between scRNA and {sample1_name}.")
    LOGGER.info("Shared genes (scRNA ∩ %s): %d", sample1_name, len(shared_genes_s1))

    if adata_vis_s2 is not None:
        shared_genes_s2, adata_sc_s2, adata_vis_s2_sh = align_shared_genes(adata_sc, adata_vis_s2)
        if len(shared_genes_s2) == 0:
            raise ValueError(f"No shared genes found between scRNA and {sample2_name}.")
        LOGGER.info("Shared genes (scRNA ∩ %s): %d", sample2_name, len(shared_genes_s2))
    else:
        shared_genes_s2 = adata_sc_s2 = adata_vis_s2_sh = None

    # ── Step 4: Treg 亚群鉴定与标注 ─────────────────────────────────────────────
    LOGGER.info("Step 4: Identifying Treg subcluster...")
    adata_t = subset_t_cells(adata_sc)

    treg_dotplot_path = output_dir / "t_cell_dotplot_horizontal.png"
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
            sc.pl.dotplot(
                adata_t,
                available_treg_markers,
                groupby="res.3",
                show=False,
            )
            plt.savefig(output_dir / "t_cell_dotplot.png", dpi=150, bbox_inches="tight")
            plt.close()
    else:
        LOGGER.warning("No Treg marker genes available; skipping Treg dotplot.")

    adata_sc_s1 = assign_treg_label(adata_sc_s1, treg_clusters=("9", "9.0"))
    LOGGER.info(
        "%s reference final cell type counts:\n%s",
        sample1_name,
        adata_sc_s1.obs["final_celltype"].value_counts().to_string(),
    )
    if adata_sc_s2 is not None:
        adata_sc_s2 = assign_treg_label(adata_sc_s2, treg_clusters=("9", "9.0"))
        LOGGER.info(
            "%s reference final cell type counts:\n%s",
            sample2_name,
            adata_sc_s2.obs["final_celltype"].value_counts().to_string(),
        )

    # ── Step 5: 训练主分析切片的 RegressionModel ──────────────────────────────────
    LOGGER.info("Step 5: Training RegressionModel for %s...", sample1_name)
    model_s1 = setup_and_train_regression_model(adata_sc_s1, regression_model_cls)

    model_s1.plot_history(50)
    history_path_s1 = output_dir / f"regression_training_history_{sample1_name}.png"
    plt.savefig(history_path_s1, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("%s 单细胞模型训练曲线已保存: %s", sample1_name, history_path_s1)

    adata_sc_s1 = export_signatures(model_s1, adata_sc_s1)
    cell_state_df_s1 = sanitize_cell_state_df(extract_cell_state_df(adata_sc_s1))

    # ── Step 6: 主分析切片的 Cell2location 空间建模 ─────────────────────────────
    LOGGER.info("Step 6: Running Cell2location deconvolution for %s...", sample1_name)
    prop_s1, adata_vis_s1_post = _run_cell2location_for_slice(
        adata_vis_raw=adata_vis_s1_sh,
        cell_state_df=cell_state_df_s1,
        cell2location_cls=cell2location_cls,
        sample_name=sample1_name,
        output_dir=output_dir,
    )
    # 同时保存标准名称（供 run_spatial_niche_analysis.py 直接读取）
    adata_vis_s1_post.write_h5ad(output_dir / "adata_vis_post.h5ad")
    LOGGER.info("%s: Step 6 完成，adata_vis_post.h5ad 已保存。", sample1_name)

    # ── Step 7: 验证切片独立建模（可选）─────────────────────────────────────────
    if adata_vis_s2_sh is not None:
        LOGGER.info("Step 7: Training RegressionModel for %s (validation)...", sample2_name)
        model_s2 = setup_and_train_regression_model(
            adata_sc_s2,
            regression_model_cls,
            max_epochs=400,
            early_stopping_patience=50,
            early_stopping_min_delta=5e-5,
        )
        model_s2.plot_history(50)
        history_path_s2 = output_dir / f"regression_training_history_{sample2_name}.png"
        plt.savefig(history_path_s2, dpi=150, bbox_inches="tight")
        plt.close()
        LOGGER.info("%s 单细胞模型训练曲线已保存: %s", sample2_name, history_path_s2)

        adata_sc_s2 = export_signatures(model_s2, adata_sc_s2)
        cell_state_df_s2 = sanitize_cell_state_df(extract_cell_state_df(adata_sc_s2))

        LOGGER.info("Step 7: Running Cell2location for %s (validation)...", sample2_name)
        prop_s2, adata_vis_s2_post = _run_cell2location_for_slice(
            adata_vis_raw=adata_vis_s2_sh,
            cell_state_df=cell_state_df_s2,
            cell2location_cls=cell2location_cls,
            sample_name=sample2_name,
            output_dir=output_dir,
        )
    else:
        LOGGER.info("Step 7: Skipped (no sample2 provided).")
        prop_s2 = adata_vis_s2_post = None
        shared_genes_s2 = None

    # ── Step 8: 多切片一致性对比图（仅在有验证切片时绘制）─────────────────────
    if prop_s2 is not None:
        LOGGER.info("Step 8: Generating cross-slice consistency comparison plots...")
        comparison_dir = output_dir / "cross_slice_comparison"
        comparison_dir.mkdir(parents=True, exist_ok=True)
        cell_types_to_compare = [
            "Hepatocyte", "Treg", "T/NK", "Myeloid", "Fibroblast", "B cell"
        ]
        # _compare_slices 内部标签固定为 CHC20/CHC23，这里用通用版本
        _compare_slices_generic(
            prop_s1, prop_s2, sample1_name, sample2_name,
            cell_types_to_compare, comparison_dir,
        )
    else:
        LOGGER.info("Step 8: Skipped (no sample2).")

    # ── Step 9: 保存共享基因列表和 scRNA AnnData ─────────────────────────────────
    adata_sc_s1.write_h5ad(output_dir / "adata_sc_post.h5ad")
    adata_sc_s1.write_h5ad(output_dir / f"adata_sc_post_{sample1_name}.h5ad")
    (output_dir / f"shared_genes_{sample1_name}.txt").write_text(
        "\n".join(shared_genes_s1), encoding="utf-8"
    )
    if adata_sc_s2 is not None:
        adata_sc_s2.write_h5ad(output_dir / f"adata_sc_post_{sample2_name}.h5ad")
        (output_dir / f"shared_genes_{sample2_name}.txt").write_text(
            "\n".join(shared_genes_s2), encoding="utf-8"
        )

    LOGGER.info("全部步骤完成。结果保存至: %s", output_dir)
    LOGGER.info("Run completed in %.2fs", perf_counter() - t0)
    LOGGER.info(
        "后续步骤：\n"
        "  Step 2 (空间生态位分析): python code/pipeline/run_spatial_niche_analysis.py "\
        "--adata %s/adata_vis_post.h5ad --out-dir %s/spatial_niche\n"
        "  Step 3 (差异表达验证):  python code/pipeline/run_de_analysis.py\n"
        "  Step 4-5 (TCGA生存分析): Rscript code/pipeline/tcga_survival_analysis.R",
        output_dir, output_dir,
    )
    return adata_sc_s1, adata_vis_s1_post, shared_genes_s1


if __name__ == "__main__":
    main()
