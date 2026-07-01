"""
================================================================================
脚本名称: run_preprocessing.py
功能概述: 数据预处理与 Cell2location 反卷积（单切片 / 联合分析两种模式）—— Step 1
================================================================================

【整体任务说明】
    本脚本是分析流程的 Step 1，串联 preprocessing.py 模块中的所有核心函数。
    支持两种运行模式：

    ━━━ 模式 A：单切片分析（默认 / 历史兼容）━━━
    对单张 Visium 切片（如 CHC20）完成预处理与反卷积，可附带一张可选验证切片。
      Step 1  加载 scRNA-seq 参考数据，执行基础质控；
      Step 2  加载主分析 Visium 空间数据；
      Step 3  取 scRNA 与切片的共享基因子集；
      Step 4  提取 T/NK 细胞亚群，将指定聚类簇重注释为 Treg；
              绘制横版 Treg 标志基因气泡图（dotplot）；
      Step 5  训练 RegressionModel（参考签名学习）；
      Step 6  Cell2location 空间建模（反卷积），保存 adata_vis_post.h5ad；
      Step 7  （可选）对第二张验证切片独立训练并运行空间建模；
      Step 8  （可选）生成多切片一致性对比图；
      Step 9  保存 scRNA AnnData 和共享基因列表。

    ━━━ 模式 B：联合分析（HCC4R + CHC20，joint_mode=True）━━━
    确认两数据集 batch effect 极小后，用同一套 scRNA 参考训练一次 RegressionModel，
    得到统一的 cell_state_df，分别对 HCC4R 和 CHC20 做 Cell2location 反卷积，
    然后合并两张切片的 spot 级结果用于联合下游分析：
      Step 1  加载 scRNA-seq 参考数据；
      Step 2  加载 HCC4R + CHC20 两张 Visium 空间数据；
      Step 3  取 scRNA 与两切片的三路共同基因子集（确保反卷积基因空间完全一致）；
      Step 4  提取 T/NK 细胞亚群，重注释 Treg；
      Step 5  用共同基因子集训练统一 RegressionModel（一次训练，共享签名）；
      Step 6  用统一签名分别对 HCC4R 和 CHC20 做 Cell2location 反卷积；
      Step 7  合并两切片的 spot 级反卷积结果（加 sample 列标记来源），
              保存 adata_vis_post_joint.h5ad 供联合 niche 分析使用；
              同时保留每张切片的单独 h5ad（adata_vis_post_HCC4R.h5ad 等）；
      Step 8  生成联合细胞组成对比图；
      Step 9  保存共享基因列表和 scRNA AnnData。

    【联合分析设计说明】
    批次效应极小时，联合分析可获得更大样本量（更多 spot），提升统计功效。
    下游 run_spatial_niche_analysis.py 在读取合并后的 adata_vis_post_joint.h5ad 时，
    需通过 adata.obs["sample"] 识别切片来源，并在建立 kNN 空间邻域时限制在
    同一切片内（通过 --per-sample-neighbors 参数激活），避免跨切片物理邻居错误。

【输入文件】
    data/scRNA_reference.h5ad    - scRNA-seq 参考数据（pre.py 生成）
    data/CHC20_Visium/           - CHC20 Visium Space Ranger 输出目录
    data/HCC4R/                  - HCC4R Visium Space Ranger 输出目录

【输出文件（联合分析模式，保存于 results/joint_HCC4R_CHC20/）】
    adata_vis_post_HCC4R.h5ad            - HCC4R 单独反卷积结果
    adata_vis_post_CHC20.h5ad            - CHC20 单独反卷积结果
    adata_vis_post_joint.h5ad            - 两切片合并结果（含 obs["sample"] 列）
    adata_vis_post.h5ad                  - 同上，标准名称副本（供下游默认读取）
    spot_cell_proportion_HCC4R.csv       - HCC4R spot 细胞类型比例
    spot_cell_proportion_CHC20.csv       - CHC20 spot 细胞类型比例
    spot_cell_proportion_joint.csv       - 合并比例表（含 sample 列）
    shared_genes_joint.txt               - 三路共享基因列表
    regression_training_history_joint.png - 统一 RegressionModel 训练曲线
    t_cell_dotplot_horizontal.png        - Treg 标志基因横版气泡图
    cross_slice_comparison/              - 两切片细胞组成对比图
    run_preprocessing.log                - 全流程运行日志

【依赖关系】
    上游：pre.py（生成 .h5ad 数据文件）
    下游：run_spatial_niche_analysis.py（读取 adata_vis_post.h5ad 或
          adata_vis_post_joint.h5ad 进行 Niche 分析）

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
        "Treg Marker Gene Expression Across T/NK Subclusters",
        fontsize=10, y=1.01,
    )
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info("横版 Treg 气泡图已保存: %s", save_path)


def plot_celltype_marker_heatmap(
    adata_sc,
    celltype_col: str,
    marker_genes: dict[str, list[str]],
    save_path: Path,
    dpi: int = 180,
) -> None:
    """
    绘制细胞类型 Marker 基因平均表达热图，风格仿照参考论文 Figure D（2.png）。

    图表特征：
      - X 轴：细胞类型（按生物学相关性排列）
      - Y 轴：marker 基因（按细胞类型分组排列）
      - 颜色：Z-score 标准化后的平均表达量，使用 Yellow-Black-Purple 渐变色
        （高表达 = 亮黄，低表达 = 深紫，零表达 = 黑色）
      - 与参考论文 2.png 配色完全一致

    【配色说明】
    使用深蓝-青绿-亮黄三色渐变（paper_ybp），与论文图 D 的 Expression
    配色一致：
      低值 → 深蓝色 (#253494)
      中值 → 青绿色 (#1FA187)
      高值 → 亮黄色 (#FDE725)

    参数
    ----
    adata_sc     : scRNA-seq AnnData（已归一化）
    celltype_col : 细胞类型列名（如 "final_celltype"）
    marker_genes : dict，键为细胞类型/组名，值为该类型的 marker 基因列表
                   例如：{"T cell": ["IL7R","CD3G","CD2"], "Myeloid": ["LYZ","AIF1"]}
    save_path    : 图片保存路径
    dpi          : 输出分辨率
    """
    import scipy.sparse as sp_mod
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.stats import zscore as scipy_zscore

    # 深蓝-青绿-亮黄渐变色（统一论文图连续信号配色）
    _PAPER_YBP = LinearSegmentedColormap.from_list(
        "paper_ybp",
        ["#253494", "#1FA187", "#FDE725"],
        N=256,
    )

    if celltype_col not in adata_sc.obs.columns:
        LOGGER.warning("Cell type column '%s' not found; skipping marker heatmap.", celltype_col)
        return

    # 收集所有 marker 基因（去重，按 marker_genes 中顺序排列）
    all_genes_ordered = []
    gene_to_group = {}
    for group, genes in marker_genes.items():
        for g in genes:
            if g not in gene_to_group:  # 防止重复
                all_genes_ordered.append(g)
                gene_to_group[g] = group

    # 过滤不在数据集中的基因
    available_genes = [g for g in all_genes_ordered if g in adata_sc.var_names]
    if not available_genes:
        LOGGER.warning("No marker genes found in adata_sc; skipping marker heatmap.")
        return

    # 提取表达矩阵（CP10K + log1p 归一化）
    gene_idx = [list(adata_sc.var_names).index(g) for g in available_genes]
    X = adata_sc.X
    if sp_mod.issparse(X):
        expr = np.asarray(X[:, gene_idx].todense())
    else:
        expr = np.array(X[:, gene_idx])

    # CP10K + log1p 归一化
    totals = expr.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    expr = np.log1p(expr / totals * 1e4)

    expr_df = pd.DataFrame(expr, index=adata_sc.obs_names, columns=available_genes)
    expr_df["celltype"] = adata_sc.obs[celltype_col].values

    # 按细胞类型分组计算均值
    mean_expr = expr_df.groupby("celltype")[available_genes].mean()

    # 对每个基因做 Z-score 标准化（跨细胞类型标准化，便于横向对比）
    mean_z = mean_expr.apply(lambda col: scipy_zscore(col) if col.std() > 0 else col, axis=0)
    mean_z = mean_z.fillna(0.0)

    # 固定细胞类型显示顺序，优先满足论文作图需求。
    preferred_order = [
        "T/NK",
        "Treg",
        "Myeloid",
        "B",
        "Endothelial",
        "Hepatocyte",
        "Fibroblast",
    ]

    def _canonical_celltype(name: str) -> str:
        normalized = name.lower().replace(" ", "").replace("-", "").replace("/", "")
        alias_map = {
            "tnk": "T/NK",
            "tcell": "T/NK",
            "nk": "T/NK",
            "treg": "Treg",
            "myeloid": "Myeloid",
            "b": "B",
            "bcell": "B",
            "bcells": "B",
            "endothelial": "Endothelial",
            "hepatocyte": "Hepatocyte",
            "malignant": "Hepatocyte",
            "fibroblast": "Fibroblast",
            "hsc": "Fibroblast",
        }
        return alias_map.get(normalized, name)

    canonical_to_ct: dict[str, str] = {}
    for ct in mean_z.index:
        canonical = _canonical_celltype(str(ct))
        canonical_to_ct.setdefault(canonical, ct)

    ct_order = [
        canonical_to_ct[canonical]
        for canonical in preferred_order
        if canonical in canonical_to_ct
    ]

    # 追加未匹配的细胞类型
    for ct in mean_z.index:
        if ct not in ct_order:
            ct_order.append(ct)

    heatmap_data = mean_z.loc[ct_order, available_genes].T  # 基因为行，细胞类型为列

    # 动态调整图片大小
    n_genes = len(available_genes)
    n_celltypes = len(ct_order)
    fig_height = max(6, n_genes * 0.28 + 2)
    fig_width  = max(8, n_celltypes * 0.9 + 3)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    im = ax.imshow(
        heatmap_data.values,
        aspect="auto",
        cmap=_PAPER_YBP,
        vmin=-2, vmax=2,
    )

    # 坐标轴刻度和标签
    ax.set_xticks(range(n_celltypes))
    ax.set_xticklabels(
        [f"{ct}" for ct in ct_order],
        rotation=45, ha="right", fontsize=8,
    )
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(available_genes, fontsize=7, fontstyle="italic")

    # 颜色条
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.04)
    cbar.set_label("Expression\n(Z-score)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title(
        "Cell Type Marker Gene Expression\n(Mean log-normalized, Z-score per gene)",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info("细胞类型 Marker 热图已保存: %s", save_path)


# ============================================================
# 辅助函数：Cell2location 缓存读取
# ============================================================

def _adata_has_regression_signatures(adata_sc) -> bool:
    return "means_per_cluster_mu_fg" in getattr(adata_sc, "varm", {})


def _load_cached_reference_signatures(
    output_dir: Path,
    sample_name: str,
    prefer_standard: bool = False,
) -> tuple[object | None, pd.DataFrame | None]:
    """Read cached RegressionModel posterior signatures from adata_sc_post*.h5ad."""
    sample_path = output_dir / f"adata_sc_post_{sample_name}.h5ad"
    standard_path = output_dir / "adata_sc_post.h5ad"
    candidates = [standard_path, sample_path] if prefer_standard else [sample_path]

    for path in candidates:
        if not path.exists():
            continue
        try:
            adata_sc_cached = sc.read_h5ad(path)
            if not _adata_has_regression_signatures(adata_sc_cached):
                LOGGER.info(
                    "Cached scRNA file lacks RegressionModel signatures, ignoring: %s",
                    path,
                )
                continue
            cell_state_df = sanitize_cell_state_df(extract_cell_state_df(adata_sc_cached))
            LOGGER.info(
                "Loaded cached RegressionModel signatures for %s from %s (%d genes x %d cell types).",
                sample_name,
                path,
                cell_state_df.shape[0],
                cell_state_df.shape[1],
            )
            return adata_sc_cached, cell_state_df
        except Exception as exc:
            LOGGER.warning("Failed to load cached RegressionModel signatures from %s: %s", path, exc)
    return None, None


def _load_cached_cell2location_slice(
    output_dir: Path,
    sample_name: str,
    prefer_standard: bool = False,
) -> tuple[pd.DataFrame, object] | None:
    """Read cached Cell2location posterior h5ad and spot proportion table if available."""
    sample_h5ad = output_dir / f"adata_vis_post_{sample_name}.h5ad"
    standard_h5ad = output_dir / "adata_vis_post.h5ad"
    h5ad_candidates = [standard_h5ad, sample_h5ad] if prefer_standard else [sample_h5ad]

    for h5ad_path in h5ad_candidates:
        if not h5ad_path.exists():
            continue
        try:
            adata_vis_post = sc.read_h5ad(h5ad_path)
            if (
                "means_cell_abundance_w_sf" not in adata_vis_post.obsm
                and "cell_abundance" not in adata_vis_post.obsm
            ):
                LOGGER.info(
                    "Cached spatial file lacks Cell2location abundance matrices, ignoring: %s",
                    h5ad_path,
                )
                continue
            if "cell_abundance" not in adata_vis_post.obsm:
                adata_vis_post.obsm["cell_abundance"] = adata_vis_post.obsm["means_cell_abundance_w_sf"]

            csv_path = output_dir / f"spot_cell_proportion_{sample_name}.csv"
            if csv_path.exists():
                proportion_df = pd.read_csv(csv_path)
            else:
                LOGGER.info(
                    "Cached spot proportion table missing for %s; recomputing from %s.",
                    sample_name,
                    h5ad_path,
                )
                proportion_df = compute_spot_cell_proportion(
                    adata_vis_post, abundance_key="means_cell_abundance_w_sf"
                )
                proportion_df.to_csv(csv_path, index=False)

            LOGGER.info(
                "Loaded cached Cell2location result for %s from %s; skipping spatial training.",
                sample_name,
                h5ad_path,
            )
            return proportion_df, adata_vis_post
        except Exception as exc:
            LOGGER.warning("Failed to load cached Cell2location result from %s: %s", h5ad_path, exc)
    return None


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
# 联合分析主流程（模式 B：共享 RegressionModel 反卷积，合并 spot 结果）
# ============================================================

def main_joint(
    path_scrna: Path,
    path_sample_a: Path,
    path_sample_b: Path,
    sample_a_name: str,
    sample_b_name: str,
    output_dir: Path,
) -> tuple:
    """
    联合分析主流程：HCC4R + CHC20（或任意两张 batch effect 极小的切片）。

    设计要点
    --------
    1. 用 scRNA 参考与两张切片的三路共同基因子集训练一次 RegressionModel，
       得到统一的 cell_state_df（不再各自独立训练）；
    2. 用统一签名分别对两张切片做 Cell2location 反卷积，保证细胞类型比例
       可直接比较（同一参考系下的估计值）；
    3. 合并两切片的 spot 级 AnnData，写入 obs["sample"] 列标记来源；
       保存为 adata_vis_post_joint.h5ad（供下游联合 niche 分析使用）；
    4. 空间邻域关系由 run_spatial_niche_analysis.py 的 --per-sample-neighbors
       参数保证"按切片内计算"，本函数不做邻域操作。

    参数
    ----
    path_scrna     : scRNA 参考 h5ad 路径
    path_sample_a  : 切片 A（如 HCC4R）Visium 目录
    path_sample_b  : 切片 B（如 CHC20）Visium 目录
    sample_a_name  : 切片 A 名称（文件命名用）
    sample_b_name  : 切片 B 名称（文件命名用）
    output_dir     : 结果输出根目录（如 results/joint_HCC4R_CHC20/）

    返回
    ----
    (adata_vis_a_post, adata_vis_b_post, adata_vis_joint, shared_genes_joint)
    """
    import anndata as ad

    global OUTPUT_DIR, LOGGER, LOG_PATH
    OUTPUT_DIR = output_dir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 重新初始化日志
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
    LOGGER.info("=== run_preprocessing.py [联合分析模式] ===")
    LOGGER.info("  样本 A : %s -> %s", sample_a_name, path_sample_a)
    LOGGER.info("  样本 B : %s -> %s", sample_b_name, path_sample_b)
    LOGGER.info("  输出目录: %s", OUTPUT_DIR)

    t0 = perf_counter()
    regression_model_cls, cell2location_cls = check_cell2location_available()

    # ── Step 1: 加载 scRNA-seq 参考数据 ─────────────────────────────────────────
    LOGGER.info("Step 1: Loading scRNA-seq reference data...")
    adata_sc = load_scrna_h5ad(path_scrna, t0=t0)

    # ── Step 2: 加载两张 Visium 空间数据 ────────────────────────────────────────
    LOGGER.info("Step 2: Loading Visium spatial data for %s and %s...", sample_a_name, sample_b_name)
    adata_vis_a = load_visium(path_sample_a, sample_name=sample_a_name)
    adata_vis_b = load_visium(path_sample_b, sample_name=sample_b_name)

    # ── Step 3: 三路共同基因对齐 ────────────────────────────────────────────────
    # 要求 scRNA × 切片A × 切片B 的三路交集，保证两切片共用同一基因空间
    LOGGER.info("Step 3: Computing three-way shared gene intersection...")
    shared_a, adata_sc_a, adata_vis_a_sh = align_shared_genes(adata_sc, adata_vis_a)
    shared_b, adata_sc_b, adata_vis_b_sh = align_shared_genes(adata_sc, adata_vis_b)

    # 取三路交集
    joint_genes = sorted(set(shared_a) & set(shared_b))
    if len(joint_genes) == 0:
        raise ValueError("No shared genes across scRNA, sample A, and sample B.")
    LOGGER.info(
        "Three-way shared genes: %d  (scRNA∩%s=%d, scRNA∩%s=%d)",
        len(joint_genes), sample_a_name, len(shared_a), sample_b_name, len(shared_b),
    )

    # 过滤到三路共同基因
    joint_gene_set = set(joint_genes)
    adata_sc_joint = adata_sc[:, [g for g in adata_sc.var_names if g in joint_gene_set]].copy()
    adata_vis_a_sh = adata_vis_a_sh[:, [g for g in adata_vis_a_sh.var_names if g in joint_gene_set]].copy()
    adata_vis_b_sh = adata_vis_b_sh[:, [g for g in adata_vis_b_sh.var_names if g in joint_gene_set]].copy()

    # ── Step 4: Treg 亚群鉴定与标注 ─────────────────────────────────────────────
    LOGGER.info("Step 4: Identifying Treg subcluster...")
    adata_t = subset_t_cells(adata_sc_joint)
    treg_markers = ["CD3D", "CD4", "FOXP3", "IL2RA"]
    available_treg_markers = [g for g in treg_markers if g in adata_t.var_names]
    missing_treg_markers = sorted(set(treg_markers) - set(available_treg_markers))
    if missing_treg_markers:
        LOGGER.warning("Missing Treg marker genes: %s", ", ".join(missing_treg_markers))

    if available_treg_markers:
        try:
            treg_dotplot_path = output_dir / "t_cell_dotplot_horizontal.png"
            plot_treg_dotplot_horizontal(
                adata_t,
                marker_genes=available_treg_markers,
                groupby="res.3",
                save_path=treg_dotplot_path,
            )
        except Exception as exc:
            LOGGER.warning("Treg dotplot failed (%s); falling back to vertical version.", exc)
            sc.pl.dotplot(adata_t, available_treg_markers, groupby="res.3", show=False)
            plt.savefig(output_dir / "t_cell_dotplot.png", dpi=150, bbox_inches="tight")
            plt.close()

    adata_sc_joint = assign_treg_label(adata_sc_joint, treg_clusters=("9", "9.0"))
    LOGGER.info(
        "Joint reference cell type counts:\n%s",
        adata_sc_joint.obs["final_celltype"].value_counts().to_string(),
    )

    # ── Step 4+: 绘制 scRNA-seq marker 热图（仿论文 Figure D 风格）──────────────
    LOGGER.info("Step 4+: Generating scRNA-seq cell type visualization plots (paper style)...")

    _celltype_col = (
        "final_celltype" if "final_celltype" in adata_sc_joint.obs.columns
        else (adata_sc_joint.obs.columns[0] if len(adata_sc_joint.obs.columns) > 0 else "celltype")
    )

    # 图：细胞类型 Marker 基因热图（仿论文 Figure D 风格：docs/plots/2.png）
    _CELLTYPE_MARKERS = {
        "T/NK":         ["IL7R", "CD3G", "CD2", "ITM2A", "CD3D"],
        "Treg":         ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
        "Myeloid":      ["LYZ", "AIF1", "RNASE1", "C1QB", "HLA-DRA"],
        "B":            ["B3GNT7", "MS4A1", "BANK1", "CD79A", "TNFRSF13C", "BCL11A"],
        "Endothelial":  ["PECAM1", "CDH5", "SPARCL1", "STC1", "SPARC", "TM4SF1"],
        "Hepatocyte":   ["APOA2", "ALB", "APOA1", "AMBP", "APOH", "TTR"],
        "Fibroblast":   ["RGS5", "COL1A1", "ACTA2", "PDGFRB"],
    }
    try:
        plot_celltype_marker_heatmap(
            adata_sc_joint,
            celltype_col=_celltype_col,
            marker_genes=_CELLTYPE_MARKERS,
            save_path=output_dir / "scrna_celltype_marker_heatmap.png",
        )
    except Exception as exc:
        LOGGER.warning("Marker heatmap failed: %s", exc)

    # ── Step 5: 统一 RegressionModel 训练（一次训练，两切片共享）────────────────
    LOGGER.info("Step 5: Training shared RegressionModel (joint reference)...")
    model_joint = setup_and_train_regression_model(adata_sc_joint, regression_model_cls)

    model_joint.plot_history(50)
    hist_path = output_dir / "regression_training_history_joint.png"
    plt.savefig(hist_path, dpi=150, bbox_inches="tight")
    plt.close()
    LOGGER.info("统一参考模型训练曲线已保存: %s", hist_path)

    adata_sc_joint = export_signatures(model_joint, adata_sc_joint)
    cell_state_df_joint = sanitize_cell_state_df(extract_cell_state_df(adata_sc_joint))
    LOGGER.info("Shared cell_state_df: %d genes × %d cell types", *cell_state_df_joint.shape)

    # ── Step 6: 用统一签名分别对两张切片做 Cell2location 反卷积 ─────────────────
    LOGGER.info("Step 6: Running Cell2location deconvolution for %s (joint ref)...", sample_a_name)
    prop_a, adata_vis_a_post = _run_cell2location_for_slice(
        adata_vis_raw=adata_vis_a_sh,
        cell_state_df=cell_state_df_joint,
        cell2location_cls=cell2location_cls,
        sample_name=sample_a_name,
        output_dir=output_dir,
    )

    LOGGER.info("Step 6: Running Cell2location deconvolution for %s (joint ref)...", sample_b_name)
    prop_b, adata_vis_b_post = _run_cell2location_for_slice(
        adata_vis_raw=adata_vis_b_sh,
        cell_state_df=cell_state_df_joint,
        cell2location_cls=cell2location_cls,
        sample_name=sample_b_name,
        output_dir=output_dir,
    )

    # ── Step 7: 合并两张切片的 spot 级结果（添加 sample 来源标记）────────────────
    LOGGER.info("Step 7: Merging spot-level results from %s and %s...", sample_a_name, sample_b_name)

    # 在 obs 中写入切片来源标记（供下游空间邻域分析按切片内计算）
    adata_vis_a_post.obs["sample"] = sample_a_name
    adata_vis_b_post.obs["sample"] = sample_b_name

    # 确保两切片使用同一基因集（应一致，保险起见取交集再对齐）
    common_vars = adata_vis_a_post.var_names.intersection(adata_vis_b_post.var_names)
    if len(common_vars) < len(adata_vis_a_post.var_names):
        LOGGER.warning(
            "After deconvolution, %d genes dropped from var intersection (%d → %d).",
            len(adata_vis_a_post.var_names) - len(common_vars),
            len(adata_vis_a_post.var_names),
            len(common_vars),
        )
    adata_a_aligned = adata_vis_a_post[:, common_vars].copy()
    adata_b_aligned = adata_vis_b_post[:, common_vars].copy()

    # obs_names 加前缀防止合并后重复
    import re as _re
    adata_a_aligned.obs_names = [f"{sample_a_name}_{n}" for n in adata_a_aligned.obs_names]
    adata_b_aligned.obs_names = [f"{sample_b_name}_{n}" for n in adata_b_aligned.obs_names]

    # 合并 obsm（means_cell_abundance_w_sf）
    def _concat_obsm_key(a, b, key):
        """
        将两个 AnnData 的同一 obsm 键拼接为一个 DataFrame/ndarray，
        行索引对应各自的 obs_names。
        """
        import numpy as _np
        val_a = a.obsm.get(key)
        val_b = b.obsm.get(key)
        if val_a is None or val_b is None:
            return None
        if isinstance(val_a, pd.DataFrame):
            return pd.concat([val_a, val_b], axis=0)
        return _np.concatenate([val_a, val_b], axis=0)

    adata_vis_joint = ad.concat(
        [adata_a_aligned, adata_b_aligned],
        axis=0,
        merge="same",
        uns_merge="same",
    )
    # concat 默认不拷贝 obsm，手动合并
    for key in ["means_cell_abundance_w_sf", "cell_abundance", "spatial"]:
        merged = _concat_obsm_key(adata_a_aligned, adata_b_aligned, key)
        if merged is not None:
            adata_vis_joint.obsm[key] = merged

    LOGGER.info(
        "Joint AnnData: %d spots (%s: %d, %s: %d), %d genes",
        adata_vis_joint.n_obs,
        sample_a_name, adata_a_aligned.n_obs,
        sample_b_name, adata_b_aligned.n_obs,
        adata_vis_joint.n_vars,
    )

    # 保存联合 h5ad
    joint_h5ad_path = output_dir / "adata_vis_post_joint.h5ad"
    adata_vis_joint.write_h5ad(joint_h5ad_path)
    LOGGER.info("联合 AnnData 已保存: %s", joint_h5ad_path)

    # 同时保存标准名称副本（供下游 run_spatial_niche_analysis.py 默认读取）
    std_path = output_dir / "adata_vis_post.h5ad"
    adata_vis_joint.write_h5ad(std_path)
    LOGGER.info("标准名称副本已保存: %s", std_path)

    # 合并比例表（含 sample 列）
    prop_a["sample"] = sample_a_name
    prop_b["sample"] = sample_b_name
    prop_joint = pd.concat([prop_a, prop_b], ignore_index=True)
    prop_joint.to_csv(output_dir / "spot_cell_proportion_joint.csv", index=False)
    LOGGER.info("联合细胞比例表已保存。")

    # ── Step 8: 生成联合细胞组成对比图 ──────────────────────────────────────────
    LOGGER.info("Step 8: Generating cross-slice comparison plots...")
    comparison_dir = output_dir / "cross_slice_comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    cell_types_to_compare = ["Hepatocyte", "Treg", "T/NK", "Myeloid", "Fibroblast", "B cell"]
    _compare_slices_generic(
        prop_a, prop_b, sample_a_name, sample_b_name,
        cell_types_to_compare, comparison_dir,
    )

    # ── Step 9: 保存共享基因列表和 scRNA AnnData ─────────────────────────────────
    adata_sc_joint.write_h5ad(output_dir / "adata_sc_post.h5ad")
    adata_sc_joint.write_h5ad(output_dir / f"adata_sc_post_joint.h5ad")
    (output_dir / "shared_genes_joint.txt").write_text(
        "\n".join(joint_genes), encoding="utf-8"
    )

    LOGGER.info("全部步骤完成。结果保存至: %s", output_dir)
    LOGGER.info("Run completed in %.2fs", perf_counter() - t0)
    LOGGER.info(
        "后续步骤：\n"
        "  Step 2 (联合空间生态位分析): "
        "python code/pipeline/run_spatial_niche_analysis.py "
        "--adata %s/adata_vis_post_joint.h5ad "
        "--out-dir %s/spatial_niche "
        "--per-sample-neighbors",
        output_dir, output_dir,
    )
    return adata_vis_a_post, adata_vis_b_post, adata_vis_joint, joint_genes


# ============================================================
# 主流程（模式 A：单切片 / 历史兼容）
# ============================================================

def main(
    path_scrna: Path = PATH_SCRNA_H5AD,
    path_sample1: Path = PATH_SAMPLE1,
    path_sample2: Path | None = PATH_SAMPLE2,
    sample1_name: str = SAMPLE1_NAME,
    sample2_name: str = SAMPLE2_NAME,
    output_dir: Path = OUTPUT_DIR,
    use_cache: bool = True,
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
    use_cache     : 若 results 中已有后验文件，则跳过耗时训练并复用缓存
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
    LOGGER.info("  use_cache    : %s", use_cache)

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

    # ── Step 4+: 绘制 scRNA-seq marker 热图（仿论文 Figure D 风格）──────────────
    LOGGER.info("Step 4+: Generating scRNA-seq cell type visualization plots (paper style)...")
    _celltype_col_s1 = (
        "final_celltype" if "final_celltype" in adata_sc_s1.obs.columns
        else (adata_sc_s1.obs.columns[0] if len(adata_sc_s1.obs.columns) > 0 else "celltype")
    )
    # 图：细胞类型 Marker 基因热图（仿论文 Figure D 风格：docs/plots/2.png）
    _CELLTYPE_MARKERS = {
        "T/NK":         ["IL7R", "CD3G", "CD2", "ITM2A", "CD3D"],
        "Treg":         ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
        "Myeloid":      ["LYZ", "AIF1", "RNASE1", "C1QB", "HLA-DRA"],
        "B":            ["B3GNT7", "MS4A1", "BANK1", "CD79A", "TNFRSF13C", "BCL11A"],
        "Endothelial":  ["PECAM1", "CDH5", "SPARCL1", "STC1", "SPARC", "TM4SF1"],
        "Hepatocyte":   ["APOA2", "ALB", "APOA1", "AMBP", "APOH", "TTR"],
        "Fibroblast":   ["RGS5", "COL1A1", "ACTA2", "PDGFRB"],
    }
    try:
        plot_celltype_marker_heatmap(
            adata_sc_s1,
            celltype_col=_celltype_col_s1,
            marker_genes=_CELLTYPE_MARKERS,
            save_path=output_dir / "scrna_celltype_marker_heatmap.png",
        )
    except Exception as exc:
        LOGGER.warning("Marker heatmap failed: %s", exc)

    # ── Step 5/6: 主分析切片 RegressionModel + Cell2location（优先复用缓存）────
    cached_s1 = (
        _load_cached_cell2location_slice(output_dir, sample1_name, prefer_standard=True)
        if use_cache else None
    )
    if cached_s1 is not None:
        prop_s1, adata_vis_s1_post = cached_s1
        standard_h5ad = output_dir / "adata_vis_post.h5ad"
        if not standard_h5ad.exists():
            adata_vis_s1_post.write_h5ad(standard_h5ad)
            LOGGER.info("%s cached result copied to standard path: %s", sample1_name, standard_h5ad)
    else:
        cached_adata_sc_s1, cell_state_df_s1 = (
            _load_cached_reference_signatures(output_dir, sample1_name, prefer_standard=True)
            if use_cache else (None, None)
        )
        if cached_adata_sc_s1 is not None and cell_state_df_s1 is not None:
            adata_sc_s1 = cached_adata_sc_s1
            LOGGER.info("Step 5: Reusing cached RegressionModel signatures for %s.", sample1_name)
        else:
            LOGGER.info("Step 5: Training RegressionModel for %s...", sample1_name)
            model_s1 = setup_and_train_regression_model(adata_sc_s1, regression_model_cls)

            model_s1.plot_history(50)
            history_path_s1 = output_dir / f"regression_training_history_{sample1_name}.png"
            plt.savefig(history_path_s1, dpi=150, bbox_inches="tight")
            plt.close()
            LOGGER.info("%s 单细胞模型训练曲线已保存: %s", sample1_name, history_path_s1)

            adata_sc_s1 = export_signatures(model_s1, adata_sc_s1)
            cell_state_df_s1 = sanitize_cell_state_df(extract_cell_state_df(adata_sc_s1))

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

    # ── Step 7: 验证切片独立建模（可选，优先复用缓存）─────────────────────────
    if adata_vis_s2_sh is not None:
        cached_s2 = (
            _load_cached_cell2location_slice(output_dir, sample2_name)
            if use_cache else None
        )
        if cached_s2 is not None:
            prop_s2, adata_vis_s2_post = cached_s2
        else:
            cached_adata_sc_s2, cell_state_df_s2 = (
                _load_cached_reference_signatures(output_dir, sample2_name)
                if use_cache else (None, None)
            )
            if cached_adata_sc_s2 is not None and cell_state_df_s2 is not None:
                adata_sc_s2 = cached_adata_sc_s2
                LOGGER.info("Step 7: Reusing cached RegressionModel signatures for %s.", sample2_name)
            else:
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
    if _adata_has_regression_signatures(adata_sc_s1):
        adata_sc_s1.write_h5ad(output_dir / "adata_sc_post.h5ad")
        adata_sc_s1.write_h5ad(output_dir / f"adata_sc_post_{sample1_name}.h5ad")
    else:
        LOGGER.info(
            "%s scRNA object has no RegressionModel signatures; keeping existing adata_sc_post cache unchanged.",
            sample1_name,
        )
    (output_dir / f"shared_genes_{sample1_name}.txt").write_text(
        "\n".join(shared_genes_s1), encoding="utf-8"
    )
    if adata_sc_s2 is not None:
        if _adata_has_regression_signatures(adata_sc_s2):
            adata_sc_s2.write_h5ad(output_dir / f"adata_sc_post_{sample2_name}.h5ad")
        else:
            LOGGER.info(
                "%s scRNA object has no RegressionModel signatures; keeping existing sample2 adata_sc_post cache unchanged.",
                sample2_name,
            )
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