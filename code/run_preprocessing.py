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
    results/scrna_tsne_celltype.png           - scRNA-seq tSNE 细胞类型分布图（仿论文 Figure C）
                                               使用高饱和度离散色板，标注各细胞类型中心位置，
                                               与 docs/plots/1.png 风格一致
    results/scrna_celltype_marker_heatmap.png - 细胞类型 Marker 基因热图（仿论文 Figure D）
                                               Yellow-Black-Purple 渐变色，Z-score 标准化，
                                               与 docs/plots/2.png 风格一致
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
# 辅助函数：tSNE 细胞类型图（仿论文 Figure C 风格）
# ============================================================

def plot_scrna_tsne_celltype(
    adata_sc,
    celltype_col: str,
    cluster_col: str,
    save_path: Path,
    dpi: int = 180,
) -> None:
    """
    绘制 scRNA-seq tSNE 细胞类型分布图，风格仿照参考论文 Figure C（1.png）。

    图表特征：
      - 每个细胞在 tSNE 坐标系中绘制为一个点（s=6，半透明）
      - 每种细胞类型使用高饱和度固定颜色，与论文 tSNE 图配色风格一致
      - 每个细胞类型在图中标注聚类编号和类型名称
      - 右侧图例展示颜色-类型对应关系
      - 坐标轴标注 tSNE 1 / tSNE 2，带箭头指示方向

    【配色说明】
    参照论文图 C（tSNE 细胞类型聚类图）的色板设计：高饱和度离散色板，
    优先高对比度颜色，避免相邻细胞类型颜色混淆。

    参数
    ----
    adata_sc     : scRNA-seq AnnData，需包含 obsm["X_tsne"] 和 obs[celltype_col]
    celltype_col : 细胞类型列名（如 "final_celltype"）
    cluster_col  : 聚类编号列名（如 "res.3"），用于标注聚类编号（可为 None 跳过）
    save_path    : 图片保存路径
    dpi          : 输出分辨率
    """
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap

    # 检查 tSNE 坐标是否存在，若不存在则尝试 UMAP 或跳过
    if "X_tsne" in adata_sc.obsm:
        embed_key = "X_tsne"
        xlabel, ylabel = "tSNE 1", "tSNE 2"
    elif "X_umap" in adata_sc.obsm:
        embed_key = "X_umap"
        xlabel, ylabel = "UMAP 1", "UMAP 2"
    else:
        LOGGER.warning("Neither X_tsne nor X_umap found; skipping tSNE celltype plot.")
        return

    if celltype_col not in adata_sc.obs.columns:
        LOGGER.warning("Cell type column '%s' not found; skipping tSNE plot.", celltype_col)
        return

    # 高饱和度离散色板（仿论文 tSNE 图颜色风格）
    PAPER_CLUSTER_COLORS = [
        "#2ca02c",  # 绿色（T cell）
        "#9467bd",  # 紫色（Myeloid）
        "#1f77b4",  # 蓝色（Malignant）
        "#d62728",  # 红色
        "#ff7f0e",  # 橙色（NK）
        "#8c564b",  # 棕色
        "#e377c2",  # 粉色
        "#7f7f7f",  # 灰色
        "#bcbd22",  # 黄绿色
        "#17becf",  # 青色
        "#aec7e8",  # 浅蓝
        "#ffbb78",  # 浅橙
        "#98df8a",  # 浅绿
        "#ff9896",  # 浅红
        "#c5b0d5",  # 浅紫
        "#c49c94",  # 浅棕
        "#f7b6d2",  # 浅粉
        "#c7c7c7",  # 浅灰
        "#dbdb8d",  # 浅黄绿
        "#9edae5",  # 浅青
        "#393b79",  # 深蓝
        "#637939",  # 深绿
        "#8c6d31",  # 深棕
        "#843c39",  # 深红棕
    ]

    celltypes = sorted(adata_sc.obs[celltype_col].unique())
    color_map = {ct: PAPER_CLUSTER_COLORS[i % len(PAPER_CLUSTER_COLORS)]
                 for i, ct in enumerate(celltypes)}

    coords = adata_sc.obsm[embed_key]
    ct_labels = adata_sc.obs[celltype_col].values
    colors = [color_map.get(ct, "#bdbdbd") for ct in ct_labels]

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(
        coords[:, 0], coords[:, 1],
        c=colors, s=5, alpha=0.6, linewidths=0,
    )

    # 标注每种细胞类型的中心位置
    for ct in celltypes:
        mask = ct_labels == ct
        if mask.sum() == 0:
            continue
        cx = float(coords[mask, 0].mean())
        cy = float(coords[mask, 1].mean())
        ax.text(cx, cy, ct, fontsize=7.5, ha="center", va="center",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))

    # 图例（右侧）
    handles = [
        mpatches.Patch(facecolor=color_map[ct], label=ct, edgecolor="none")
        for ct in celltypes
    ]
    ax.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
        frameon=False, fontsize=7, ncol=1, title="Cell type",
        title_fontsize=8,
    )

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title("Single-cell Landscape: Cell Type Distribution", fontsize=11)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # 添加坐标轴箭头（仿论文风格）
    ax.annotate("", xy=(0.08, 0.0), xytext=(0.0, 0.0),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))
    ax.annotate("", xy=(0.0, 0.08), xytext=(0.0, 0.0),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2))

    fig.tight_layout()
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    LOGGER.info("tSNE 细胞类型图已保存: %s", save_path)


# ============================================================
# 辅助函数：细胞类型 Marker 基因热图（仿论文 Figure D 风格）
# ============================================================

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
    使用 Yellow-Black-Purple 三色渐变（paper_ybp），与论文图 D 的 Expression
    配色一致：
      低值 → 深紫色 (#3a0063)
      中值 → 纯黑色 (#000000)
      高值 → 亮黄色 (#f5e642)

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

    # Yellow-Black-Purple 渐变色（仿论文 Figure D 配色）
    _PAPER_YBP = LinearSegmentedColormap.from_list(
        "paper_ybp",
        ["#3a0063", "#000000", "#f5e642"],   # purple → black → yellow
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
    from scipy.stats import zscore as scipy_zscore
    mean_z = mean_expr.apply(lambda col: scipy_zscore(col) if col.std() > 0 else col, axis=0)
    mean_z = mean_z.fillna(0.0)

    # 确定细胞类型列的顺序（按 marker_genes 键中出现的顺序排列）
    # 先找 marker_genes 中出现的细胞类型，再追加其余类型
    ct_order = []
    for group in marker_genes:
        # 对 group 名尝试模糊匹配（忽略大小写和空格）
        for ct in mean_z.index:
            if group.lower().replace(" ", "") in ct.lower().replace(" ", "") and ct not in ct_order:
                ct_order.append(ct)
    # 追加未匹配的细胞类型
    for ct in mean_z.index:
        if ct not in ct_order:
            ct_order.append(ct)
    # 过滤不在 mean_z 中的类型
    ct_order = [ct for ct in ct_order if ct in mean_z.index]

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

    # ── Step 4+: 绘制 scRNA-seq 细胞类型图（仿论文 Figure C/D 风格）──────────────
    LOGGER.info("Step 4+: Generating scRNA-seq cell type visualization plots (paper style)...")

    # 图1：tSNE/UMAP 细胞类型分布图（仿论文 Figure C 风格：docs/plots/1.png）
    try:
        plot_scrna_tsne_celltype(
            adata_sc,
            celltype_col="final_celltype" if "final_celltype" in adata_sc.obs.columns
            else (adata_sc.obs.columns[0] if len(adata_sc.obs.columns) > 0 else "celltype"),
            cluster_col="res.3" if "res.3" in adata_sc.obs.columns else None,
            save_path=OUTPUT_DIR / "scrna_tsne_celltype.png",
        )
    except Exception as exc:
        LOGGER.warning("tSNE cell type plot failed: %s", exc)

    # 图2：细胞类型 Marker 基因热图（仿论文 Figure D 风格：docs/plots/2.png）
    # Marker 基因按照参考论文中各细胞类型的特异性基因组织
    CELLTYPE_MARKERS = {
        "T cell":       ["IL7R", "CD3G", "CD2", "ITM2A", "CD3D"],
        "Myeloid":      ["LYZ", "AIF1", "RNASE1", "C1QB", "HLA-DRA"],
        "NK":           ["GNLY", "GZMB", "KLRD1", "KLRF1"],
        "B cell":       ["B3GNT7", "MS4A1", "BANK1", "CD79A", "TNFRSF13C", "BCL11A"],
        "Malignant":    ["APOA2", "ALB", "APOA1", "AMBP", "APOH", "TTR"],
        "Endothelial":  ["PECAM1", "CDH5", "SPARCL1", "STC1", "SPARC", "TM4SF1"],
        "Epithelial":   ["INSR", "KRT18", "KRT19", "DEFB1", "CTSK", "EPCAM", "SOX4"],
        "Plasma cell":  ["JCHAIN", "TCF4", "TCL1A", "IGLL1", "MZB1", "IGLL5", "SSR4"],
        "HSC":          ["RGS5", "COL1A1", "ACTA2", "PDGFRB"],
    }
    try:
        plot_celltype_marker_heatmap(
            adata_sc,
            celltype_col="final_celltype" if "final_celltype" in adata_sc.obs.columns
            else (adata_sc.obs.columns[0] if len(adata_sc.obs.columns) > 0 else "celltype"),
            marker_genes=CELLTYPE_MARKERS,
            save_path=OUTPUT_DIR / "scrna_celltype_marker_heatmap.png",
        )
    except Exception as exc:
        LOGGER.warning("Marker heatmap failed: %s", exc)

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