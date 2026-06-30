"""
================================================================================
模块名称: paper_plot_functions.py
功能概述: 论文图表绘图函数库（与数据准备逻辑解耦）
================================================================================

【用途说明】
    本模块集中存放所有论文绘图函数，以及配色常量和通用辅助函数。
    修改绘图参数（颜色、字体、大小、样式等）只需在本文件中操作，
    无需接触数据准备逻辑（run_paper_figures.py）。

【函数列表】
    ── 通用工具 ──
      _save(fig, path, dpi)
      _expression_frame(adata) → pd.DataFrame
      _get_abundance(adata) → pd.DataFrame
      _proportions(abundance) → pd.DataFrame
      _wilcoxon_test(high_vals, low_vals) → float
      _p_label(p) → str

    ── Figure 1 ──
      plot_fig1A_workflow(out_dir, dpi)

    ── Figure 2 ──
      plot_fig2A_HE(visium_dir, spatial_df, out_dir, dpi)
      plot_fig2B_spatial_abundance(spatial_df, cell_types, out_dir, sample_name, dpi)
      plot_fig2C_coloc_heatmap(proportions, cell_types, out_dir, sample_name, dpi)

    ── Figure 3 ──
      plot_fig3A_spatial_domain(spatial_df, out_dir, domain_col, sample_name, dpi)
      plot_fig3B_domain_composition(spatial_df, cell_types, out_dir, domain_col, sample_name, dpi)
      plot_fig3C_score_violins(spatial_df, out_dir, domain_col, sample_name, dpi)

    ── Figure 4 ──
      plot_fig4A_volcano(niche_scores_df, adata, out_dir, sample_name, fc_threshold, fdr_threshold, dpi)
      plot_fig4B_niche_marker_heatmap(adata, spatial_df, ranked_genes_path, out_dir, top_n, domain_col, sample_name, dpi)
      plot_fig4C_pathway_bubble(adata, spatial_df, out_dir, domain_col, sample_name, dpi)

    ── Figure 5 ──
      compute_signature_score(adata, signature_genes, score_name, min_genes)
      plot_fig5A_signature_projection(chc20_adata, signature_score, out_dir, signature_name, matched_genes, dpi)
      plot_fig5B_validation_violin(hcc4r_score_df, chc20_score_df, out_dir, score_cols, dpi)
      plot_fig5C_sensitivity(hcc4r_niche_dir, out_dir, dpi)
      _plot_sensitivity_paper(sens_df, path, dpi)         # 内部子函数

【修改建议】
    - 修改全局配色：直接修改 PAPER_YBP / PAPER_CLUSTER_COLORS / DOMAIN_COLORS / CELLTYPE_COLORS
    - 修改 DPI：修改 PAPER_DPI 常量（或在调用时传入 dpi 参数）
    - 修改单张图的样式：直接找到对应函数，修改其中的 figsize / fontsize / alpha 等参数
    - 添加新图：在本文件中新增函数，再在 run_paper_figures.py 中调用即可
================================================================================
"""
from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.sparse as sp
import scipy.stats as scipy_stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ============================================================
# ★ 全局配色与样式常量（修改这里来统一调整论文风格）
# ============================================================

# 主色板：深蓝-青绿-亮黄，适合空间丰度、signature score 和表达热图
PAPER_YBP = LinearSegmentedColormap.from_list(
    "paper_ybp",
    ["#253494", "#1FA187", "#FDE725"],
    N=256,
)

# 高饱和度离散色板（用于细胞类型着色）
PAPER_CLUSTER_COLORS = [
    "#2ca02c", "#9467bd", "#1f77b4", "#d62728", "#ff7f0e",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]

# Spatial Domain 专用配色（键名与 niche CSV 中的 spatial_region 列值对应）
DOMAIN_COLORS = {
    "tumor_core":              "#d62728",
    "tumor_edge":              "#ff7f0e",
    "stroma_immune":           "#1f77b4",
    "immunosuppressive_niche": "#9467bd",
    "other":                   "#bdbdbd",
}

# 细胞类型配色
CELLTYPE_COLORS = {
    "Hepatocyte":  "#45496a",
    "Malignant":   "#7d8bae",
    "T/NK":        "#e5857b",
    "Treg":        "#f1b2b2",
    "Myeloid":     "#e8ccc7",
    "Fibroblast":  "#edce7a",
    "HSC":         "#4fb19d",
    "B cell":      "#9ac5e5",
    "Endothelial": "#b7bda0",
    "Plasma cell": "#45958e",
    "Epithelial":  "#fbe7ab",
}

# 全局 DPI（期刊投稿建议 ≥ 300）
PAPER_DPI: int = 300

# Domain 的显示名称映射
DOMAIN_DISPLAY_NAMES = {
    "tumor_core":    "Tumor Core",
    "tumor_edge":    "Tumor Edge",
    "stroma_immune": "Stroma/\nImmune",
    "other":         "Other",
}

# Spot 大小（散点图中每个 spot 的点大小）
SPOT_SIZE: float = 9.0

# 小提琴图 alpha 值
VIOLIN_ALPHA: float = 0.55


# ============================================================
# 通用工具函数
# ============================================================

def _save(fig: plt.Figure, path: Path, dpi: int = PAPER_DPI) -> None:
    """统一保存图片，设置白底，自动关闭 Figure 释放内存。"""
    fig.patch.set_facecolor("white")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logging.info("Saved: %s", path)


def _expression_frame(adata: ad.AnnData) -> pd.DataFrame:
    """从 AnnData 中提取 log1p-normalized (10k) 表达量 DataFrame。"""
    x = adata.X
    if sp.issparse(x):
        x = x.copy().astype(float).tocsr()
        totals = np.asarray(x.sum(axis=1)).ravel()
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        x = sp.diags(scale) @ x
        x.data = np.log1p(x.data)
        x = x.toarray()
    else:
        x = np.array(x, dtype=float)
        totals = x.sum(axis=1)
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        x = np.log1p(x * scale[:, None])
    return pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names)


def compute_signature_score(
    adata: ad.AnnData,
    signature_genes: list[str],
    score_name: str = "hcc4r_signature_score",
    min_genes: int = 5,
) -> tuple[pd.Series, list[str]]:
    """Score each spot by the mean log-normalized expression of signature genes.

    The score is deliberately computed with the same formula in discovery and
    validation slides, so Figure 5A and 5B are directly comparable.
    """
    clean_genes = [str(g).strip() for g in signature_genes if str(g).strip()]
    var_names = pd.Index(adata.var_names.astype(str))
    upper_lookup = {g.upper(): g for g in var_names}

    matched: list[str] = []
    seen: set[str] = set()
    for gene in clean_genes:
        hit = gene if gene in var_names else upper_lookup.get(gene.upper())
        if hit is not None and hit not in seen:
            matched.append(hit)
            seen.add(hit)

    if len(matched) < min_genes:
        logging.warning(
            "Only %d/%d signature genes found in AnnData; need at least %d.",
            len(matched), len(clean_genes), min_genes,
        )
        return pd.Series(dtype=float, name=score_name), matched

    expr = _expression_frame(adata[:, matched].copy())
    score = expr.mean(axis=1).rename(score_name)
    logging.info(
        "Computed %s using %d/%d signature genes.",
        score_name, len(matched), len(clean_genes),
    )
    return score, matched


def _get_abundance(adata: ad.AnnData) -> pd.DataFrame:
    """从 AnnData.obsm 中读取 cell2location 丰度矩阵。"""
    import re
    for key in ["means_cell_abundance_w_sf", "cell_abundance"]:
        if key in adata.obsm:
            raw = adata.obsm[key]
            if isinstance(raw, pd.DataFrame):
                df = raw.copy()
            else:
                factors = adata.uns.get("mod", {}).get("factor_names")
                df = pd.DataFrame(raw, index=adata.obs_names, columns=factors)
            df.index = adata.obs_names
            df.columns = [
                re.sub(r"^(means|q\d+|median)_?cell_abundance_w_sf_", "", str(c))
                for c in df.columns
            ]
            return df
    raise KeyError("No abundance key found in adata.obsm")


def _proportions(abundance: pd.DataFrame) -> pd.DataFrame:
    """将丰度矩阵归一化为各 spot 细胞类型比例（行和=1）。"""
    return abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)


def _wilcoxon_test(high_vals: np.ndarray, low_vals: np.ndarray) -> float:
    """单侧 Mann-Whitney U 检验（high > low），返回 p 值。"""
    try:
        _, p = scipy_stats.mannwhitneyu(
            high_vals[~np.isnan(high_vals)],
            low_vals[~np.isnan(low_vals)],
            alternative="greater",
        )
        return float(p)
    except Exception:
        return 1.0


def _p_label(p: float) -> str:
    """将 p 值转换为显著性标记符号。"""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


def _resolve_domain_col(spatial_df: pd.DataFrame, preferred: str = "spatial_region") -> Optional[str]:
    """自动寻找可用的 domain 列名，返回 None 表示未找到。"""
    if preferred in spatial_df.columns:
        return preferred
    for fallback in ["niche_semantic_label", "neighborhood_cluster"]:
        if fallback in spatial_df.columns:
            return fallback
    return None


# ============================================================
# Figure 1A: 课题技术路线图
# ============================================================

def plot_fig1A_workflow(out_dir: Path, dpi: int = PAPER_DPI) -> None:
    """绘制三列布局的课题技术路线图。

    ★ 可调参数：
      figsize=(15, 9)     — 图幅大小
      c_scrna / c_spatial / c_tcga — 三列主色
      _box: w=3.0/3.6, h=0.85 — 方块宽高
      fontsize — 各处字号
    """
    fig, ax = plt.subplots(figsize=(15, 9))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 9)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── 三列主色 ────────────────────────────────────────────
    c_scrna   = "#007a8b"
    c_spatial = "#f93800"
    c_tcga    = "#ffb500"
    c_bg      = "#f8f8f8"

    # ── 背景框 ────────────────────────────────────────────────────────
    for x0, w, color, title in [
        (0.2,  4.1, c_scrna,   "① scRNA-seq Reference"),
        (5.0,  4.8, c_spatial, "② Spatial Transcriptomics"),
        (10.5, 4.3, c_tcga,    "③ Clinical Validation"),
    ]:
        rect = mpatches.FancyBboxPatch(
            (x0, 0.3), w, 8.3,
            boxstyle="round,pad=0.2",
            facecolor=c_bg, edgecolor=color, linewidth=2.8, zorder=1,
        )
        ax.add_patch(rect)
        ax.text(x0 + w / 2, 8.78, title, ha="center", va="bottom",
                fontsize=11.5, fontweight="bold", color=color)

    # ── 内部方块与箭头辅助函数 ─────────────────────────────────────────
    def _box(ax_, x, y, text, fc, w=3.0, h=0.85):
        rect = mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.12",
            facecolor=fc, edgecolor="white", linewidth=0, alpha=0.88, zorder=2,
        )
        ax_.add_patch(rect)
        ax_.text(x, y, text, ha="center", va="center", fontsize=8.5,
                 color="white", fontweight="bold", zorder=3, wrap=True,
                 multialignment="center")

    def _arrow(ax_, x, y1, y2, color):
        ax_.annotate("", xy=(x, y2), xytext=(x, y1),
                     arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6), zorder=3)

    # ── scRNA-seq 左列 ──────────────────────────────────────────────
    scrna_x = 2.25
    for y, text in [
        (7.2, "Liver Cancer scRNA-seq\n(Multi-patient Atlas)"),
        (5.7, "Cell Clustering &\nType Annotation"),
        (4.2, "Treg Subcluster\nIdentification (FOXP3⁺)"),
        (2.7, "Cell-type Reference\nSignature Matrix"),
    ]:
        _box(ax, scrna_x, y, text, c_scrna)
    for y1, y2 in [(6.78, 6.15), (5.28, 4.65), (3.78, 3.15)]:
        _arrow(ax, scrna_x, y1, y2, c_scrna)

    # ── 空间转录组中列 ─────────────────────────────────────────────
    sp_x = 7.4
    sp_items = [
        (7.2, "HCC4R Visium Slide\n(Discovery Cohort)", "#f93800"),
        (5.7, "Cell2location\nDeconvolution", "#888"),
        (4.2, "Spatial Niche Discovery\n(kNN + Leiden)", "#888"),
        (2.7, "CHC20 Visium Slide\n(Validation, Label Transfer)", "#ffb500"),
    ]
    for y, text, color in sp_items:
        _box(ax, sp_x, y, text, color)
    for y1, y2 in [(6.78, 6.15), (5.28, 4.65), (3.78, 3.15)]:
        _arrow(ax, sp_x, y1, y2, c_spatial)

    # ── 临床队列右列 ───────────────────────────────────────────────
    tc_x = 12.65
    for y, text in [
        (7.2, "TCGA-LIHC / ICGC\n(n > 300 patients)"),
        (5.7, "ssGSEA Scoring\n(Niche Signature Genes)"),
        (4.2, "Survival Stratification\n(High vs Low)"),
        (2.7, "Kaplan-Meier Curve\n& Cox Regression"),
    ]:
        _box(ax, tc_x, y, text, c_tcga, w=3.6)
    for y1, y2 in [(6.78, 6.15), (5.28, 4.65), (3.78, 3.15)]:
        _arrow(ax, tc_x, y1, y2, c_tcga)

    # ── 跨列箭头 ──────────────────────────────────────────────────
    # scRNA → 空间（参考签名）
    ax.annotate("", xy=(5.75, 2.7), xytext=(3.75, 2.7),
                arrowprops=dict(arrowstyle="-|>", color="#666", lw=2.0), zorder=3)
    ax.text(4.75, 2.95, "cell_state_df", ha="center", fontsize=7.5,
            color="#666", style="italic")

    # 空间 → TCGA（签名基因）
    ax.annotate("", xy=(10.7, 4.2), xytext=(9.05, 4.2),
                arrowprops=dict(arrowstyle="-|>", color="#666", lw=2.0), zorder=3)
    ax.text(9.875, 4.45, "Signature Genes", ha="center", fontsize=7.5,
            color="#666", style="italic")

    ax.set_title(
        "Study Design: Spatial Immunosuppressive Niche in Hepatocellular Carcinoma",
        fontsize=13, fontweight="bold", pad=14, color="#333",
    )

    _save(fig, out_dir / "fig1A_workflow_diagram.png", dpi=dpi)


# ============================================================
# Figure 2A: H&E 染色图
# ============================================================

def plot_fig2A_HE(
    visium_dir: Optional[Path],
    spatial_df: pd.DataFrame,
    out_dir: Path,
    dpi: int = PAPER_DPI,
) -> None:
    """展示 HCC4R H&E 图像，若无原始图则用 spot 空间分布图代替。

    ★ 可调参数：
      figsize=(8, 8)  — 图幅大小
      s=2 / s=9       — 覆盖 spot / placeholder spot 大小
      ax.set_facecolor("#f5e6d3") — placeholder 背景色
    """
    he_image = None

    if visium_dir and visium_dir.exists():
        for fname in ["tissue_hires_image.png", "tissue_lowres_image.png",
                      "tissue_hires_image.jpg"]:
            p = visium_dir / "spatial" / fname
            if p.exists():
                try:
                    he_image = plt.imread(str(p))
                    logging.info("H&E image loaded: %s", p)
                    break
                except Exception as e:
                    logging.warning("Failed to load H&E image: %s", e)

    fig, ax = plt.subplots(figsize=(8, 8))

    if he_image is not None:
        ax.imshow(he_image, origin="upper")
        if "spatial_x" in spatial_df.columns:
            ax.scatter(spatial_df["spatial_x"], spatial_df["spatial_y"],
                       s=2, alpha=0.2, c="white", linewidths=0)
    else:
        logging.warning("H&E image unavailable. Using spatial region plot as placeholder.")
        ax.set_facecolor("#f5e6d3")
        if "spatial_x" in spatial_df.columns and "spatial_region" in spatial_df.columns:
            for region, color in DOMAIN_COLORS.items():
                mask = spatial_df["spatial_region"] == region
                if mask.sum() > 0:
                    ax.scatter(spatial_df.loc[mask, "spatial_x"],
                               spatial_df.loc[mask, "spatial_y"],
                               c=color, s=SPOT_SIZE, alpha=0.8, linewidths=0, label=region)
            ax.legend(loc="upper right", fontsize=9, frameon=False,
                      title="Spatial Region", title_fontsize=10)
        elif "spatial_x" in spatial_df.columns:
            ax.scatter(spatial_df["spatial_x"], spatial_df["spatial_y"],
                       c="#c49c94", s=SPOT_SIZE, alpha=0.7, linewidths=0)
        ax.text(0.5, 0.02,
                "Placeholder: set --hcc4r-visium-dir to display actual H&E image",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
                color="#888888", style="italic")

    ax.set_title("HCC4R: H&E Histology (Spatial Slide Overview)", fontsize=13, fontweight="bold")
    ax.axis("off")
    _save(fig, out_dir / "fig2A_HE_image.png", dpi=dpi)


# ============================================================
# Figure 2B: 核心细胞类型空间丰度图（每种细胞类型单独一张）
# ============================================================

def plot_fig2B_spatial_abundance(
    spatial_df: pd.DataFrame,
    cell_types: list[str],
    out_dir: Path,
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """为每种细胞类型单独生成空间丰度散点图（paper_ybp 配色）。

    ★ 可调参数：
      figsize=(7, 7)                  — 单张图幅
      s=9, alpha=0.9                  — spot 大小与透明度
      cmap=PAPER_YBP                  — 颜色映射
      vmin/vmax: percentile(2/98)     — 颜色映射范围（裁剪极端值）
      cbar shrink=0.65                — 色条缩放比例
    """
    for ct in cell_types:
        col_name = ct if ct in spatial_df.columns else ct.replace("/", "_")
        if col_name not in spatial_df.columns:
            logging.warning("Cell type '%s' not found in spatial_df; skipping.", ct)
            continue

        x = spatial_df["spatial_x"].to_numpy()
        y = spatial_df["spatial_y"].to_numpy()
        v = spatial_df[col_name].fillna(0).to_numpy()

        fig, ax = plt.subplots(figsize=(7, 7))
        sc = ax.scatter(
            x, y, c=v, s=SPOT_SIZE, cmap=PAPER_YBP,
            vmin=np.percentile(v, 2), vmax=np.percentile(v, 98),
            linewidths=0, alpha=0.9,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, shrink=0.65)
        cbar.set_label("Proportion", fontsize=10)
        cbar.ax.tick_params(labelsize=8)

        safe_name = ct.replace("/", "_").replace(" ", "_")
        ax.set_title(f"{ct} Spatial Abundance ({sample_name})", fontsize=13, fontweight="bold")
        ax.set_aspect("equal")
        ax.axis("off")

        _save(fig, out_dir / f"fig2B_spatial_{safe_name}.png", dpi=dpi)


# ============================================================
# Figure 2C: 细胞空间共定位相关性热图
# ============================================================

def plot_fig2C_coloc_heatmap(
    proportions: pd.DataFrame,
    cell_types: list[str],
    out_dir: Path,
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """计算并可视化细胞类型丰度间的 Pearson 相关性热图。

    ★ 可调参数：
      cmap="vlag"                 — 色板（center=0 的发散色板）
      annot_kws={"size": 8}       — 注释字号
      linewidths=0.4              — 单元格线宽
      fig_size: max(7, n*0.9+2)  — 图幅自适应列数
    """
    cols = [c for c in cell_types if c in proportions.columns]
    if len(cols) < 2:
        logging.warning("Not enough cell types for correlation heatmap; skipping 2C.")
        return

    corr = proportions[cols].corr(method="pearson")

    # 推荐排序（与生物含义一致）
    pref = ["Hepatocyte", "Malignant", "Treg", "T/NK", "Myeloid", "Fibroblast",
            "HSC", "B cell", "Endothelial", "Plasma cell", "Epithelial"]
    ordered = [c for c in pref if c in corr.columns] + \
              [c for c in corr.columns if c not in pref]
    corr = corr.loc[ordered, ordered]

    n = len(ordered)
    fig_size = max(7, n * 0.9 + 2)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.85))

    sns.heatmap(
        corr, ax=ax,
        cmap="vlag", center=0, vmin=-1, vmax=1,
        annot=True, fmt=".2f", annot_kws={"size": 8},
        linewidths=0.4, square=True,
        cbar_kws={"label": "Pearson r", "shrink": 0.75},
    )
    ax.set_title(
        f"Cell Type Spatial Co-localization ({sample_name})",
        fontsize=12, fontweight="bold", pad=10,
    )
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)

    _save(fig, out_dir / "fig2C_coloc_heatmap.png", dpi=dpi)


# ============================================================
# Figure 3A: 空间生态位聚类图
# ============================================================

def plot_fig3A_spatial_domain(
    spatial_df: pd.DataFrame,
    out_dir: Path,
    domain_col: str = "spatial_region",
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """以颜色展示各 Spatial Domain 的空间分布。

    ★ 可调参数：
      figsize=(8, 8)       — 图幅大小
      s=10, alpha=0.88     — spot 大小与透明度
      label fontsize=8.5   — 域标签字号
      legend markerscale=1.5 — 图例 marker 缩放
    """
    domain_col = _resolve_domain_col(spatial_df, domain_col)
    if domain_col is None:
        logging.warning("No domain column found; skipping 3A.")
        return

    x = spatial_df["spatial_x"].to_numpy()
    y = spatial_df["spatial_y"].to_numpy()
    labels = spatial_df[domain_col].astype(str).to_numpy()
    unique_labels = sorted(set(labels))

    color_map: dict[str, str] = {}
    for lbl in unique_labels:
        if lbl in DOMAIN_COLORS:
            color_map[lbl] = DOMAIN_COLORS[lbl]
        else:
            idx = len(color_map) % len(PAPER_CLUSTER_COLORS)
            color_map[lbl] = PAPER_CLUSTER_COLORS[idx]

    fig, ax = plt.subplots(figsize=(8, 8))
    for lbl in unique_labels:
        mask = labels == lbl
        ax.scatter(x[mask], y[mask], c=color_map[lbl], s=10, alpha=0.88,
                   linewidths=0, label=lbl)

    display_names = {
        "tumor_core":    "Tumor Core",
        "tumor_edge":    "Tumor Edge",
        "stroma_immune": "Stroma/Immune",
    }
    for lbl in unique_labels:
        mask = labels == lbl
        if mask.sum() == 0:
            continue
        display = display_names.get(lbl, lbl)
        if display and display != "other":
            cx, cy = float(x[mask].mean()), float(y[mask].mean())
            ax.text(cx, cy, display, fontsize=8.5, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

    ax.legend(loc="upper right", fontsize=9, frameon=False,
              title="Spatial Domain", title_fontsize=10, markerscale=1.5)
    ax.set_title(f"Spatial Niche Map ({sample_name})", fontsize=13, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")

    _save(fig, out_dir / "fig3A_spatial_domain_map.png", dpi=dpi)


# ============================================================
# Figure 3B: 空间生态位细胞组成 100% 堆叠条形图
# ============================================================

def plot_fig3B_domain_composition(
    spatial_df: pd.DataFrame,
    cell_types: list[str],
    out_dir: Path,
    domain_col: str = "spatial_region",
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """绘制各 Spatial Domain 的细胞类型组成 100% 堆叠条形图。

    ★ 可调参数：
      figsize: max(7, n_domain*1.8+3) × 6.5 — 图幅宽度自适应
      width=0.65          — 条形宽度
      label threshold: 5% — 低于此比例不显示文字标注
      fontsize 11/7.5     — 轴标签/内嵌文字字号
    """
    domain_col = _resolve_domain_col(spatial_df, domain_col)
    if domain_col is None:
        logging.warning("No domain column; skipping 3B.")
        return

    valid_cols = [c for c in cell_types if c in spatial_df.columns]
    if not valid_cols:
        valid_cols = [c.replace("/", "_") for c in cell_types
                      if c.replace("/", "_") in spatial_df.columns]
    if not valid_cols:
        logging.warning("No cell type columns found; skipping 3B.")
        return

    grouped = spatial_df.groupby(domain_col)[valid_cols].mean()
    row_sums = grouped.sum(axis=1)
    grouped_pct = grouped.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0) * 100

    domain_order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    domain_order = [d for d in domain_order if d in grouped_pct.index]
    domain_order += [d for d in grouped_pct.index if d not in domain_order]
    grouped_pct = grouped_pct.loc[domain_order]

    ct_colors = [CELLTYPE_COLORS.get(ct, PAPER_CLUSTER_COLORS[i % len(PAPER_CLUSTER_COLORS)])
                 for i, ct in enumerate(valid_cols)]

    fig, ax = plt.subplots(figsize=(max(7, len(domain_order) * 1.8 + 3), 6.5))

    bottom = np.zeros(len(grouped_pct))
    for ct, color in zip(valid_cols, ct_colors):
        vals = grouped_pct[ct].to_numpy()
        bars = ax.bar(range(len(grouped_pct)), vals, bottom=bottom,
                      color=color, label=ct.replace("_", "/"), width=0.65,
                      edgecolor="white", linewidth=0.3)
        for i, (bar, val) in enumerate(zip(bars, vals)):
            if val > 5:
                ax.text(bar.get_x() + bar.get_width() / 2, bottom[i] + val / 2,
                        f"{val:.0f}%", ha="center", va="center",
                        fontsize=7.5, color="white", fontweight="bold")
        bottom += vals

    bar_display_names = {
        "tumor_core":    "Tumor\nCore",
        "tumor_edge":    "Tumor\nEdge",
        "stroma_immune": "Stroma /\nImmune",
        "other":         "Other",
    }
    ax.set_xticks(range(len(grouped_pct)))
    ax.set_xticklabels([bar_display_names.get(d, d) for d in grouped_pct.index],
                        fontsize=11, fontweight="bold")
    ax.set_ylabel("Cell Type Proportion (%)", fontsize=11)
    ax.set_ylim(0, 110)
    ax.set_title(f"Cell Type Composition per Spatial Domain ({sample_name})",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8.5, frameon=False,
              bbox_to_anchor=(1.2, 1.0), title="Cell Type", title_fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=9)

    _save(fig, out_dir / "fig3B_domain_composition_barplot.png", dpi=dpi)


# ============================================================
# Figure 3C: 多层次评分小提琴图（每个评分单独一张图）
# ============================================================

def plot_fig3C_score_violins(
    spatial_df: pd.DataFrame,
    out_dir: Path,
    domain_col: str = "spatial_region",
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """为三个评分分别生成小提琴图，每张图单独保存。

    ★ 可调参数：
      score_configs — 评分列名、显示名称、输出文件名、颜色
      VIOLIN_ALPHA=0.55          — 小提琴透明度（全局常量）
      figsize: max(7, n*1.7+2)  — 图幅自适应 domain 数
      widths=0.7 / boxplot widths=0.12 — 小提琴/箱线宽度
      significance: tumor_core vs tumor_edge — 显著性标注对象
    """
    # ★ 修改这里来调整各评分的颜色或新增评分
    score_configs = [
        ("Treg_like_score",               "Treg-like Score",
         "fig3C_violin_treg_like_score.png",               "#007a8b"),
        ("immune_stroma_score",           "Immune-Stroma Score",
         "fig3C_violin_immune_stroma_score.png",           "#f93800"),
        ("immunosuppressive_niche_score", "Immunosuppressive Niche Score",
         "fig3C_violin_immunosuppressive_niche_score.png", "#ffb500"),
    ]

    domain_col = _resolve_domain_col(spatial_df, domain_col)
    if domain_col is None:
        logging.warning("No domain column; skipping 3C.")
        return

    domain_order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    domain_order = [d for d in domain_order if d in spatial_df[domain_col].unique()]
    domain_order += sorted(d for d in spatial_df[domain_col].unique()
                           if d not in domain_order)

    violin_display = {
        "tumor_core":    "Tumor Core",
        "tumor_edge":    "Tumor Edge",
        "stroma_immune": "Stroma/\nImmune",
        "other":         "Other",
    }

    for score_col, score_label, fname, color in score_configs:
        if score_col not in spatial_df.columns:
            logging.warning("Score '%s' not found; skipping.", score_col)
            continue

        groups = [spatial_df.loc[spatial_df[domain_col] == d, score_col].dropna().to_numpy()
                  for d in domain_order]

        fig, ax = plt.subplots(figsize=(max(7, len(domain_order) * 1.7 + 2), 5.5))

        vp = ax.violinplot(groups, positions=range(len(domain_order)),
                           widths=0.7, showmedians=True, showextrema=True)
        for pc in vp["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(VIOLIN_ALPHA)
            pc.set_edgecolor("white")
        vp["cmedians"].set_color("#222222")
        vp["cmedians"].set_linewidth(2.2)
        for part_key in ["cbars", "cmins", "cmaxes"]:
            if part_key in vp:
                vp[part_key].set_color("#888888")
                vp[part_key].set_linewidth(1)

        # 内嵌箱线
        ax.boxplot(groups, positions=range(len(domain_order)),
                   widths=0.12, patch_artist=True,
                   medianprops=dict(color="black", linewidth=2),
                   boxprops=dict(facecolor="white", alpha=0.9, linewidth=1),
                   whiskerprops=dict(linewidth=0), capprops=dict(linewidth=0),
                   flierprops=dict(marker=""))

        # 显著性标注（tumor_core vs tumor_edge）
        if "tumor_core" in domain_order and "tumor_edge" in domain_order:
            i0 = domain_order.index("tumor_core")
            i1 = domain_order.index("tumor_edge")
            g0, g1 = groups[i0], groups[i1]
            if len(g0) > 0 and len(g1) > 0:
                p_val = _wilcoxon_test(g1, g0)
                all_vals = np.concatenate([g for g in groups if len(g) > 0])
                y_max = np.percentile(all_vals, 99)
                y_span = np.percentile(all_vals, 99) - np.percentile(all_vals, 1)
                bh = y_max + y_span * 0.08
                ax.plot([i0, i0, i1, i1],
                        [bh, bh + y_span * 0.03, bh + y_span * 0.03, bh],
                        c="#333", lw=1)
                ax.text((i0 + i1) / 2, bh + y_span * 0.04,
                        _p_label(p_val), ha="center", va="bottom", fontsize=10)

        ax.set_xticks(range(len(domain_order)))
        ax.set_xticklabels([violin_display.get(d, d) for d in domain_order],
                            fontsize=10.5, fontweight="bold")
        ax.set_ylabel(score_label, fontsize=11)
        ax.set_title(f"{score_label} by Spatial Domain ({sample_name})",
                     fontsize=13, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=9)

        _save(fig, out_dir / fname, dpi=dpi)


# ============================================================
# Figure 4A: 差异基因火山图
# ============================================================

def plot_fig4A_volcano(
    niche_scores_df: pd.DataFrame,
    adata: ad.AnnData,
    out_dir: Path,
    sample_name: str = "HCC4R",
    fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    dpi: int = PAPER_DPI,
) -> None:
    """绘制 niche_high vs niche_low 的全基因火山图（论文版）。

    ★ 可调参数：
      highlight_genes    — 需要标注的基因列表
      fc_threshold=0.5   — log2FC 阈值（竖虚线位置）
      fdr_threshold=0.05 — FDR 阈值（横虚线位置）
      figsize=(9, 7)     — 图幅大小
      s=10, alpha=0.65   — spot 大小与透明度
      c_up="#d62728", c_dn="#4575b4", c_ns="#bdbdbd" — 三类点颜色
    """
    # ★ 修改这里来调整需要标注的基因
    highlight_genes = [
        "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
        "TGFB1", "IL10", "CXCL12", "CCL22",
        "FAP", "ACTA2", "POSTN", "COL1A1",
        "CD163", "MRC1", "SPP1", "VEGFA",
        "LAG3", "PDCD1", "HAVCR2",
    ]
    # 颜色
    c_up = "#d62728"
    c_dn = "#4575b4"
    c_ns = "#bdbdbd"

    if "niche_high" not in niche_scores_df.columns:
        logging.warning("'niche_high' not found; skipping 4A.")
        return

    logging.info("Figure 4A: Computing differential expression...")
    expr_df = _expression_frame(adata)

    if "spot_id" in niche_scores_df.columns:
        nh = niche_scores_df.set_index("spot_id")["niche_high"]
    else:
        nh = niche_scores_df["niche_high"]
    niche_high = nh.reindex(expr_df.index).fillna(False).astype(bool)

    high_expr = expr_df.loc[niche_high]
    low_expr  = expr_df.loc[~niche_high]
    mean_h = high_expr.mean()
    mean_l = low_expr.mean()
    lfc    = np.log2((mean_h + 1) / (mean_l + 1))
    dfrac  = (high_expr > 0).mean() - (low_expr > 0).mean()

    pvals = np.ones(len(expr_df.columns))
    cand_mask = (lfc.to_numpy() > 0.1) | (dfrac.to_numpy() > 0.03)
    gene_list = list(expr_df.columns)
    for gi in np.where(cand_mask)[0]:
        gene = gene_list[gi]
        try:
            _, p = scipy_stats.mannwhitneyu(
                high_expr[gene].to_numpy(), low_expr[gene].to_numpy(),
                alternative="greater",
            )
            pvals[gi] = p
        except Exception:
            pass

    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    neg_log10_fdr = -np.log10(np.clip(fdr, 1e-300, 1))

    c_list = [
        c_up if (lfc.iloc[i] > fc_threshold and fdr[i] < fdr_threshold)
        else c_dn if (lfc.iloc[i] < -fc_threshold and fdr[i] < fdr_threshold)
        else c_ns
        for i in range(len(lfc))
    ]

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(lfc.to_numpy(), neg_log10_fdr, c=c_list, s=10, alpha=0.65, linewidths=0)
    ax.axvline(fc_threshold,  c="black", ls="--", lw=0.9)
    ax.axvline(-fc_threshold, c="black", ls="--", lw=0.9)
    ax.axhline(-np.log10(fdr_threshold), c="black", ls="--", lw=0.9)

    labeled = set()
    for i, gene in enumerate(gene_list):
        if gene not in highlight_genes or gene in labeled:
            continue
        lfc_v, fdr_v = float(lfc.iloc[i]), fdr[i]
        if lfc_v > fc_threshold * 0.5 and fdr_v < fdr_threshold * 5:
            ax.annotate(gene, (lfc_v, neg_log10_fdr[i]),
                        textcoords="offset points", xytext=(5, 2),
                        fontsize=7.5, color="#222",
                        arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.5))
            labeled.add(gene)

    n_up = int(((lfc > fc_threshold) & (fdr < fdr_threshold)).sum())
    n_dn = int(((lfc < -fc_threshold) & (fdr < fdr_threshold)).sum())
    legend_els = [
        mpatches.Patch(facecolor=c_up, label=f"Up-regulated (n={n_up})"),
        mpatches.Patch(facecolor=c_dn, label=f"Down-regulated (n={n_dn})"),
        mpatches.Patch(facecolor=c_ns, label="Not significant"),
    ]
    ax.legend(handles=legend_els, fontsize=9, frameon=False, loc="upper left")
    ax.set_xlabel("log₂FC (Niche-High / Niche-Low)", fontsize=11)
    ax.set_ylabel("−log₁₀(FDR)", fontsize=11)
    ax.set_title(
        f"Spatial Differential Expression: Immunosuppressive Niche ({sample_name})",
        fontsize=12, fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    _save(fig, out_dir / "fig4A_volcano.png", dpi=dpi)


# ============================================================
# Figure 4B: Niche Marker 热图（签名基因 × Spatial Domain）
# ============================================================

def plot_fig4B_niche_marker_heatmap(
    adata: ad.AnnData,
    spatial_df: pd.DataFrame,
    ranked_genes_path: Optional[Path],
    out_dir: Path,
    top_n: int = 30,
    domain_col: str = "spatial_region",
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """展示 Top-N 签名基因在各 Spatial Domain 中的平均表达热图。

    ★ 可调参数：
      top_n=30           — 展示的签名基因数量
      cmap=PAPER_YBP     — 热图颜色映射
      vmin=-2, vmax=2    — Z-score 颜色范围
      fig_size: max(8, n_d*1.6+2) × max(8, n_g*0.32+2) — 自适应图幅
      fallback sig_genes — 当无签名基因文件时使用的预设基因列表
    """
    from scipy.stats import zscore as scipy_zscore

    if ranked_genes_path and ranked_genes_path.exists():
        ranked = pd.read_csv(ranked_genes_path)
        sort_col = "composite_score" if "composite_score" in ranked.columns else "log2_fc"
        sig_genes = (ranked.sort_values(sort_col, ascending=False)
                     .head(top_n)["gene"].tolist())
    else:
        # ★ 修改这里来更换默认备用基因列表
        sig_genes = [
            "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
            "TGFB1", "IL10", "CXCL12", "CCL22",
            "FAP", "ACTA2", "POSTN", "COL1A1",
            "CD163", "MRC1", "SPP1",
            "ALB", "APOA2", "APOA1",
        ]
        logging.warning("Signature gene file not found; using preset list.")

    available = [g for g in sig_genes if g in adata.var_names]
    if not available:
        logging.warning("No signature genes in adata; skipping 4B.")
        return

    expr_df = _expression_frame(adata)[available]

    domain_col = _resolve_domain_col(spatial_df, domain_col)
    if domain_col is None:
        logging.warning("Domain column not found; skipping 4B.")
        return

    if "spot_id" in spatial_df.columns:
        domain_map = spatial_df.set_index("spot_id")[domain_col]
    else:
        domain_map = spatial_df[domain_col]

    expr_df = expr_df.copy()
    expr_df["domain"] = domain_map.reindex(expr_df.index).fillna("other")

    mean_by_domain = expr_df.groupby("domain")[available].mean()
    domain_order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    domain_order = [d for d in domain_order if d in mean_by_domain.index]
    domain_order += [d for d in mean_by_domain.index if d not in domain_order]
    mean_by_domain = mean_by_domain.loc[domain_order]

    mean_z = mean_by_domain.apply(
        lambda col: scipy_zscore(col) if col.std() > 0 else col, axis=0
    ).fillna(0.0).T  # rows=genes, cols=domains

    n_g, n_d = mean_z.shape
    fig, ax = plt.subplots(figsize=(max(6, n_d * 1.6 + 2), max(8, n_g * 0.32 + 2)))
    im = ax.imshow(mean_z.values, aspect="auto", cmap=PAPER_YBP, vmin=-2, vmax=2)

    domain_display = {
        "tumor_core": "Tumor Core", "tumor_edge": "Tumor Edge",
        "stroma_immune": "Stroma/\nImmune", "other": "Other",
    }
    ax.set_xticks(range(n_d))
    ax.set_xticklabels([domain_display.get(d, d) for d in domain_order],
                        rotation=30, ha="right", fontsize=10, fontweight="bold")
    ax.set_yticks(range(n_g))
    ax.set_yticklabels(mean_z.index.tolist(), fontsize=8.5, fontstyle="italic")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04, shrink=0.6)
    cbar.set_label("Expression (Z-score)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(f"Niche Signature Genes × Spatial Domain ({sample_name}, Top-{top_n})",
                 fontsize=12, fontweight="bold")

    _save(fig, out_dir / "fig4B_niche_marker_heatmap.png", dpi=dpi)


# ============================================================
# Figure 4C: 通路富集气泡图（简化 ssGSEA）
# ============================================================

def plot_fig4C_pathway_bubble(
    adata: ad.AnnData,
    spatial_df: pd.DataFrame,
    out_dir: Path,
    domain_col: str = "spatial_region",
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """对各 Spatial Domain 进行预定义通路的基因集评分，以气泡图展示。

    X 轴：Spatial Domain
    Y 轴：通路名称
    气泡大小：富集评分（越大越富集）
    气泡颜色：红=高于全局均值，蓝=低于全局均值

    ★ 可调参数：
      PATHWAYS          — 通路-基因字典（可增删通路或修改基因集）
      bubble_size: 30 + 350 * ... — 气泡大小缩放公式
      figsize: max(8, n_d*1.8+2) × max(6, n_p*0.7+2) — 自适应图幅
      c_up="#d62728", c_dn="#4575b4" — 富集方向颜色
    """
    # ★ 修改这里来增删通路或调整基因集
    PATHWAYS = {
        "TGF-β Signaling":             ["TGFB1", "TGFB2", "TGFBR1", "SMAD2", "SMAD3", "SMAD4"],
        "EMT":                          ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ACTA2"],
        "Angiogenesis":                 ["VEGFA", "VEGFB", "KDR", "FLT1", "ANGPT2", "PECAM1"],
        "IL6-JAK-STAT3":               ["IL6", "JAK1", "JAK2", "STAT3", "SOCS3", "BCL2"],
        "Treg Immune Suppression":      ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IL10", "TGFB1"],
        "TAM Polarization (M2)":        ["CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CCL22"],
        "CAF Activation":               ["FAP", "ACTA2", "POSTN", "COL1A1", "COL1A2", "FN1"],
        "Immune Checkpoint":            ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
        "Hypoxia":                      ["HIF1A", "VEGFA", "LDHA", "PGK1", "SLC2A1", "BNIP3"],
        "TNF-α Signaling via NF-κB":   ["TNFAIP3", "ICAM1", "IL6", "MMP9", "BIRC3", "NFKBIA"],
    }

    domain_col = _resolve_domain_col(spatial_df, domain_col)
    if domain_col is None:
        logging.warning("No domain column; skipping 4C.")
        return

    expr_df = _expression_frame(adata)

    if "spot_id" in spatial_df.columns:
        domain_map = spatial_df.set_index("spot_id")[domain_col]
    else:
        domain_map = spatial_df[domain_col]

    expr_df = expr_df.copy()
    expr_df["domain"] = domain_map.reindex(expr_df.index).fillna("other")

    domain_order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    domain_order = [d for d in domain_order if d in expr_df["domain"].unique()]
    domain_order += [d for d in expr_df["domain"].unique() if d not in domain_order]

    records = []
    for pathway, genes in PATHWAYS.items():
        avail_genes = [g for g in genes if g in expr_df.columns]
        if not avail_genes:
            continue
        scores = expr_df[avail_genes].mean(axis=1)
        expr_df_temp = expr_df.copy()
        expr_df_temp["_score"] = scores
        mean_scores = expr_df_temp.groupby("domain")["_score"].mean()
        global_mean = float(mean_scores.mean())
        for domain in domain_order:
            if domain not in mean_scores:
                continue
            s = float(mean_scores[domain])
            records.append({
                "pathway": pathway,
                "domain": domain,
                "score": s,
                "delta": s - global_mean,
                "n_genes": len(avail_genes),
            })

    if not records:
        logging.warning("No pathway scores computed; skipping 4C.")
        return

    result_df = pd.DataFrame(records)

    # 排序通路（按 tumor_edge 评分降序）
    if "tumor_edge" in domain_order:
        edge_scores = (result_df[result_df["domain"] == "tumor_edge"]
                       .set_index("pathway")["score"])
        pathway_order = edge_scores.sort_values(ascending=False).index.tolist()
    else:
        pathway_order = list(PATHWAYS.keys())
    pathway_order = [p for p in pathway_order if p in result_df["pathway"].unique()]

    n_pathways = len(pathway_order)
    n_domains  = len(domain_order)

    fig, ax = plt.subplots(figsize=(max(8, n_domains * 1.8 + 2),
                                    max(6, n_pathways * 0.7 + 2)))

    all_scores = result_df["score"].to_numpy()
    s_min, s_max = all_scores.min(), all_scores.max()
    s_range = s_max - s_min if s_max > s_min else 1.0

    bubble_display = {
        "tumor_core": "Tumor\nCore",
        "tumor_edge": "Tumor\nEdge",
        "stroma_immune": "Stroma/\nImmune",
        "other": "Other",
    }

    for row_i, pathway in enumerate(pathway_order):
        for col_i, domain in enumerate(domain_order):
            sub = result_df[(result_df["pathway"] == pathway) &
                            (result_df["domain"] == domain)]
            if sub.empty:
                continue
            score = float(sub["score"].iloc[0])
            delta = float(sub["delta"].iloc[0])
            bubble_size = 30 + 350 * (score - s_min) / s_range
            color = "#d62728" if delta > 0 else "#4575b4"
            ax.scatter(col_i, row_i, s=bubble_size, c=color, alpha=0.8,
                       linewidths=0.3, edgecolors="white")

    ax.set_xticks(range(n_domains))
    ax.set_xticklabels([bubble_display.get(d, d) for d in domain_order],
                        fontsize=10.5, fontweight="bold")
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(pathway_order, fontsize=9)
    ax.set_xlim(-0.5, n_domains - 0.5)
    ax.set_ylim(-0.5, n_pathways - 0.5)
    ax.invert_yaxis()

    legend_els = [
        mpatches.Patch(facecolor="#d62728", alpha=0.8, label="Above global mean"),
        mpatches.Patch(facecolor="#4575b4", alpha=0.8, label="Below global mean"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
               markersize=6,  label="Low score"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888",
               markersize=11, label="High score"),
    ]
    ax.legend(handles=legend_els, loc="upper right", fontsize=8,
              frameon=False, bbox_to_anchor=(1.28, 1.0))
    ax.set_title(f"Pathway Enrichment by Spatial Domain ({sample_name})\n"
                 "(bubble size = mean expression score; color = relative to global mean)",
                 fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)
    ax.set_facecolor("#f9f9f9")
    for y_grid in range(n_pathways):
        ax.axhline(y_grid, color="white", lw=0.8, zorder=0)
    for x_grid in range(n_domains):
        ax.axvline(x_grid, color="white", lw=0.8, zorder=0)

    _save(fig, out_dir / "fig4C_pathway_bubble.png", dpi=dpi)


# ============================================================
# Figure 5A: HCC4R-derived signature projection in CHC20
# ============================================================

def plot_fig5A_signature_projection(
    chc20_adata: ad.AnnData,
    signature_score: pd.Series,
    out_dir: Path,
    signature_name: str = "HCC4R-derived niche signature",
    matched_genes: Optional[list[str]] = None,
    dpi: int = PAPER_DPI,
) -> None:
    """Project an HCC4R-derived niche signature score onto CHC20 coordinates."""
    if "spatial" not in chc20_adata.obsm:
        logging.warning("CHC20 AnnData has no obsm['spatial']; skipping Figure 5A.")
        return

    score = signature_score.reindex(chc20_adata.obs_names).dropna()
    if score.empty:
        logging.warning("CHC20 signature score is empty; skipping Figure 5A.")
        return

    coords = np.asarray(chc20_adata.obsm["spatial"])
    coord_df = pd.DataFrame(coords, index=chc20_adata.obs_names, columns=["x", "y"])
    coord_df = coord_df.loc[score.index]
    values = score.to_numpy(dtype=float)

    if len(values) > 2:
        vmin, vmax = np.percentile(values, [2, 98])
        if np.isclose(vmin, vmax):
            vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
    else:
        vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    sc = ax.scatter(
        coord_df["x"], coord_df["y"],
        c=values, s=SPOT_SIZE, cmap=PAPER_YBP,
        vmin=vmin, vmax=vmax, linewidths=0, alpha=0.92,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, shrink=0.68)
    cbar.set_label("Signature score", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    gene_note = f"n={len(matched_genes)} genes" if matched_genes else "matched genes"
    ax.set_title(
        f"CHC20: {signature_name}\n"
        f"Projected from HCC4R ({gene_note})",
        fontsize=12, fontweight="bold",
    )
    ax.set_aspect("equal")
    ax.axis("off")

    _save(fig, out_dir / "fig5A_hcc4r_signature_projection_CHC20.png", dpi=dpi)


# ============================================================
# Figure 5A: CHC20 域标签迁移预测图（随机森林分类器）
# ============================================================

def plot_fig5A_label_transfer(
    hcc4r_niche_df: pd.DataFrame,
    hcc4r_adata: ad.AnnData,
    chc20_adata: ad.AnnData,
    out_dir: Path,
    domain_col: str = "spatial_region",
    dpi: int = PAPER_DPI,
) -> None:
    """用 HCC4R 训练随机森林分类器，在 CHC20 上预测 Spatial Domain 标签并可视化。

    特征 = 各 spot 的细胞类型比例（来自 cell2location 丰度矩阵）。

    ★ 可调参数：
      n_estimators=200    — 随机森林树数量
      class_weight="balanced" — 类别权重（解决类别不平衡）
      figsize=(8, 8)      — 图幅大小
      s=10, alpha=0.85    — spot 大小与透明度
      label fontsize=8.5  — 域标签字号
    """
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.preprocessing import LabelEncoder
    except ImportError:
        logging.warning("sklearn not available; skipping Figure 5A.")
        return

    if domain_col not in hcc4r_niche_df.columns:
        logging.warning("Domain col '%s' not in HCC4R niche df; skipping 5A.", domain_col)
        return

    try:
        hcc4r_abund = _get_abundance(hcc4r_adata)
        hcc4r_prop  = _proportions(hcc4r_abund)
    except Exception as e:
        logging.warning("Failed to get HCC4R abundance: %s; skipping 5A.", e)
        return

    try:
        chc20_abund = _get_abundance(chc20_adata)
        chc20_prop  = _proportions(chc20_abund)
    except Exception as e:
        logging.warning("Failed to get CHC20 abundance: %s; skipping 5A.", e)
        return

    common_ct = [c for c in hcc4r_prop.columns if c in chc20_prop.columns]
    if not common_ct:
        logging.warning("No common cell types between HCC4R and CHC20; skipping 5A.")
        return

    if "spot_id" in hcc4r_niche_df.columns:
        label_map = hcc4r_niche_df.set_index("spot_id")[domain_col]
    else:
        label_map = hcc4r_niche_df[domain_col]

    y_train = label_map.reindex(hcc4r_prop.index).fillna("other")
    X_train = hcc4r_prop[common_ct].fillna(0)

    le = LabelEncoder()
    y_enc = le.fit_transform(y_train.to_numpy())
    clf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1,
                                 class_weight="balanced")
    clf.fit(X_train.to_numpy(), y_enc)

    X_test = chc20_prop[common_ct].fillna(0)
    y_pred_enc = clf.predict(X_test.to_numpy())
    y_pred = le.inverse_transform(y_pred_enc)

    chc20_coords = np.asarray(chc20_adata.obsm["spatial"])
    x = chc20_coords[:, 0]
    y = chc20_coords[:, 1]

    unique_labels = sorted(set(y_pred))
    color_map = {
        lbl: DOMAIN_COLORS.get(lbl, PAPER_CLUSTER_COLORS[i % len(PAPER_CLUSTER_COLORS)])
        for i, lbl in enumerate(unique_labels)
    }

    fig, ax = plt.subplots(figsize=(8, 8))
    for lbl in unique_labels:
        mask = y_pred == lbl
        ax.scatter(x[mask], y[mask], c=color_map[lbl], s=10, alpha=0.85,
                   linewidths=0, label=lbl)

    label_display_5a = {
        "tumor_core": "Tumor Core", "tumor_edge": "Tumor Edge",
        "stroma_immune": "Stroma/Immune",
    }
    for lbl in unique_labels:
        mask = y_pred == lbl
        if mask.sum() == 0:
            continue
        display = label_display_5a.get(lbl, lbl)
        if display and display != "other":
            cx, cy = float(x[mask].mean()), float(y[mask].mean())
            ax.text(cx, cy, display, fontsize=8.5, ha="center", va="center",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

    ax.legend(loc="upper right", fontsize=9, frameon=False,
              title="Predicted Domain", title_fontsize=10, markerscale=1.5)
    ax.set_title("CHC20: Spatial Domain Label Transfer\n"
                 "(Predicted by HCC4R-trained Random Forest)",
                 fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")

    _save(fig, out_dir / "fig5A_label_transfer.png", dpi=dpi)


# ============================================================
# Figure 5B: HCC4R vs CHC20 验证对比小提琴图
# ============================================================

def plot_fig5B_validation_violin(
    hcc4r_score_df: pd.DataFrame,
    chc20_score_df: pd.DataFrame,
    out_dir: Path,
    score_cols: Optional[list[str]] = None,
    dpi: int = PAPER_DPI,
) -> None:
    """Compare HCC4R and CHC20 distributions of the same signature score.

    Figure 5B now uses the HCC4R-derived signature score in both datasets,
    matching Figure 5A's projected CHC20 score.
    """
    if score_cols is None:
        score_cols = ["hcc4r_signature_score"]

    display_names = {
        "hcc4r_signature_score": "HCC4R-derived Niche Signature Score",
        "immunosuppressive_niche_score": "Immunosuppressive Niche Score",
    }

    for score in score_cols:
        hcc4r_vals = (
            hcc4r_score_df[score].dropna().to_numpy(dtype=float)
            if score in hcc4r_score_df.columns else np.array([])
        )
        chc20_vals = (
            chc20_score_df[score].dropna().to_numpy(dtype=float)
            if score in chc20_score_df.columns else np.array([])
        )

        if len(hcc4r_vals) == 0 or len(chc20_vals) == 0:
            logging.warning(
                "Score '%s' missing or empty in HCC4R/CHC20 signature scores; skipping 5B.",
                score,
            )
            continue

        if len(hcc4r_vals) > 1 and len(chc20_vals) > 1:
            _, p_two = scipy_stats.mannwhitneyu(
                hcc4r_vals, chc20_vals, alternative="two-sided"
            )
        else:
            p_two = 1.0

        fig, ax = plt.subplots(figsize=(7, 5.5))

        groups = [hcc4r_vals, chc20_vals]
        positions = [0, 1]
        colors_vl = ["#f4b41a", "#143d59"]
        labels = ["HCC4R\n(Discovery)", "CHC20\n(Validation)"]

        vp = ax.violinplot(
            groups, positions=positions, widths=0.6,
            showmedians=True, showextrema=True,
        )
        for i, pc in enumerate(vp["bodies"]):
            pc.set_facecolor(colors_vl[i])
            pc.set_alpha(VIOLIN_ALPHA)
            pc.set_edgecolor("white")
        vp["cmedians"].set_color("#222")
        vp["cmedians"].set_linewidth(2.2)

        ax.boxplot(
            groups, positions=positions, widths=0.1, patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            boxprops=dict(facecolor="white", alpha=0.9, linewidth=1),
            whiskerprops=dict(linewidth=0), capprops=dict(linewidth=0),
            flierprops=dict(marker=""),
        )

        all_vals = np.concatenate(groups)
        y_max = np.percentile(all_vals, 99)
        y_span = np.percentile(all_vals, 99) - np.percentile(all_vals, 1)
        if np.isclose(y_span, 0):
            y_span = max(abs(float(y_max)), 1.0)
        bh = y_max + y_span * 0.08
        ax.plot([0, 0, 1, 1],
                [bh, bh + y_span * 0.03, bh + y_span * 0.03, bh],
                c="#333", lw=1)
        p_str = f"p={p_two:.3f}" if p_two >= 0.001 else "p<0.001"
        ax.text(0.5, bh + y_span * 0.04, p_str,
                ha="center", va="bottom", fontsize=10, color="#333")

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=11, fontweight="bold")
        score_name = display_names.get(score, score)
        ax.set_ylabel(score_name, fontsize=11)
        ax.set_title(
            f"{score_name}\nHCC4R vs CHC20 Distribution",
            fontsize=12, fontweight="bold",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=9)

        safe = score.replace("/", "_")
        _save(fig, out_dir / f"fig5B_signature_score_violin_{safe}.png", dpi=dpi)


# ============================================================
# Figure 5C: 敏感性分析与参数扫描（两张独立图）
# ============================================================

def plot_fig5C_sensitivity(
    hcc4r_niche_dir: Path,
    out_dir: Path,
    dpi: int = PAPER_DPI,
) -> None:
    """直接使用 Step 2 已生成的敏感性分析结果，重新以论文风格绘制。"""
    sens_path = hcc4r_niche_dir / "sensitivity_analysis.csv"
    if sens_path.exists():
        sens_df = pd.read_csv(sens_path)
        _plot_sensitivity_paper(sens_df, out_dir / "fig5C_sensitivity_stability.png", dpi=dpi)
    else:
        logging.warning("sensitivity_analysis.csv not found; skipping Figure 5C.")


def _plot_sensitivity_paper(sens_df: pd.DataFrame, path: Path, dpi: int = PAPER_DPI) -> None:
    """论文版敏感性分析图：2×2 布局（kNN/radius × Spearman ρ/Weighted Jaccard）。

    ★ 可调参数：
      figsize=(11, 8)          — 图幅大小
      thresholds               — 各指标的合格阈值（橙色虚线）
      bar alpha=0.85           — 条形透明度
      colors: "#d62728"=参考, "#4575b4"=其他 — 条形颜色
    """
    req_cols = ["param_mode", "param_label", "spearman_rho", "weighted_jaccard", "is_reference"]
    if not all(c in sens_df.columns for c in req_cols):
        logging.warning("sensitivity_analysis.csv missing required columns; skipping.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle("Sensitivity Analysis: Niche Score Stability",
                 fontsize=13, fontweight="bold", y=1.01)

    metrics = ["spearman_rho", "weighted_jaccard"]
    modes   = ["knn", "radius"]
    y_labels = {"spearman_rho": "Spearman ρ", "weighted_jaccard": "Weighted Jaccard"}
    # ★ 修改这里来调整合格阈值线
    thresholds = {"spearman_rho": 0.90, "weighted_jaccard": 0.80}
    mode_labels = {"knn": "kNN Neighbors", "radius": "Radius Multiplier"}

    for row, metric in enumerate(metrics):
        for col, mode in enumerate(modes):
            ax = axes[row, col]
            sub = sens_df[sens_df["param_mode"] == mode].copy()
            if sub.empty:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes)
                continue

            colors = ["#d62728" if bool(r) else "#4575b4"
                      for r in sub["is_reference"]]
            ax.bar(range(len(sub)), sub[metric], color=colors, alpha=0.85,
                   edgecolor="white", linewidth=0.5)
            ax.axhline(thresholds[metric], color="#ff7f0e", linestyle="--",
                       linewidth=1.5, label=f"Threshold ({thresholds[metric]})")
            ax.axhline(1.0, color="#d62728", linestyle=":", linewidth=1.0,
                       label="Reference (=1.0)")
            ax.set_xticks(range(len(sub)))
            ax.set_xticklabels(sub["param_label"].tolist(), rotation=30, ha="right",
                                fontsize=8)
            ax.set_ylim(0, 1.12)
            ax.set_ylabel(y_labels[metric], fontsize=10)
            ax.set_title(f"{mode_labels[mode]} — {y_labels[metric]}", fontsize=10)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if row == 0 and col == 0:
                ax.legend(fontsize=8, frameon=False)

    plt.tight_layout()
    _save(fig, path, dpi=dpi)


