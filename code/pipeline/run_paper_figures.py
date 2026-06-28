"""
================================================================================
脚本名称: run_paper_figures.py
功能概述: 论文图表专用生成脚本
================================================================================

【整体任务说明】
    本脚本基于 run_preprocessing.py (Step 1) 和 run_spatial_niche_analysis.py
    (Step 2) 的输出结果，专门生成适合插入论文的高质量图表。

    分析架构：
      - HCC4R    → 主分析切片（Discovery Cohort）
      - CHC20    → 独立验证切片（Validation Cohort）
      - TCGA/ICGC → 外部临床队列（Translational Validation）

    所有子图单独生成（不拼接），DPI=300，适合期刊投稿。

    ── Figure 1: 课题总体设计与单细胞细胞图谱 ──
      1A: 课题技术路线图（Workflow Diagram）
      1B: scRNA-seq 细胞分群 tSNE/UMAP 图（使用 Step 1 已生成的图）
      1C: 细胞类型特异性 Marker 基因热图（使用 Step 1 已生成的图）
      1D: T/NK 亚群 Treg 鉴定气泡图（使用 Step 1 已生成的图）

    ── Figure 2: 空间解卷积与多细胞成分联合映射 ──
      2A: HCC4R H&E 染色图（fig2A_HE_image.png）
      2B: 核心细胞类型空间丰度图（4 张：Malignant/HSC/Treg/Myeloid）
      2C: 细胞空间共定位相关性热图（fig2C_coloc_heatmap.png）

    ── Figure 3: 主分析切片 HCC4R 空间生态位聚类与功能定量 ──
      3A: 空间生态位聚类图（fig3A_spatial_domain_map.png）
      3B: 空间生态位细胞组成堆叠条形图（fig3B_domain_composition_barplot.png）
      3C: 多层次评分小提琴图（3 张独立文件）

    ── Figure 4: 免疫抑制生态位特征基因提取与机制挖掘 ──
      4A: 空间差异基因火山图（fig4A_volcano.png）
      4B: Niche Marker 热图（fig4B_niche_marker_heatmap.png）
      4C: 通路富集气泡图（fig4C_pathway_bubble.png）

    ── Figure 5: CHC20 独立验证 ──
      5A: CHC20 域标签迁移预测图（fig5A_label_transfer.png）
      5B: HCC4R vs CHC20 验证对比小提琴图（fig5B_validation_violin.png）
      5C: 敏感性分析与参数扫描图（fig5C_sensitivity_params.png）

【使用方式】
    # 先确保 Step 1 和 Step 2 已运行（HCC4R 为主分析）：
    python code/run_hcc4r.py

    # 然后运行本脚本：
    python code/pipeline/run_paper_figures.py

    # 可指定不同的输入输出路径：
    python code/pipeline/run_paper_figures.py \\
        --hcc4r-niche-dir results/HCC4R/spatial_niche \\
        --chc20-niche-dir results/CHC20/spatial_niche \\
        --hcc4r-adata results/HCC4R/adata_vis_post.h5ad \\
        --chc20-adata results/CHC20/adata_vis_post.h5ad \\
        --scrna-adata results/HCC4R/adata_sc_post.h5ad \\
        --out-dir results/paper_figures

【依赖关系】
    上游：
      results/HCC4R/spatial_niche/spatial_niche_scores.csv
      results/HCC4R/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv
      results/HCC4R/adata_vis_post.h5ad
      results/HCC4R/adata_sc_post.h5ad
      results/CHC20/adata_vis_post.h5ad  (可选，用于 Figure 5)
      results/CHC20/spatial_niche/spatial_niche_scores.csv  (可选，用于 Figure 5B)
================================================================================
"""
from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Optional
import shutil

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
# 全局配色：Yellow-Black-Purple（仿论文风格）
# ============================================================
PAPER_YBP = LinearSegmentedColormap.from_list(
    "paper_ybp",
    ["#3a0063", "#000000", "#f5e642"],
    N=256,
)

# 高饱和度离散色板（用于细胞类型着色）
PAPER_CLUSTER_COLORS = [
    "#2ca02c", "#9467bd", "#1f77b4", "#d62728", "#ff7f0e",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]

# Spatial Domain 专用配色
DOMAIN_COLORS = {
    "tumor_core":              "#d62728",
    "tumor_edge":              "#ff7f0e",
    "stroma_immune":           "#1f77b4",
    "immunosuppressive_niche": "#9467bd",
    "other":                   "#bdbdbd",
}

# 细胞类型配色
CELLTYPE_COLORS = {
    "Hepatocyte":  "#1f77b4",
    "Malignant":   "#d62728",
    "T/NK":        "#2ca02c",
    "Treg":        "#9467bd",
    "Myeloid":     "#ff7f0e",
    "Fibroblast":  "#8c564b",
    "HSC":         "#8c564b",
    "B cell":      "#e377c2",
    "Endothelial": "#17becf",
    "Plasma cell": "#bcbd22",
    "Epithelial":  "#aec7e8",
}

PAPER_DPI = 300

_REPO_ROOT = Path(__file__).resolve().parents[2]


# ============================================================
# 参数解析
# ============================================================

def _parse_args() -> argparse.Namespace:
    results = _REPO_ROOT / "results"
    parser = argparse.ArgumentParser(
        description="Generate all paper figures for the spatial immunosuppressive niche study."
    )
    parser.add_argument("--hcc4r-niche-dir", type=Path,
                        default=results / "HCC4R" / "spatial_niche")
    parser.add_argument("--chc20-niche-dir", type=Path,
                        default=results / "CHC20" / "spatial_niche")
    parser.add_argument("--hcc4r-adata", type=Path,
                        default=results / "HCC4R" / "adata_vis_post.h5ad")
    parser.add_argument("--chc20-adata", type=Path,
                        default=results / "CHC20" / "adata_vis_post.h5ad")
    parser.add_argument("--scrna-adata", type=Path,
                        default=results / "HCC4R" / "adata_sc_post.h5ad")
    parser.add_argument("--hcc4r-visium-dir", type=Path, default=None,
                        help="HCC4R Visium Space Ranger 原始目录（用于读取 H&E 图）")
    parser.add_argument("--out-dir", type=Path,
                        default=results / "paper_figures")
    parser.add_argument("--dpi", type=int, default=PAPER_DPI)
    return parser.parse_args()


# ============================================================
# 日志配置
# ============================================================

def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


# ============================================================
# 通用工具函数
# ============================================================

def _save(fig: plt.Figure, path: Path, dpi: int = PAPER_DPI) -> None:
    fig.patch.set_facecolor("white")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logging.info("Saved: %s", path)


def _expression_frame(adata: ad.AnnData) -> pd.DataFrame:
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


def _get_abundance(adata: ad.AnnData) -> pd.DataFrame:
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
    return abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)


def _wilcoxon_test(high_vals: np.ndarray, low_vals: np.ndarray) -> float:
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
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


# ============================================================
# Figure 1A: 课题技术路线图
# ============================================================

def plot_fig1A_workflow(out_dir: Path, dpi: int = PAPER_DPI) -> None:
    """绘制三列布局的课题技术路线图。"""
    fig, ax = plt.subplots(figsize=(15, 9))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 9)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    c_scrna   = "#2166ac"
    c_spatial = "#d6604d"
    c_tcga    = "#4dac26"
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

    # ── scRNA-seq 左列 ──────────────────────────────────────────────
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
        (7.2, "HCC4R Visium Slide\n(Discovery Cohort)", "#d6604d"),
        (5.7, "Cell2location\nDeconvolution", "#888"),
        (4.2, "Spatial Niche Discovery\n(kNN + Leiden)", "#888"),
        (2.7, "CHC20 Visium Slide\n(Validation, Label Transfer)", "#f4a442"),
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
    """展示 HCC4R H&E 图像，若无原始图则用 spot 空间分布图代替。"""
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
                               c=color, s=9, alpha=0.8, linewidths=0, label=region)
            ax.legend(loc="upper right", fontsize=9, frameon=False,
                      title="Spatial Region", title_fontsize=10)
        elif "spatial_x" in spatial_df.columns:
            ax.scatter(spatial_df["spatial_x"], spatial_df["spatial_y"],
                       c="#c49c94", s=9, alpha=0.7, linewidths=0)
        ax.text(0.5, 0.02,
                "Placeholder: set --hcc4r-visium-dir to display actual H&E image",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
                color="#888888", style="italic")

    ax.set_title("HCC4R: H&E Histology (Spatial Slide Overview)", fontsize=13, fontweight="bold")
    ax.axis("off")
    _save(fig, out_dir / "fig2A_HE_image.png", dpi=dpi)


# ============================================================
# Figure 2B: 核心细胞类型空间丰度图（单张切片版）
# ============================================================

def plot_fig2B_spatial_abundance(
    spatial_df: pd.DataFrame,
    cell_types: list[str],
    out_dir: Path,
    sample_name: str = "HCC4R",
    dpi: int = PAPER_DPI,
) -> None:
    """为每种细胞类型单独生成空间丰度图（paper_ybp 配色）。"""
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
            x, y, c=v, s=9, cmap=PAPER_YBP,
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
    """计算并可视化细胞类型丰度间的 Pearson 相关性热图。"""
    cols = [c for c in cell_types if c in proportions.columns]
    if len(cols) < 2:
        logging.warning("Not enough cell types for correlation heatmap; skipping 2C.")
        return

    corr = proportions[cols].corr(method="pearson")

    # 推荐排序
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
    """以颜色展示各 Spatial Domain 的空间分布。"""
    if domain_col not in spatial_df.columns:
        for fallback in ["niche_semantic_label", "neighborhood_cluster"]:
            if fallback in spatial_df.columns:
                domain_col = fallback
                break
        else:
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
    """绘制各 Spatial Domain 的细胞类型组成 100% 堆叠条形图。"""
    if domain_col not in spatial_df.columns:
        for fallback in ["niche_semantic_label", "neighborhood_cluster"]:
            if fallback in spatial_df.columns:
                domain_col = fallback
                break
        else:
            logging.warning("No domain column; skipping 3B.")
            return

    valid_cols = [c for c in cell_types if c in spatial_df.columns]
    if not valid_cols:
        # 尝试 T/NK → T_NK
        valid_cols = [c.replace("/", "_") for c in cell_types
                      if c.replace("/", "_") in spatial_df.columns]
    if not valid_cols:
        logging.warning("No cell type columns found; skipping 3B.")
        return

    grouped = spatial_df.groupby(domain_col)[valid_cols].mean()
    row_sums = grouped.sum(axis=1)
    grouped_pct = grouped.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0) * 100

    # 排序 domain
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

    display_names = {
        "tumor_core":    "Tumor\nCore",
        "tumor_edge":    "Tumor\nEdge",
        "stroma_immune": "Stroma /\nImmune",
        "other":         "Other",
    }
    ax.set_xticks(range(len(grouped_pct)))
    ax.set_xticklabels([display_names.get(d, d) for d in grouped_pct.index],
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
    """为三个评分分别生成小提琴图，每张图单独保存。"""
    score_configs = [
        ("Treg_like_score",               "Treg-like Score",
         "fig3C_violin_treg_like_score.png",               "#9467bd"),
        ("immune_stroma_score",           "Immune-Stroma Score",
         "fig3C_violin_immune_stroma_score.png",           "#1f77b4"),
        ("immunosuppressive_niche_score", "Immunosuppressive Niche Score",
         "fig3C_violin_immunosuppressive_niche_score.png", "#d62728"),
    ]

    if domain_col not in spatial_df.columns:
        for fallback in ["niche_semantic_label", "neighborhood_cluster"]:
            if fallback in spatial_df.columns:
                domain_col = fallback
                break
        else:
            logging.warning("No domain column; skipping 3C.")
            return

    domain_order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    domain_order = [d for d in domain_order if d in spatial_df[domain_col].unique()]
    domain_order += sorted(d for d in spatial_df[domain_col].unique()
                           if d not in domain_order)

    display_names = {
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
        groups = [g for g in groups]   # keep all (even empty)

        fig, ax = plt.subplots(figsize=(max(7, len(domain_order) * 1.7 + 2), 5.5))

        vp = ax.violinplot(groups, positions=range(len(domain_order)),
                           widths=0.7, showmedians=True, showextrema=True)
        for pc in vp["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.55)
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
            g0 = groups[i0]
            g1 = groups[i1]
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
        ax.set_xticklabels([display_names.get(d, d) for d in domain_order],
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
    """绘制 niche_high vs niche_low 的全基因火山图（论文版）。"""
    highlight_genes = [
        "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
        "TGFB1", "IL10", "CXCL12", "CCL22",
        "FAP", "ACTA2", "POSTN", "COL1A1",
        "CD163", "MRC1", "SPP1", "VEGFA",
        "LAG3", "PDCD1", "HAVCR2",
    ]

    if "niche_high" not in niche_scores_df.columns:
        logging.warning("'niche_high' not found; skipping 4A.")
        return

    logging.info("Figure 4A: Computing differential expression...")
    expr_df = _expression_frame(adata)

    # 对齐 niche_high
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
        "#d62728" if (lfc.iloc[i] > fc_threshold and fdr[i] < fdr_threshold)
        else "#4575b4" if (lfc.iloc[i] < -fc_threshold and fdr[i] < fdr_threshold)
        else "#bdbdbd"
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
        mpatches.Patch(facecolor="#d62728", label=f"Up-regulated (n={n_up})"),
        mpatches.Patch(facecolor="#4575b4", label=f"Down-regulated (n={n_dn})"),
        mpatches.Patch(facecolor="#bdbdbd", label="Not significant"),
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
    """展示 Top-N 签名基因在各 Spatial Domain 中的平均表达热图。"""
    from scipy.stats import zscore as scipy_zscore

    # 读取签名基因
    if ranked_genes_path and ranked_genes_path.exists():
        ranked = pd.read_csv(ranked_genes_path)
        sort_col = "composite_score" if "composite_score" in ranked.columns else "log2_fc"
        sig_genes = (ranked.sort_values(sort_col, ascending=False)
                     .head(top_n)["gene"].tolist())
    else:
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

    if domain_col not in spatial_df.columns:
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
    """
    对各 Spatial Domain 进行预定义通路的基因集评分，以气泡图展示。

    X 轴：Spatial Domain
    Y 轴：通路名称
    气泡大小：富集评分（越大越富集）
    气泡颜色：富集方向（红=高于全局均值，蓝=低于全局均值）
    """
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

    if domain_col not in spatial_df.columns:
        for fallback in ["niche_semantic_label", "neighborhood_cluster"]:
            if fallback in spatial_df.columns:
                domain_col = fallback
                break
        else:
            logging.warning("No domain column; skipping 4C.")
            return

    # 表达矩阵
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

    # 计算各通路在各 domain 的平均评分
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
                "delta": s - global_mean,   # 相对全局均值的差
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

    # 气泡图
    n_pathways = len(pathway_order)
    n_domains  = len(domain_order)

    fig, ax = plt.subplots(figsize=(max(8, n_domains * 1.8 + 2),
                                    max(6, n_pathways * 0.7 + 2)))

    # 归一化气泡大小
    all_scores = result_df["score"].to_numpy()
    s_min, s_max = all_scores.min(), all_scores.max()
    s_range = s_max - s_min if s_max > s_min else 1.0

    domain_display = {
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
    ax.set_xticklabels([domain_display.get(d, d) for d in domain_order],
                        fontsize=10.5, fontweight="bold")
    ax.set_yticks(range(n_pathways))
    ax.set_yticklabels(pathway_order, fontsize=9)
    ax.set_xlim(-0.5, n_domains - 0.5)
    ax.set_ylim(-0.5, n_pathways - 0.5)
    ax.invert_yaxis()

    # 图例
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
    for y in range(n_pathways):
        ax.axhline(y, color="white", lw=0.8, zorder=0)
    for x in range(n_domains):
        ax.axvline(x, color="white", lw=0.8, zorder=0)

    _save(fig, out_dir / "fig4C_pathway_bubble.png", dpi=dpi)


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
    """
    用 HCC4R 训练随机森林分类器（特征=细胞类型比例），
    在 CHC20 上预测 Spatial Domain 标签，并可视化空间分布。
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

    # 获取 HCC4R 丰度特征
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

    # 对齐特征列
    common_ct = [c for c in hcc4r_prop.columns if c in chc20_prop.columns]
    if not common_ct:
        logging.warning("No common cell types between HCC4R and CHC20; skipping 5A.")
        return

    # 构建训练数据
    if "spot_id" in hcc4r_niche_df.columns:
        label_map = hcc4r_niche_df.set_index("spot_id")[domain_col]
    else:
        label_map = hcc4r_niche_df[domain_col]

    y_train = label_map.reindex(hcc4r_prop.index).fillna("other")
    X_train = hcc4r_prop[common_ct].fillna(0)

    # 训练
    le = LabelEncoder()
    y_enc = le.fit_transform(y_train.to_numpy())
    clf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1,
                                 class_weight="balanced")
    clf.fit(X_train.to_numpy(), y_enc)

    # 预测 CHC20
    X_test = chc20_prop[common_ct].fillna(0)
    y_pred_enc = clf.predict(X_test.to_numpy())
    y_pred = le.inverse_transform(y_pred_enc)

    # 空间坐标（CHC20）
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

    display_names = {
        "tumor_core": "Tumor Core", "tumor_edge": "Tumor Edge",
        "stroma_immune": "Stroma/Immune",
    }
    for lbl in unique_labels:
        mask = y_pred == lbl
        if mask.sum() == 0:
            continue
        display = display_names.get(lbl, lbl)
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
    hcc4r_niche_df: pd.DataFrame,
    chc20_niche_df: pd.DataFrame,
    out_dir: Path,
    score_cols: Optional[list[str]] = None,
    dpi: int = PAPER_DPI,
) -> None:
    """
    并排小提琴图：比较 HCC4R 和 CHC20 在关键评分上的分布一致性。

    若两个 CSV 都有评分数据，对每个 score 生成单独图。
    """
    if score_cols is None:
        score_cols = [
            "immunosuppressive_niche_score",
            "Treg_like_score",
            "immune_stroma_score",
        ]

    display_names = {
        "immunosuppressive_niche_score": "Immunosuppressive Niche Score",
        "Treg_like_score":               "Treg-like Score",
        "immune_stroma_score":           "Immune-Stroma Score",
    }

    for score in score_cols:
        hcc4r_vals = hcc4r_niche_df[score].dropna().to_numpy() \
            if score in hcc4r_niche_df.columns else np.array([])
        chc20_vals = chc20_niche_df[score].dropna().to_numpy() \
            if score in chc20_niche_df.columns else np.array([])

        if len(hcc4r_vals) == 0 and len(chc20_vals) == 0:
            logging.warning("Score '%s' not in either niche df; skipping.", score)
            continue

        # Wilcoxon 检验（双侧，比较两组分布）
        if len(hcc4r_vals) > 1 and len(chc20_vals) > 1:
            _, p_two = scipy_stats.mannwhitneyu(
                hcc4r_vals, chc20_vals, alternative="two-sided"
            )
        else:
            p_two = 1.0

        fig, ax = plt.subplots(figsize=(7, 5.5))

        groups = []
        positions = []
        colors_vl = []
        labels = []
        if len(hcc4r_vals) > 0:
            groups.append(hcc4r_vals)
            positions.append(0)
            colors_vl.append("#d6604d")
            labels.append("HCC4R\n(Discovery)")
        if len(chc20_vals) > 0:
            groups.append(chc20_vals)
            positions.append(1)
            colors_vl.append("#4dac26")
            labels.append("CHC20\n(Validation)")

        if not groups:
            plt.close(fig)
            continue

        vp = ax.violinplot(groups, positions=positions,
                           widths=0.6, showmedians=True, showextrema=True)
        for i, pc in enumerate(vp["bodies"]):
            pc.set_facecolor(colors_vl[i])
            pc.set_alpha(0.6)
            pc.set_edgecolor("white")
        vp["cmedians"].set_color("#222")
        vp["cmedians"].set_linewidth(2.2)

        ax.boxplot(groups, positions=positions, widths=0.1, patch_artist=True,
                   medianprops=dict(color="black", linewidth=2),
                   boxprops=dict(facecolor="white", alpha=0.9, linewidth=1),
                   whiskerprops=dict(linewidth=0), capprops=dict(linewidth=0),
                   flierprops=dict(marker=""))

        # 显著性标注
        if len(groups) == 2:
            all_vals = np.concatenate(groups)
            y_max = np.percentile(all_vals, 99)
            y_span = np.percentile(all_vals, 99) - np.percentile(all_vals, 1)
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
        ax.set_title(f"{score_name}\nHCC4R vs CHC20 Cross-Cohort Consistency",
                     fontsize=12, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=9)

        safe = score.replace("/", "_")
        _save(fig, out_dir / f"fig5B_validation_violin_{safe}.png", dpi=dpi)


# ============================================================
# Figure 5C: 敏感性分析与参数扫描整合图（两张独立图）
# ============================================================

def plot_fig5C_sensitivity(
    hcc4r_niche_dir: Path,
    out_dir: Path,
    dpi: int = PAPER_DPI,
) -> None:
    """
    直接使用 Step 2 已生成的敏感性分析和参数扫描结果，
    生成论文版整洁图（重新绘制，统一 DPI 和风格）。
    """
    # ── 5C-1: 敏感性分析图（Spearman ρ + Weighted Jaccard）─────────────────
    sens_path = hcc4r_niche_dir / "sensitivity_analysis.csv"
    if sens_path.exists():
        sens_df = pd.read_csv(sens_path)
        _plot_sensitivity_paper(sens_df, out_dir / "fig5C_sensitivity_stability.png", dpi=dpi)
    else:
        logging.warning("sensitivity_analysis.csv not found; skipping 5C-1.")

    # ── 5C-2: 参数扫描热图（k × quantile）───────────────────────────────────
    scan_path = hcc4r_niche_dir / "param_scan_deg_stability.csv"
    if scan_path.exists():
        scan_df = pd.read_csv(scan_path)
        _plot_param_scan_paper(scan_df, out_dir / "fig5C_param_scan_heatmap.png", dpi=dpi)
    else:
        logging.warning("param_scan_deg_stability.csv not found; skipping 5C-2.")


def _plot_sensitivity_paper(sens_df: pd.DataFrame, path: Path, dpi: int = PAPER_DPI) -> None:
    """论文版敏感性分析图：2×2 布局（kNN + radius × Spearman ρ + Weighted Jaccard）。"""
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


def _plot_param_scan_paper(scan_df: pd.DataFrame, path: Path, dpi: int = PAPER_DPI) -> None:
    """论文版参数扫描热图（k × quantile，颜色编码 DEG 数量）。"""
    if scan_df.empty or "n_sig_deg" not in scan_df.columns:
        logging.warning("Empty param scan df; skipping 5C-2.")
        return

    pivot = scan_df.pivot(index="k", columns="quantile", values="n_sig_deg")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    sns.heatmap(pivot, cmap="YlOrRd", annot=True, fmt="d",
                linewidths=0.5, ax=ax,
                cbar_kws={"label": "# Significant DEGs (FDR<0.05, log₂FC>0.5)",
                           "shrink": 0.8})
    ax.set_xlabel("niche_high quantile threshold", fontsize=10)
    ax.set_ylabel("kNN neighbors (k)", fontsize=10)
    ax.set_title("Parameter Stability: k × quantile Grid Search\n"
                 "(darker = more stable DEGs; optimal = darkest plateau region)",
                 fontsize=11, fontweight="bold")

    _save(fig, path, dpi=dpi)


# ============================================================
# 将 Step 1 已生成的图复制/链接到 paper_figures 目录
# ============================================================

def _copy_existing_figure(src: Path, dst: Path) -> None:
    """如果源文件存在，复制到目标路径。"""
    if src.exists():
        shutil.copy2(src, dst)
        logging.info("Copied existing figure: %s → %s", src.name, dst.name)
    else:
        logging.warning("Expected figure not found (will be generated by Step 1): %s", src)


# ============================================================
# 主流程
# ============================================================

def main() -> int:
    _setup_logging()
    args = _parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    logging.info("Paper figures output directory: %s", out_dir)

    dpi = args.dpi

    # ============================================================
    # 加载 HCC4R 核心数据
    # ============================================================
    hcc4r_niche_csv = args.hcc4r_niche_dir / "spatial_niche_scores.csv"
    ranked_genes_csv = args.hcc4r_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv"

    hcc4r_niche_df: Optional[pd.DataFrame] = None
    if hcc4r_niche_csv.exists():
        hcc4r_niche_df = pd.read_csv(hcc4r_niche_csv)
        logging.info("Loaded HCC4R niche scores: %d spots", len(hcc4r_niche_df))
    else:
        logging.warning(
            "HCC4R niche scores not found: %s\n"
            "Please run Step 2 first: python code/pipeline/run_spatial_niche_analysis.py "
            "--adata results/HCC4R/adata_vis_post.h5ad --out-dir results/HCC4R/spatial_niche",
            hcc4r_niche_csv,
        )

    hcc4r_adata: Optional[ad.AnnData] = None
    if args.hcc4r_adata.exists():
        logging.info("Loading HCC4R AnnData: %s", args.hcc4r_adata)
        hcc4r_adata = ad.read_h5ad(args.hcc4r_adata)
        logging.info("HCC4R AnnData: %d spots × %d genes", hcc4r_adata.n_obs, hcc4r_adata.n_vars)
    else:
        logging.warning("HCC4R AnnData not found: %s", args.hcc4r_adata)

    hcc4r_prop_df: Optional[pd.DataFrame] = None
    if hcc4r_adata is not None:
        try:
            hcc4r_abund = _get_abundance(hcc4r_adata)
            hcc4r_prop_df = _proportions(hcc4r_abund)
        except Exception as e:
            logging.warning("Could not get HCC4R abundance: %s", e)

    # CHC20 数据（可选）
    chc20_adata: Optional[ad.AnnData] = None
    chc20_niche_df: Optional[pd.DataFrame] = None
    chc20_niche_csv = args.chc20_niche_dir / "spatial_niche_scores.csv"

    if args.chc20_adata.exists():
        logging.info("Loading CHC20 AnnData: %s", args.chc20_adata)
        chc20_adata = ad.read_h5ad(args.chc20_adata)

    if chc20_niche_csv.exists():
        chc20_niche_df = pd.read_csv(chc20_niche_csv)
        logging.info("Loaded CHC20 niche scores: %d spots", len(chc20_niche_df))

    # 决定要使用的细胞类型列表
    if hcc4r_prop_df is not None:
        all_cell_types = hcc4r_prop_df.columns.tolist()
    elif hcc4r_niche_df is not None:
        skip_cols = {"spot_id", "sample", "spatial_x", "spatial_y", "neighbor_radius",
                     "immunosuppressive_gene_score", "neighborhood_cluster",
                     "niche_semantic_label", "hep_high", "distance_to_hep_high",
                     "has_hep_high_neighbor", "Treg_like_score", "immune_stroma_score",
                     "immunosuppressive_niche_score", "niche_high", "spatial_region"}
        all_cell_types = [c for c in hcc4r_niche_df.columns
                          if c not in skip_cols and not c.startswith("neighbor_")
                          and not c.startswith("rank_")]
    else:
        all_cell_types = ["Hepatocyte", "Malignant", "Treg", "T/NK",
                          "Myeloid", "Fibroblast", "HSC", "B cell"]

    # 2B 展示的重点细胞类型
    key_cell_types_2B = ["Malignant", "HSC", "Treg", "Myeloid"]

    # ============================================================
    # ── Figure 1 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 1 ===")

    # 1A: 技术路线图
    logging.info("Figure 1A: Workflow diagram...")
    plot_fig1A_workflow(out_dir, dpi=dpi)

    # 1B/1C/1D: 复制 Step 1 已生成的图
    if hcc4r_niche_df is not None:
        hcc4r_results_dir = args.hcc4r_niche_dir.parent
        for src_name, dst_name in [
            ("scrna_tsne_celltype.png",          "fig1B_scrna_tsne_celltype.png"),
            ("scrna_celltype_marker_heatmap.png","fig1C_scrna_marker_heatmap.png"),
            ("t_cell_dotplot_horizontal.png",    "fig1D_treg_dotplot.png"),
        ]:
            _copy_existing_figure(hcc4r_results_dir / src_name, out_dir / dst_name)

    # ============================================================
    # ── Figure 2 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 2 ===")

    # 2A: H&E 图
    logging.info("Figure 2A: H&E image...")
    plot_fig2A_HE(args.hcc4r_visium_dir, hcc4r_niche_df, out_dir, dpi=dpi)

    # 2B: 空间丰度图
    if hcc4r_niche_df is not None:
        logging.info("Figure 2B: Spatial abundance maps...")
        plot_fig2B_spatial_abundance(
            hcc4r_niche_df,
            cell_types=key_cell_types_2B,
            out_dir=out_dir,
            sample_name="HCC4R",
            dpi=dpi,
        )
    else:
        logging.warning("HCC4R niche df not available; skipping Figure 2B.")

    # 2C: 共定位热图
    if hcc4r_prop_df is not None:
        logging.info("Figure 2C: Co-localization heatmap...")
        plot_fig2C_coloc_heatmap(
            hcc4r_prop_df,
            cell_types=all_cell_types,
            out_dir=out_dir,
            sample_name="HCC4R",
            dpi=dpi,
        )
    elif hcc4r_niche_df is not None:
        # 用 niche_df 中的细胞类型列代替
        logging.info("Figure 2C: Co-localization heatmap (from niche df)...")
        prop_cols = [c for c in all_cell_types if c in hcc4r_niche_df.columns]
        if prop_cols:
            plot_fig2C_coloc_heatmap(
                hcc4r_niche_df[prop_cols],
                cell_types=prop_cols,
                out_dir=out_dir,
                sample_name="HCC4R",
                dpi=dpi,
            )

    # ============================================================
    # ── Figure 3 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 3 ===")

    if hcc4r_niche_df is not None:
        logging.info("Figure 3A: Spatial domain map...")
        plot_fig3A_spatial_domain(hcc4r_niche_df, out_dir, dpi=dpi)

        logging.info("Figure 3B: Domain composition barplot...")
        plot_fig3B_domain_composition(
            hcc4r_niche_df, cell_types=all_cell_types,
            out_dir=out_dir, dpi=dpi,
        )

        logging.info("Figure 3C: Score violin plots...")
        plot_fig3C_score_violins(hcc4r_niche_df, out_dir, dpi=dpi)
    else:
        logging.warning("HCC4R niche df not available; skipping Figure 3.")

    # ============================================================
    # ── Figure 4 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 4 ===")

    if hcc4r_niche_df is not None and hcc4r_adata is not None:
        logging.info("Figure 4A: Volcano plot...")
        plot_fig4A_volcano(hcc4r_niche_df, hcc4r_adata, out_dir, dpi=dpi)

        logging.info("Figure 4B: Niche marker heatmap...")
        plot_fig4B_niche_marker_heatmap(
            hcc4r_adata, hcc4r_niche_df, ranked_genes_csv,
            out_dir=out_dir, top_n=30, dpi=dpi,
        )

        logging.info("Figure 4C: Pathway bubble plot...")
        plot_fig4C_pathway_bubble(
            hcc4r_adata, hcc4r_niche_df,
            out_dir=out_dir, dpi=dpi,
        )
    else:
        logging.warning("HCC4R data not available; skipping Figure 4.")

    # ============================================================
    # ── Figure 5 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 5 ===")

    # 5A: 标签迁移
    if hcc4r_niche_df is not None and hcc4r_adata is not None and chc20_adata is not None:
        logging.info("Figure 5A: Label transfer to CHC20...")
        plot_fig5A_label_transfer(
            hcc4r_niche_df, hcc4r_adata, chc20_adata,
            out_dir=out_dir, dpi=dpi,
        )
    else:
        logging.warning("Missing data for Figure 5A; skipping.")

    # 5B: 跨队列验证小提琴图
    if hcc4r_niche_df is not None and chc20_niche_df is not None:
        logging.info("Figure 5B: Cross-cohort validation violin plots...")
        plot_fig5B_validation_violin(
            hcc4r_niche_df, chc20_niche_df,
            out_dir=out_dir, dpi=dpi,
        )
    else:
        logging.warning("CHC20 niche df not available; skipping Figure 5B.")

    # 5C: 敏感性分析
    logging.info("Figure 5C: Sensitivity and parameter scan...")
    plot_fig5C_sensitivity(args.hcc4r_niche_dir, out_dir, dpi=dpi)

    # ============================================================
    # 汇总
    # ============================================================
    generated = sorted(out_dir.glob("*.png"))
    logging.info(
        "\n=== Paper figures generation complete ===\n"
        "Output directory : %s\n"
        "Total figures    : %d\n"
        "File list:\n  %s",
        out_dir,
        len(generated),
        "\n  ".join(f.name for f in generated),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
