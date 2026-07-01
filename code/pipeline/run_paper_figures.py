"""
================================================================================
脚本名称: run_paper_figures.py
功能概述: 论文图表数据准备与流程编排（Figure 1–5）
================================================================================

【职责说明】
    本脚本只负责：
      1. 解析命令行参数
      2. 加载各类结果数据（AnnData、niche CSV 等）
      3. 按 Figure 1 → 5 顺序调用 paper_plot_functions.py 中的绘图函数

    ★ 若需修改绘图样式、颜色、字号、图幅大小，请编辑：
      code/pipeline/paper_plot_functions.py

【Figure 5 说明（新版：双模块评分验证）】
    Figure 5 已从单一 signature score 升级为双模块评分体系：
      · Fig.5A: CHC20 免疫模块空间评分（immune_signature_score 投影到空间坐标）
      · Fig.5B: CHC20 基质/ECM 模块空间评分（stromal_ECM_signature_score 投影）
      · Fig.5C: CHC20 组合 niche score（z_immune + z_stromal，反映免疫抑制微环境）
      · Fig.5D: HCC4R vs CHC20 combined niche score 分布对比（小提琴+箱线）

    模块基因来源（优先级由高到低）：
      1. ranked_genes_csv 中存在 module 列（"immune"/"stromal"/"ecm"）→ 自动读取
      2. 使用脚本内预设的先验基因列表（_IMMUNE_GENES_DEFAULT / _STROMAL_GENES_DEFAULT）

    输出 CSV 文件：
      fig5_dual_module_scores_HCC4R.csv  — HCC4R 三列评分（免疫/基质/组合）
      fig5_dual_module_scores_CHC20.csv  — CHC20 三列评分（免疫/基质/组合）
      fig5_immune_genes_used_CHC20.txt   — CHC20 实际匹配到的免疫模块基因
      fig5_stromal_genes_used_CHC20.txt  — CHC20 实际匹配到的基质模块基因

【上游依赖】
    results/HCC4R/spatial_niche/spatial_niche_scores.csv                 (Step 2)
    results/HCC4R/spatial_niche/immunosuppressive_niche_signature_genes_ranked.csv
    results/HCC4R/adata_vis_post.h5ad                                    (Step 1)
    results/HCC4R/adata_sc_post.h5ad                                     (Step 1)
    results/CHC20/adata_vis_post.h5ad                            (可选，Figure 5A/5B/5C)
    results/HCC4R/spatial_niche/immunosuppressive_niche_signature_genes.txt

【使用方式】
    # 先运行 Step 1 + Step 2（以 HCC4R 为主）：
    python code/run_hcc4r.py

    # 再生成论文图表：
    python code/pipeline/run_paper_figures.py

    # 或通过入口脚本一键运行：
    python code/run_hcc4r.py --paper-figures
    python code/run_hcc4r.py --figures-only
================================================================================
"""
from __future__ import annotations

import argparse
import logging
import shutil
import sys
import warnings
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd

warnings.filterwarnings("ignore")

# 导入所有绘图函数（修改样式请去 paper_plot_functions.py）
from paper_plot_functions import (  # noqa: E402
    PAPER_DPI,
    plot_fig1A_workflow,
    plot_fig2A_HE,
    plot_fig2B_spatial_abundance,
    plot_fig2C_coloc_heatmap,
    plot_fig3A_spatial_domain,
    plot_fig3B_domain_composition,
    plot_fig3C_score_violins,
    plot_fig4A_volcano,
    plot_fig4B_niche_marker_heatmap,
    plot_fig4C_pathway_bubble,
    # ── Figure 5 新版：双模块评分 ──────────────────────────────────
    compute_dual_module_scores,
    plot_fig5A_immune_spatial,
    plot_fig5B_stromal_spatial,
    plot_fig5C_combined_niche_spatial,
    plot_fig5D_comparison_violin,
    # ── Figure 5 旧版（向后兼容）───────────────────────────────────
    compute_signature_score,
    plot_fig5A_signature_projection,
    plot_fig5B_validation_violin,
    plot_fig5C_sensitivity,
    _get_abundance,
    _proportions,
)

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
    parser.add_argument("--hcc4r-visium-dir", type=Path,
                        default=_REPO_ROOT / "data" / "HCC4R",
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
# 辅助：复制 Step 1 已生成的图
# ============================================================

def _copy_existing_figure(src: Path, dst: Path) -> None:
    """如果源文件存在则复制到目标路径，否则仅打印警告。"""
    if src.exists():
        shutil.copy2(src, dst)
        logging.info("Copied existing figure: %s → %s", src.name, dst.name)
    else:
        logging.warning("Expected figure not found (will be generated by Step 1): %s", src)


# ============================================================
# 数据加载
# ============================================================

def _load_data(args: argparse.Namespace):
    """加载所有上游结果文件，返回各数据对象（不存在时为 None）。"""

    # ── HCC4R niche scores CSV ──
    hcc4r_niche_csv = args.hcc4r_niche_dir / "spatial_niche_scores.csv"
    hcc4r_niche_df: Optional[pd.DataFrame] = None
    if hcc4r_niche_csv.exists():
        hcc4r_niche_df = pd.read_csv(hcc4r_niche_csv)
        logging.info("Loaded HCC4R niche scores: %d spots", len(hcc4r_niche_df))
    else:
        logging.warning(
            "HCC4R niche scores not found: %s\n"
            "  → Please run Step 2 first:\n"
            "    python code/pipeline/run_spatial_niche_analysis.py "
            "--adata results/HCC4R/adata_vis_post.h5ad "
            "--out-dir results/HCC4R/spatial_niche",
            hcc4r_niche_csv,
        )

    # ── HCC4R AnnData（空间转录组）──
    hcc4r_adata: Optional[ad.AnnData] = None
    if args.hcc4r_adata.exists():
        logging.info("Loading HCC4R AnnData: %s", args.hcc4r_adata)
        hcc4r_adata = ad.read_h5ad(args.hcc4r_adata)
        logging.info("HCC4R AnnData: %d spots × %d genes",
                     hcc4r_adata.n_obs, hcc4r_adata.n_vars)
    else:
        logging.warning("HCC4R AnnData not found: %s", args.hcc4r_adata)

    # ── HCC4R 细胞类型比例矩阵 ──
    hcc4r_prop_df: Optional[pd.DataFrame] = None
    if hcc4r_adata is not None:
        try:
            hcc4r_prop_df = _proportions(_get_abundance(hcc4r_adata))
        except Exception as e:
            logging.warning("Could not get HCC4R cell abundance: %s", e)

    # ── CHC20 AnnData（可选，用于 Figure 5）──
    chc20_adata: Optional[ad.AnnData] = None
    if args.chc20_adata.exists():
        logging.info("Loading CHC20 AnnData: %s", args.chc20_adata)
        chc20_adata = ad.read_h5ad(args.chc20_adata)

    # ── CHC20 niche scores CSV（可选）──
    chc20_niche_df: Optional[pd.DataFrame] = None
    chc20_niche_csv = args.chc20_niche_dir / "spatial_niche_scores.csv"
    if chc20_niche_csv.exists():
        chc20_niche_df = pd.read_csv(chc20_niche_csv)
        logging.info("Loaded CHC20 niche scores: %d spots", len(chc20_niche_df))

    return hcc4r_niche_df, hcc4r_adata, hcc4r_prop_df, chc20_adata, chc20_niche_df


def _load_signature_genes(args: argparse.Namespace, ranked_genes_csv: Path) -> list[str]:
    """Load HCC4R-derived immunosuppressive niche signature genes."""
    text_candidates = [
        args.hcc4r_niche_dir / "immunosuppressive_niche_signature_genes.txt",
        args.hcc4r_niche_dir.parent / "spatial_signature_genes.txt",
    ]
    for path in text_candidates:
        if path.exists():
            genes = [g.strip() for g in path.read_text(encoding="utf-8").splitlines()]
            genes = [g for g in genes if g]
            if genes:
                logging.info("Loaded HCC4R signature genes: %d from %s", len(genes), path)
                return genes

    if ranked_genes_csv.exists():
        ranked = pd.read_csv(ranked_genes_csv)
        for col in ["gene", "gene_symbol", "names", "index"]:
            if col in ranked.columns:
                genes = ranked[col].dropna().astype(str).str.strip().tolist()
                genes = [g for g in genes if g]
                if genes:
                    logging.info(
                        "Loaded HCC4R signature genes: %d from ranked CSV column '%s'",
                        len(genes), col,
                    )
                    return genes

    logging.warning(
        "No HCC4R signature gene list found. Checked: %s",
        ", ".join(str(p) for p in text_candidates + [ranked_genes_csv]),
    )
    return []


def _resolve_cell_types(
    hcc4r_prop_df: Optional[pd.DataFrame],
    hcc4r_niche_df: Optional[pd.DataFrame],
) -> tuple[list[str], list[str]]:
    """
    确定细胞类型列表。
    返回 (all_cell_types, key_cell_types_2B)
      - all_cell_types   : 用于 2C / 3B 等所有细胞类型
      - key_cell_types_2B: 用于 2B 空间丰度图的重点细胞类型
    """
    if hcc4r_prop_df is not None:
        all_cell_types = hcc4r_prop_df.columns.tolist()
    elif hcc4r_niche_df is not None:
        skip_cols = {
            "spot_id", "sample", "spatial_x", "spatial_y", "neighbor_radius",
            "immunosuppressive_gene_score", "neighborhood_cluster",
            "niche_semantic_label", "hep_high", "distance_to_hep_high",
            "has_hep_high_neighbor", "Treg_like_score", "immune_stroma_score",
            "immunosuppressive_niche_score", "niche_high", "spatial_region",
        }
        all_cell_types = [
            c for c in hcc4r_niche_df.columns
            if c not in skip_cols
            and not c.startswith("neighbor_")
            and not c.startswith("rank_")
        ]
    else:
        all_cell_types = [
            "Hepatocyte", "Malignant", "Treg", "T/NK",
            "Myeloid", "Fibroblast", "HSC", "B cell",
        ]

    # 2B 只展示最核心的 4 种
    key_cell_types_2B = ["Malignant", "HSC", "Treg", "Myeloid"]

    return all_cell_types, key_cell_types_2B


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

    # ── 数据加载 ──────────────────────────────────────────────
    (hcc4r_niche_df,
     hcc4r_adata,
     hcc4r_prop_df,
     chc20_adata,
     chc20_niche_df) = _load_data(args)

    ranked_genes_csv = args.hcc4r_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv"
    all_cell_types, key_cell_types_2B = _resolve_cell_types(hcc4r_prop_df, hcc4r_niche_df)

    # ============================================================
    # ── Figure 1 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 1 ===")

    logging.info("Figure 1A: Study workflow schematic...")
    plot_fig1A_workflow(out_dir, dpi=dpi)

    # 1C / 1D: 复制 Step 1 已生成的图
    if hcc4r_niche_df is not None:
        hcc4r_results_dir = args.hcc4r_niche_dir.parent
        for src_name, dst_name in [
            ("scrna_celltype_marker_heatmap.png", "fig1C_scrna_marker_heatmap.png"),
            ("t_cell_dotplot_horizontal.png",     "fig1D_treg_dotplot.png"),
        ]:
            _copy_existing_figure(hcc4r_results_dir / src_name, out_dir / dst_name)

    # ============================================================
    # ── Figure 2 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 2 ===")

    # 2A: H&E 图（需要 Visium 原始目录；无则用 placeholder）
    logging.info("Figure 2A: H&E image...")
    if hcc4r_niche_df is not None:
        plot_fig2A_HE(args.hcc4r_visium_dir, hcc4r_niche_df, out_dir, dpi=dpi)
    else:
        logging.warning("HCC4R niche df not available; skipping Figure 2A.")

    # 2B: 各细胞类型空间丰度散点图（每种独立文件）
    if hcc4r_niche_df is not None:
        logging.info("Figure 2B: Spatial abundance maps for %s...", key_cell_types_2B)
        plot_fig2B_spatial_abundance(
            hcc4r_niche_df,
            cell_types=key_cell_types_2B,
            out_dir=out_dir,
            sample_name="HCC4R",
            dpi=dpi,
        )
    else:
        logging.warning("HCC4R niche df not available; skipping Figure 2B.")

    # 2C: 细胞类型共定位相关性热图
    if hcc4r_prop_df is not None:
        logging.info("Figure 2C: Co-localization heatmap (from abundance)...")
        plot_fig2C_coloc_heatmap(
            hcc4r_prop_df,
            cell_types=all_cell_types,
            out_dir=out_dir,
            sample_name="HCC4R",
            dpi=dpi,
        )
    elif hcc4r_niche_df is not None:
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
    else:
        logging.warning("No abundance data; skipping Figure 2C.")

    # ============================================================
    # ── Figure 3 ──────────────────────────────────────────────
    # ============================================================
    logging.info("=== Figure 3 ===")

    if hcc4r_niche_df is not None:
        logging.info("Figure 3A: Spatial domain map...")
        plot_fig3A_spatial_domain(hcc4r_niche_df, out_dir, dpi=dpi)

        logging.info("Figure 3B: Domain composition barplot...")
        plot_fig3B_domain_composition(
            hcc4r_niche_df,
            cell_types=all_cell_types,
            out_dir=out_dir,
            dpi=dpi,
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

        logging.info("Figure 4B: Niche marker heatmap (top 30)...")
        plot_fig4B_niche_marker_heatmap(
            hcc4r_adata,
            hcc4r_niche_df,
            ranked_genes_csv,
            out_dir=out_dir,
            top_n=30,
            dpi=dpi,
        )

        logging.info("Figure 4C: Pathway bubble plot...")
        plot_fig4C_pathway_bubble(hcc4r_adata, hcc4r_niche_df, out_dir, dpi=dpi)
    else:
        logging.warning("HCC4R data not available; skipping Figure 4.")

    # ============================================================
    # ── Figure 5（新版：双模块评分验证）─────────────────────────
    # ============================================================
    logging.info("=== Figure 5 ===")

    # ── 从签名基因文件中拆分免疫模块 / 基质模块基因 ──────────────
    # 规则：先尝试从 ranked_genes_csv 读取 module 列（immune/stromal）；
    #       若不存在 module 列，则用预设的先验基因列表作为两个子模块。
    def _split_module_genes(
        ranked_csv: Path,
    ) -> tuple[list[str], list[str]]:
        """从 ranked CSV 中拆分免疫/基质模块基因，或返回预设列表。"""
        # ★ 预设先验基因：可在此直接修改
        _IMMUNE_GENES_DEFAULT = [
            "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF18",
            "IL10", "TGFB1", "PDCD1", "LAG3", "HAVCR2", "ENTPD1",
            "CD163", "MRC1", "ARG1", "CCL22", "CCL17", "VSIR",
            "IDO1", "IDO2", "CXCL12",
        ]
        _STROMAL_GENES_DEFAULT = [
            "FAP", "ACTA2", "POSTN", "COL1A1", "COL1A2", "COL3A1",
            "COL4A1", "VCAN", "FN1", "LOXL2", "THBS2", "SPP1",
            "MMP2", "MMP9", "MMP11", "SPARC", "LUM", "DCN",
            "PDGFRA", "PDGFRB", "PECAM1",
        ]
        if ranked_csv.exists():
            try:
                ranked = pd.read_csv(ranked_csv)
                if "module" in ranked.columns and "gene" in ranked.columns:
                    immune_genes = (
                        ranked.loc[ranked["module"].str.lower() == "immune", "gene"]
                        .dropna().str.strip().tolist()
                    )
                    stromal_genes = (
                        ranked.loc[ranked["module"].str.lower().isin(["stromal", "ecm", "stroma"]),
                                   "gene"]
                        .dropna().str.strip().tolist()
                    )
                    if immune_genes and stromal_genes:
                        logging.info(
                            "Module genes loaded from ranked CSV: "
                            "immune=%d, stromal=%d",
                            len(immune_genes), len(stromal_genes),
                        )
                        return immune_genes, stromal_genes
            except Exception as _e:
                logging.warning("Could not parse module column from ranked CSV: %s", _e)

        logging.info(
            "Using preset module gene lists: immune=%d, stromal=%d",
            len(_IMMUNE_GENES_DEFAULT), len(_STROMAL_GENES_DEFAULT),
        )
        return _IMMUNE_GENES_DEFAULT, _STROMAL_GENES_DEFAULT

    immune_genes, stromal_genes = _split_module_genes(ranked_genes_csv)

    # ── HCC4R 双模块评分（用于 Fig.5D 对比）──────────────────────
    hcc4r_combined_score: pd.Series = pd.Series(dtype=float)
    if hcc4r_adata is not None:
        logging.info("Computing dual-module scores for HCC4R (discovery)...")
        (hcc4r_immune, hcc4r_stromal, hcc4r_combined_score,
         hcc4r_m_immune, hcc4r_m_stromal) = compute_dual_module_scores(
            hcc4r_adata, immune_genes, stromal_genes,
        )
        if not hcc4r_combined_score.empty:
            # 保存 HCC4R 双模块评分 CSV（供后续分析复用）
            hcc4r_module_csv = out_dir / "fig5_dual_module_scores_HCC4R.csv"
            pd.DataFrame({
                "spot_id":                  hcc4r_combined_score.index,
                "sample":                   "HCC4R",
                "immune_signature_score":   hcc4r_immune.reindex(hcc4r_combined_score.index),
                "stromal_ECM_signature_score": hcc4r_stromal.reindex(hcc4r_combined_score.index),
                "immune_stromal_niche_score": hcc4r_combined_score.to_numpy(dtype=float),
            }).to_csv(hcc4r_module_csv, index=False)
            logging.info("Saved HCC4R dual-module scores: %s", hcc4r_module_csv)
        else:
            logging.warning("HCC4R combined score is empty; Fig.5D will be skipped.")
    else:
        logging.warning("HCC4R AnnData not available; skipping HCC4R dual-module scores.")

    # ── CHC20 双模块评分 + 空间投影图（Fig.5A / 5B / 5C）────────
    chc20_combined_score: pd.Series = pd.Series(dtype=float)
    if chc20_adata is not None:
        logging.info("Computing dual-module scores for CHC20 (validation)...")
        (chc20_immune, chc20_stromal, chc20_combined_score,
         chc20_m_immune, chc20_m_stromal) = compute_dual_module_scores(
            chc20_adata, immune_genes, stromal_genes,
        )

        if not chc20_immune.empty:
            logging.info("Figure 5A: CHC20 immune module spatial score...")
            plot_fig5A_immune_spatial(
                chc20_adata, chc20_immune, out_dir=out_dir,
                matched_genes=chc20_m_immune, dpi=dpi,
            )
        else:
            logging.warning("CHC20 immune score empty; skipping Figure 5A.")

        if not chc20_stromal.empty:
            logging.info("Figure 5B: CHC20 stromal/ECM module spatial score...")
            plot_fig5B_stromal_spatial(
                chc20_adata, chc20_stromal, out_dir=out_dir,
                matched_genes=chc20_m_stromal, dpi=dpi,
            )
        else:
            logging.warning("CHC20 stromal score empty; skipping Figure 5B.")

        if not chc20_combined_score.empty:
            logging.info("Figure 5C: CHC20 combined immune-stromal niche score...")
            plot_fig5C_combined_niche_spatial(
                chc20_adata, chc20_combined_score, out_dir=out_dir, dpi=dpi,
            )
            # 保存 CHC20 双模块评分 CSV
            chc20_module_csv = out_dir / "fig5_dual_module_scores_CHC20.csv"
            pd.DataFrame({
                "spot_id":                     chc20_combined_score.index,
                "sample":                      "CHC20",
                "immune_signature_score":      chc20_immune.reindex(chc20_combined_score.index),
                "stromal_ECM_signature_score": chc20_stromal.reindex(chc20_combined_score.index),
                "immune_stromal_niche_score":  chc20_combined_score.to_numpy(dtype=float),
            }).to_csv(chc20_module_csv, index=False)
            logging.info("Saved CHC20 dual-module scores: %s", chc20_module_csv)

            # 保存使用的基因列表
            (out_dir / "fig5_immune_genes_used_CHC20.txt").write_text(
                "\n".join(chc20_m_immune), encoding="utf-8"
            )
            (out_dir / "fig5_stromal_genes_used_CHC20.txt").write_text(
                "\n".join(chc20_m_stromal), encoding="utf-8"
            )
        else:
            logging.warning("CHC20 combined score empty; skipping Figure 5C.")
    else:
        logging.warning(
            "CHC20 AnnData not available. To generate Figure 5A/5B/5C, provide:\n"
            "  --chc20-adata results/CHC20/adata_vis_post.h5ad"
        )

    # ── Fig.5D: HCC4R vs CHC20 combined niche score 对比小提琴图 ──
    if not hcc4r_combined_score.empty and not chc20_combined_score.empty:
        logging.info("Figure 5D: HCC4R vs CHC20 combined niche score comparison...")
        plot_fig5D_comparison_violin(
            hcc4r_combined_score, chc20_combined_score,
            out_dir=out_dir, dpi=dpi,
        )
    else:
        logging.warning(
            "HCC4R or CHC20 combined score empty; skipping Figure 5D.\n"
            "  HCC4R combined empty: %s | CHC20 combined empty: %s",
            hcc4r_combined_score.empty,
            chc20_combined_score.empty,
        )

    # ── 敏感性分析（旧版 Fig.5C 功能，保留为补充图）──────────────
    logging.info("Figure 5-Supplement: Sensitivity analysis...")
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
