"""
================================================================================
脚本名称: run_chc23_validation.py
功能概述: CHC23 切片独立验证分析 —— Step 1.5（介于 Step 1 与 Step 2 之间）
================================================================================

【整体任务说明】
    本脚本以 CHC23 切片的 Cell2location 反卷积结果为输入，执行三个阶段的
    独立验证分析，证明在 CHC20 上发现的免疫抑制空间生态位具有跨患者重现性：

    阶段一（Step 2 重跑）：
        对 CHC23 运行与主分析完全相同的 Niche 分析流程，输出平行的结果目录
        results/spatial_niche_chc23/，用于与 CHC20 逐图对比。

    阶段二（跨切片一致性量化）：
        1. Signature gene 重叠度（Jaccard 相似度）：
           读取两张切片各自 Top-50 signature gene 列表，计算 Jaccard 系数。
           验收标准：Jaccard ≥ 0.4（50 个基因中 ≥ 20 个重叠）。
        2. 细胞类型比例 Spearman 相关：
           对两个切片的各细胞类型均值比例向量做 Spearman 相关分析。
        3. niche_high 占比一致性：
           从 spatial_niche_parameters.csv 读取 niche_high 占比，
           判断两切片差异是否 < 5%。

    阶段三（Step 3 DE 验证）：
        使用 CHC23 自己的 niche_high 标签作为分组变量，调用 scanpy Wilcoxon
        检验完成差异表达分析，输出 CHC23 的 signature gene 列表。
        最终生成 CHC20/CHC23 Top-50 signature gene 重叠可视化图（UpSet 风格）。

【输入文件】
    results/adata_vis_post_CHC23.h5ad         - CHC23 Cell2location 反卷积后的空间数据
                                                （由 run_preprocessing.py Step 7 生成）
    results/spatial_niche/                    - CHC20 主分析 niche 结果目录
      ├── immunosuppressive_niche_signature_genes_ranked.csv  - CHC20 Top 基因表
      └── spatial_niche_parameters.csv                        - CHC20 分析参数

【输出文件】
    results/spatial_niche_chc23/              - CHC23 平行 Niche 分析结果目录
      ├── spatial_niche_scores.csv            - CHC23 完整评分表
      ├── spatial_niche_parameters.csv        - CHC23 分析参数
      ├── immunosuppressive_niche_signature_genes_ranked.csv
      ├── immunosuppressive_niche_signature_genes.txt
      └── plots/                              - CHC23 各类可视化图表
    results/chc23_validation/                 - 跨切片一致性量化结果
      ├── cross_slice_consistency_report.csv  - 一致性量化指标汇总表
      ├── signature_gene_overlap_venn.png     - 两切片 Top-50 基因重叠韦恩图
      ├── celltype_spearman_correlation.png   - 细胞类型比例 Spearman 相关图
      └── gene_overlap_heatmap.png            - 共同基因表达热图
    results/spatial_signature_genes_chc23.txt - CHC23 DE 验证的 signature gene 列表

【参考文献】
    - Schürch et al., Cell, 2020 (Cellular Neighborhoods 方法)
    - Kleshchevnikov et al., Nature Biotechnology, 2022 (Cell2location)
================================================================================
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats


# ============================================================
# 日志配置
# ============================================================

def _setup_logging() -> None:
    """初始化全局日志系统（格式：时间戳 | 级别 | 消息，输出到 stdout）。"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


# ============================================================
# 参数解析
# ============================================================

def _parse_args() -> argparse.Namespace:
    """
    解析命令行参数，为整个验证流程提供配置选项。

    所有参数均有默认值（基于项目标准目录结构），可直接运行无需传参。

    返回
    ----
    argparse.Namespace，包含所有参数值的命名空间对象。
    """
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    parser = argparse.ArgumentParser(
        description="CHC23 independent validation analysis (Step 1.5)."
    )
    # 输入路径
    parser.add_argument(
        "--chc23-adata", type=Path,
        default=results / "adata_vis_post_CHC23.h5ad",
        help="CHC23 Cell2location 反卷积后的 AnnData 文件路径",
    )
    parser.add_argument(
        "--chc20-niche-dir", type=Path,
        default=results / "spatial_niche",
        help="CHC20 主分析 niche 结果目录",
    )
    # 输出路径
    parser.add_argument(
        "--chc23-niche-dir", type=Path,
        default=results / "spatial_niche_chc23",
        help="CHC23 niche 分析输出目录",
    )
    parser.add_argument(
        "--validation-dir", type=Path,
        default=results / "chc23_validation",
        help="跨切片一致性量化结果输出目录",
    )
    parser.add_argument(
        "--chc23-sig-out", type=Path,
        default=results / "spatial_signature_genes_chc23.txt",
        help="CHC23 DE 验证 signature gene 列表输出路径",
    )
    # Niche 分析参数（与主分析保持一致）
    parser.add_argument("--abundance-key",  default="means_cell_abundance_w_sf")
    parser.add_argument("--hepatocyte-col", default="Hepatocyte")
    parser.add_argument("--treg-col",       default="Treg")
    parser.add_argument("--myeloid-col",    default="Myeloid")
    parser.add_argument("--fibroblast-col", default="Fibroblast")
    parser.add_argument("--tnk-col",        default="T/NK")
    parser.add_argument("--n-neighbors",    type=int,   default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--niche-high-quantile", type=float, default=0.80)
    parser.add_argument("--top-niche-genes", type=int, default=50)
    parser.add_argument("--top-n-overlap",   type=int, default=50,
                        help="用于 Jaccard 计算的 Top-N 基因数量")
    parser.add_argument("--jaccard-threshold", type=float, default=0.4,
                        help="Jaccard 相似度验收标准（默认 0.4）")
    parser.add_argument("--seed", type=int, default=1234)
    return parser.parse_args()


# ============================================================
# 阶段一：调用 run_spatial_niche_analysis.py 重跑 CHC23 Niche 分析
# ============================================================

def run_chc23_niche_analysis(
    chc23_adata: Path,
    chc23_niche_dir: Path,
    chc23_sig_out: Path,
    n_neighbors: int,
    leiden_resolution: float,
    niche_high_quantile: float,
    top_niche_genes: int,
    seed: int,
) -> int:
    """
    调用 run_spatial_niche_analysis.py 对 CHC23 执行与主分析完全相同的 Niche 分析。

    通过 subprocess 调用，确保环境一致性，同时将参数与主分析对齐（相同的
    n_neighbors、leiden_resolution、niche_high_quantile、seed），
    使两切片的结果具有可比性。

    参数
    ----
    chc23_adata        : CHC23 反卷积后的 AnnData 文件路径
    chc23_niche_dir    : CHC23 niche 结果输出目录
    chc23_sig_out      : CHC23 signature gene 输出路径
    n_neighbors        : kNN 邻居数（与主分析一致）
    leiden_resolution  : Leiden 聚类分辨率（与主分析一致）
    niche_high_quantile: niche_high 分位数阈值（与主分析一致）
    top_niche_genes    : 输出 signature gene 数量
    seed               : 随机数种子

    返回
    ----
    int，subprocess 返回码（0 = 成功）。
    """
    script = Path(__file__).resolve().parent / "run_spatial_niche_analysis.py"
    cmd = [
        sys.executable, str(script),
        "--adata",               str(chc23_adata),
        "--out-dir",             str(chc23_niche_dir),
        "--signature-out",       str(chc23_sig_out),
        "--n-neighbors",         str(n_neighbors),
        "--leiden-resolution",   str(leiden_resolution),
        "--niche-high-quantile", str(niche_high_quantile),
        "--top-niche-genes",     str(top_niche_genes),
        "--seed",                str(seed),
    ]
    logging.info("Phase 1: Running CHC23 niche analysis...")
    logging.info("Command: %s", " ".join(cmd))
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        logging.error(
            "CHC23 niche analysis failed with return code %d.", result.returncode
        )
    else:
        logging.info("Phase 1: CHC23 niche analysis completed successfully.")
    return result.returncode


# ============================================================
# 阶段二：跨切片一致性量化
# ============================================================

def _load_ranked_genes(ranked_csv: Path, top_n: int) -> list[str]:
    """
    从 signature gene 排名 CSV 文件中读取 Top-N 基因列表。

    CSV 文件由 _rank_niche_genes 生成，包含 gene / mean_high / mean_low /
    log2_fc 等列，已按 log2_fc 降序排列。

    参数
    ----
    ranked_csv : immunosuppressive_niche_signature_genes_ranked.csv 文件路径
    top_n      : 取前 top_n 个基因

    返回
    ----
    list[str]，基因名列表（最多 top_n 个）。
    """
    if not ranked_csv.exists():
        logging.warning("Ranked gene CSV not found: %s", ranked_csv)
        return []
    df = pd.read_csv(ranked_csv)
    if "gene" not in df.columns:
        logging.warning("Column 'gene' not found in %s", ranked_csv)
        return []
    return df["gene"].head(top_n).tolist()


def compute_gene_jaccard(genes_chc20: list[str], genes_chc23: list[str]) -> float:
    """
    计算两张切片 Top-N signature gene 列表之间的 Jaccard 相似度。

    Jaccard = |A ∩ B| / |A ∪ B|
    值域 [0, 1]，越接近 1 说明两切片发现的特征基因越一致。
    验收标准：≥ 0.4（即 50 个基因中至少 20 个重叠）认为一致性良好。

    参数
    ----
    genes_chc20 : CHC20 Top-N 基因列表
    genes_chc23 : CHC23 Top-N 基因列表

    返回
    ----
    float，Jaccard 系数（0.0 表示无重叠，1.0 表示完全重叠）。
    """
    set20 = set(genes_chc20)
    set23 = set(genes_chc23)
    intersection = set20 & set23
    union = set20 | set23
    if not union:
        return 0.0
    jaccard = len(intersection) / len(union)
    logging.info(
        "Signature gene Jaccard (Top-%d): %.3f  (overlap=%d, CHC20=%d, CHC23=%d)",
        max(len(genes_chc20), len(genes_chc23)),
        jaccard, len(intersection), len(set20), len(set23),
    )
    return jaccard


def compute_celltype_spearman(
    prop_chc20_csv: Path,
    prop_chc23_csv: Path,
) -> pd.DataFrame:
    """
    计算 CHC20 与 CHC23 各细胞类型均值比例向量之间的 Spearman 相关系数。

    每种细胞类型对应一个数值对（CHC20 均值，CHC23 均值），
    Spearman 相关衡量两切片在细胞组成模式上的相似程度。
    r > 0.8 认为组成模式高度一致。

    参数
    ----
    prop_chc20_csv : CHC20 spot_cell_proportion.csv 文件路径
    prop_chc23_csv : CHC23 spot_cell_proportion.csv 文件路径

    返回
    ----
    pd.DataFrame，包含 cell_type / mean_chc20 / mean_chc23 / spearman_r / p_value 列。
    """
    if not prop_chc20_csv.exists() or not prop_chc23_csv.exists():
        logging.warning("Cell proportion CSVs not found; skipping Spearman analysis.")
        return pd.DataFrame()

    prop20 = pd.read_csv(prop_chc20_csv)
    prop23 = pd.read_csv(prop_chc23_csv)

    # 排除 spot_id 列，只保留细胞类型比例列
    exclude = {"spot_id"}
    cols20 = [c for c in prop20.columns if c not in exclude]
    cols23 = [c for c in prop23.columns if c not in exclude]
    common_cols = [c for c in cols20 if c in cols23]

    if not common_cols:
        logging.warning("No common cell type columns found between CHC20 and CHC23.")
        return pd.DataFrame()

    mean20 = prop20[common_cols].mean()
    mean23 = prop23[common_cols].mean()

    rho, pval = stats.spearmanr(mean20.values, mean23.values)
    logging.info(
        "Cell type proportion Spearman r=%.3f (p=%.4f, n_celltypes=%d)",
        rho, pval, len(common_cols),
    )

    result_df = pd.DataFrame({
        "cell_type":  common_cols,
        "mean_chc20": mean20.values,
        "mean_chc23": mean23.values,
    })
    result_df.attrs["spearman_r"] = rho
    result_df.attrs["p_value"] = pval
    return result_df


def check_niche_high_consistency(
    params_chc20_csv: Path,
    params_chc23_csv: Path,
    tolerance: float = 0.05,
) -> dict:
    """
    比较两张切片的 niche_high spot 占比，评估生态位规模一致性。

    从各自的 spatial_niche_parameters.csv 中读取
    n_niche_high（niche_high spot 数量）和总 spot 数，
    计算各自的占比及差异。
    差异 < tolerance（默认 5%）认为规模一致。

    参数
    ----
    params_chc20_csv : CHC20 spatial_niche_parameters.csv 路径
    params_chc23_csv : CHC23 spatial_niche_parameters.csv 路径
    tolerance        : 允许的最大占比差异（默认 0.05 即 5%）

    返回
    ----
    dict，包含 pct_chc20、pct_chc23、diff、is_consistent 等字段。
    """
    result = {"pct_chc20": None, "pct_chc23": None, "diff": None, "is_consistent": None}

    def _read_param(csv_path: Path, key: str):
        """从参数 CSV 中按 parameter 列读取指定 key 的 value。"""
        if not csv_path.exists():
            return None
        df = pd.read_csv(csv_path)
        row = df[df["parameter"] == key]
        if row.empty:
            return None
        return float(row["value"].iloc[0])

    n_high_20  = _read_param(params_chc20_csv, "n_niche_high")
    n_high_23  = _read_param(params_chc23_csv, "n_niche_high")

    # 通过 n_niche_high / niche_high_quantile 反算总 spot 数（或直接读取）
    q20 = _read_param(params_chc20_csv, "niche_high_quantile") or 0.80
    q23 = _read_param(params_chc23_csv, "niche_high_quantile") or 0.80

    if n_high_20 is not None and q20 > 0:
        n_total_20 = n_high_20 / (1 - q20)
        pct_20 = n_high_20 / n_total_20
        result["pct_chc20"] = round(pct_20, 4)

    if n_high_23 is not None and q23 > 0:
        n_total_23 = n_high_23 / (1 - q23)
        pct_23 = n_high_23 / n_total_23
        result["pct_chc23"] = round(pct_23, 4)

    if result["pct_chc20"] is not None and result["pct_chc23"] is not None:
        diff = abs(result["pct_chc20"] - result["pct_chc23"])
        result["diff"] = round(diff, 4)
        result["is_consistent"] = diff < tolerance
        logging.info(
            "niche_high pct: CHC20=%.1f%%, CHC23=%.1f%%, diff=%.1f%% → %s",
            result["pct_chc20"] * 100,
            result["pct_chc23"] * 100,
            diff * 100,
            "PASS" if result["is_consistent"] else "FAIL",
        )

    return result


# ============================================================
# 阶段二：可视化
# ============================================================

def _plot_gene_overlap_venn(
    genes_chc20: list[str],
    genes_chc23: list[str],
    jaccard: float,
    path: Path,
) -> None:
    """
    绘制两切片 Top-N signature gene 重叠的韦恩图（简化版，使用 matplotlib 绘制圆圈）。

    图示：CHC20 专有基因（左）、共有基因（中间）、CHC23 专有基因（右），
    标注重叠基因名称（不超过 20 个以保持可读性）。

    参数
    ----
    genes_chc20 : CHC20 Top-N 基因列表
    genes_chc23 : CHC23 Top-N 基因列表
    jaccard     : 预计算的 Jaccard 系数（标注在图标题中）
    path        : 输出图片路径
    """
    set20 = set(genes_chc20)
    set23 = set(genes_chc23)
    only20 = sorted(set20 - set23)
    only23 = sorted(set23 - set20)
    common = sorted(set20 & set23)

    fig, ax = plt.subplots(figsize=(10, 6))
    # 绘制两个半透明圆圈模拟韦恩图
    circle20 = plt.Circle((0.35, 0.5), 0.28, color="#5b9bd5", alpha=0.35, transform=ax.transAxes, clip_on=False)
    circle23 = plt.Circle((0.65, 0.5), 0.28, color="#ed7d31", alpha=0.35, transform=ax.transAxes, clip_on=False)
    ax.add_patch(circle20)
    ax.add_patch(circle23)

    # 标注数量
    ax.text(0.22, 0.50, f"CHC20 only\n(n={len(only20)})",
            ha="center", va="center", fontsize=11, transform=ax.transAxes, fontweight="bold")
    ax.text(0.78, 0.50, f"CHC23 only\n(n={len(only23)})",
            ha="center", va="center", fontsize=11, transform=ax.transAxes, fontweight="bold")
    ax.text(0.50, 0.50, f"Shared\n(n={len(common)})",
            ha="center", va="center", fontsize=11, transform=ax.transAxes, fontweight="bold", color="#333333")

    # 标注共有基因名（最多 15 个）
    if common:
        gene_text = "\n".join(common[:15]) + ("..." if len(common) > 15 else "")
        ax.text(0.50, 0.10, gene_text, ha="center", va="bottom",
                fontsize=7, transform=ax.transAxes, color="#333333",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title(
        f"Signature Gene Overlap: CHC20 vs CHC23 (Top-{len(genes_chc20)})\n"
        f"Jaccard = {jaccard:.3f}",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Gene overlap Venn diagram saved: %s", path)


def _plot_celltype_spearman(
    celltype_df: pd.DataFrame,
    path: Path,
) -> None:
    """
    绘制细胞类型比例散点图（CHC20 均值 vs CHC23 均值），叠加 Spearman 相关系数。

    每个点代表一种细胞类型，x 轴为 CHC20 的均值比例，y 轴为 CHC23 的均值比例。
    对角线为完美一致参考线。Spearman r 越高，说明两切片组成模式越相似。

    参数
    ----
    celltype_df : 含 cell_type / mean_chc20 / mean_chc23 列的 DataFrame
                  （由 compute_celltype_spearman 返回，attrs 中含 spearman_r）
    path        : 输出图片路径
    """
    if celltype_df.empty:
        return

    rho = celltype_df.attrs.get("spearman_r", float("nan"))
    pval = celltype_df.attrs.get("p_value", float("nan"))

    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.scatter(celltype_df["mean_chc20"], celltype_df["mean_chc23"],
               s=60, color="#5b9bd5", edgecolors="white", linewidths=0.5, zorder=3)

    # 标注细胞类型名称
    for _, row in celltype_df.iterrows():
        ax.annotate(
            row["cell_type"],
            (row["mean_chc20"], row["mean_chc23"]),
            textcoords="offset points", xytext=(5, 3),
            fontsize=8, color="#444444",
        )

    # 对角线参考线
    lim_max = max(celltype_df["mean_chc20"].max(), celltype_df["mean_chc23"].max()) * 1.1
    ax.plot([0, lim_max], [0, lim_max], "k--", linewidth=0.8, alpha=0.5, label="y = x")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Mean proportion (CHC20)")
    ax.set_ylabel("Mean proportion (CHC23)")
    ax.set_title(
        f"Cell Type Proportion: CHC20 vs CHC23\n"
        f"Spearman r = {rho:.3f}  (p = {pval:.4f})",
        fontsize=11,
    )
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    logging.info("Cell type Spearman correlation plot saved: %s", path)


def _plot_common_gene_heatmap(
    adata_chc20: ad.AnnData,
    adata_chc23: ad.AnnData,
    common_genes: list[str],
    df_chc20: pd.DataFrame,
    df_chc23: pd.DataFrame,
    path: Path,
    max_genes: int = 30,
) -> None:
    """
    绘制两切片 niche_high vs niche_low 各共同基因均值表达量的热图。

    热图行为共同基因（最多 max_genes 个），列为四组：
    CHC20_niche_high、CHC20_niche_low、CHC23_niche_high、CHC23_niche_low。
    颜色编码 Z-score 标准化后的均值表达量，便于跨切片视觉对比。

    参数
    ----
    adata_chc20  : CHC20 AnnData（提供基因表达矩阵）
    adata_chc23  : CHC23 AnnData（提供基因表达矩阵）
    common_genes : 两切片共有 signature gene 列表
    df_chc20     : CHC20 spatial_niche_scores.csv 内容（含 niche_high 列）
    df_chc23     : CHC23 spatial_niche_scores.csv 内容（含 niche_high 列）
    path         : 输出图片路径
    max_genes    : 热图最大显示基因数（避免过于拥挤）
    """
    if not common_genes:
        logging.warning("No common genes to plot heatmap.")
        return

    use_genes = common_genes[:max_genes]

    def _mean_expr(adata: ad.AnnData, genes: list[str], mask: pd.Series) -> pd.Series:
        """计算给定 mask（bool Series）下各基因的均值表达量。"""
        import scipy.sparse as sp
        x = adata.X
        if sp.issparse(x):
            x = x.toarray()
        x = np.array(x, dtype=float)
        mask_arr = mask.reindex(adata.obs_names).fillna(False).astype(bool).to_numpy()
        gene_idx = [
            list(adata.var_names).index(g)
            for g in genes if g in adata.var_names
        ]
        present = [g for g in genes if g in adata.var_names]
        if not gene_idx:
            return pd.Series(dtype=float)
        subset = x[mask_arr, :][:, gene_idx]
        return pd.Series(subset.mean(axis=0), index=present)

    high_mask20 = df_chc20.set_index("spot_id")["niche_high"].astype(bool) if "spot_id" in df_chc20.columns \
        else df_chc20["niche_high"].astype(bool)
    high_mask23 = df_chc23.set_index("spot_id")["niche_high"].astype(bool) if "spot_id" in df_chc23.columns \
        else df_chc23["niche_high"].astype(bool)

    col_20h = _mean_expr(adata_chc20, use_genes, high_mask20)
    col_20l = _mean_expr(adata_chc20, use_genes, ~high_mask20)
    col_23h = _mean_expr(adata_chc23, use_genes, high_mask23)
    col_23l = _mean_expr(adata_chc23, use_genes, ~high_mask23)

    hm = pd.DataFrame({
        "CHC20_High": col_20h,
        "CHC20_Low":  col_20l,
        "CHC23_High": col_23h,
        "CHC23_Low":  col_23l,
    }).dropna()

    if hm.empty:
        logging.warning("Empty heatmap data; skipping.")
        return

    # Z-score 行归一化
    hm_z = hm.subtract(hm.mean(axis=1), axis=0).divide(
        hm.std(axis=1).replace(0, 1.0), axis=0
    )

    fig, ax = plt.subplots(figsize=(7, max(4, len(hm_z) * 0.35 + 1.5)))
    sns.heatmap(
        hm_z, cmap="RdBu_r", center=0,
        annot=False, linewidths=0.3,
        ax=ax, cbar_kws={"label": "Z-score"},
    )
    ax.set_title(
        f"Shared Signature Genes: CHC20 vs CHC23\n"
        f"(niche_high vs niche_low, top {len(hm_z)} genes)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Common gene expression heatmap saved: %s", path)


def run_cross_slice_consistency(
    chc20_niche_dir: Path,
    chc23_niche_dir: Path,
    results_dir: Path,
    validation_dir: Path,
    top_n: int,
    jaccard_threshold: float,
    chc20_adata_path: Path | None = None,
    chc23_adata_path: Path | None = None,
) -> dict:
    """
    执行阶段二：跨切片一致性量化分析，生成量化指标表和可视化图。

    按顺序完成：
      1. 计算 signature gene Jaccard 相似度，判断是否达到验收标准；
      2. 计算细胞类型比例 Spearman 相关，评估组成模式一致性；
      3. 检查 niche_high 占比差异是否 < 5%；
      4. 绘制韦恩图、散点图、基因表达热图；
      5. 汇总所有指标写入 cross_slice_consistency_report.csv。

    参数
    ----
    chc20_niche_dir  : CHC20 niche 结果目录（含 CSV 文件）
    chc23_niche_dir  : CHC23 niche 结果目录（含 CSV 文件）
    results_dir      : 项目 results 根目录（用于定位 proportion CSV）
    validation_dir   : 验证结果输出目录
    top_n            : 用于 Jaccard 计算的 Top-N 基因数
    jaccard_threshold: Jaccard 验收阈值
    chc20_adata_path : CHC20 AnnData 路径（用于热图，可选）
    chc23_adata_path : CHC23 AnnData 路径（用于热图，可选）

    返回
    ----
    dict，包含所有量化指标（jaccard、spearman_r、niche_pct_diff 等）。
    """
    validation_dir.mkdir(parents=True, exist_ok=True)
    logging.info("Phase 2: Cross-slice consistency analysis → %s", validation_dir)

    report = {}

    # ── 1. Signature gene Jaccard ──────────────────────────────────────────────
    genes_chc20 = _load_ranked_genes(
        chc20_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv", top_n
    )
    genes_chc23 = _load_ranked_genes(
        chc23_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv", top_n
    )
    jaccard = compute_gene_jaccard(genes_chc20, genes_chc23)
    report["jaccard"] = round(jaccard, 4)
    report["jaccard_pass"] = jaccard >= jaccard_threshold
    report["n_overlap"] = len(set(genes_chc20) & set(genes_chc23))

    _plot_gene_overlap_venn(
        genes_chc20, genes_chc23, jaccard,
        validation_dir / "signature_gene_overlap_venn.png",
    )

    # ── 2. 细胞类型比例 Spearman 相关 ─────────────────────────────────────────
    prop_csv_chc20 = results_dir / "spot_cell_proportion_CHC20.csv"
    # 兼容两种命名约定
    if not prop_csv_chc20.exists():
        prop_csv_chc20 = results_dir / "spot_cell_proportion.csv"
    prop_csv_chc23 = results_dir / "spot_cell_proportion_CHC23.csv"

    celltype_df = compute_celltype_spearman(prop_csv_chc20, prop_csv_chc23)
    if not celltype_df.empty:
        report["spearman_r"] = round(celltype_df.attrs.get("spearman_r", float("nan")), 4)
        report["spearman_p"] = round(celltype_df.attrs.get("p_value", float("nan")), 6)
        _plot_celltype_spearman(
            celltype_df,
            validation_dir / "celltype_spearman_correlation.png",
        )

    # ── 3. niche_high 占比一致性 ──────────────────────────────────────────────
    niche_pct_result = check_niche_high_consistency(
        chc20_niche_dir / "spatial_niche_parameters.csv",
        chc23_niche_dir / "spatial_niche_parameters.csv",
    )
    report.update({
        "pct_chc20": niche_pct_result.get("pct_chc20"),
        "pct_chc23": niche_pct_result.get("pct_chc23"),
        "pct_diff":  niche_pct_result.get("diff"),
        "pct_consistent": niche_pct_result.get("is_consistent"),
    })

    # ── 4. 共同基因表达热图（可选，需要 AnnData）─────────────────────────────
    common_genes = sorted(set(genes_chc20) & set(genes_chc23))
    if common_genes and chc20_adata_path and chc23_adata_path:
        try:
            adata20 = ad.read_h5ad(chc20_adata_path)
            adata23 = ad.read_h5ad(chc23_adata_path)
            scores20 = pd.read_csv(chc20_niche_dir / "spatial_niche_scores.csv")
            scores23 = pd.read_csv(chc23_niche_dir / "spatial_niche_scores.csv")
            _plot_common_gene_heatmap(
                adata20, adata23, common_genes,
                scores20, scores23,
                validation_dir / "gene_overlap_heatmap.png",
            )
        except Exception as exc:
            logging.warning("Could not generate gene overlap heatmap: %s", exc)

    # ── 5. 汇总报告 ───────────────────────────────────────────────────────────
    report_df = pd.DataFrame([{
        "metric": k, "value": str(v)
    } for k, v in report.items()])
    report_csv = validation_dir / "cross_slice_consistency_report.csv"
    report_df.to_csv(report_csv, index=False)
    logging.info("Consistency report saved: %s", report_csv)

    # 终端打印摘要
    logging.info(
        "\n========== Cross-Slice Consistency Summary ==========\n"
        "  Jaccard (Top-%d genes): %.3f → %s (threshold=%.2f)\n"
        "  Spearman r:             %.3f  (p=%.4f)\n"
        "  niche_high pct CHC20:   %.1f%%\n"
        "  niche_high pct CHC23:   %.1f%%\n"
        "  pct diff:               %.1f%% → %s\n"
        "=====================================================",
        top_n,
        report.get("jaccard", float("nan")),
        "PASS ✓" if report.get("jaccard_pass") else "FAIL ✗",
        jaccard_threshold,
        report.get("spearman_r", float("nan")),
        report.get("spearman_p", float("nan")),
        (report.get("pct_chc20") or 0) * 100,
        (report.get("pct_chc23") or 0) * 100,
        (report.get("pct_diff") or 0) * 100,
        "PASS ✓" if report.get("pct_consistent") else "FAIL ✗",
    )

    return report


# ============================================================
# 阶段三：CHC23 Step 3 DE 验证
# ============================================================

def run_chc23_de_validation(
    chc23_adata_path: Path,
    chc23_scores_csv: Path,
    chc23_sig_out: Path,
    top_n: int = 50,
) -> None:
    """
    对 CHC23 执行差异表达分析（Step 3 验证）。

    直接在 Python 内完成 scanpy Wilcoxon 检验，无需调用外部脚本。
    使用 CHC23 niche 分析输出的 niche_high 布尔列作为分组标签，
    在 niche_high=True 的 spot 中寻找特异性高表达基因。

    参数
    ----
    chc23_adata_path : CHC23 AnnData 文件路径
    chc23_scores_csv : CHC23 spatial_niche_scores.csv 文件路径（含 niche_high 列）
    chc23_sig_out    : 输出的 signature gene 列表文件路径
    top_n            : 导出 Top-N 个差异最显著的基因
    """
    if not chc23_adata_path.exists():
        logging.error("CHC23 AnnData not found: %s", chc23_adata_path)
        return
    if not chc23_scores_csv.exists():
        logging.error("CHC23 niche scores CSV not found: %s", chc23_scores_csv)
        logging.error("Please run Phase 1 (run_chc23_niche_analysis) first.")
        return

    logging.info("Phase 3: Running DE validation for CHC23...")

    # 加载数据
    adata = ad.read_h5ad(chc23_adata_path)
    scores_df = pd.read_csv(chc23_scores_csv)

    # 对齐 niche_high 标签
    if "spot_id" in scores_df.columns:
        scores_df = scores_df.set_index("spot_id")
    if "niche_high" not in scores_df.columns:
        logging.error("Column 'niche_high' not found in %s", chc23_scores_csv)
        return

    labels = scores_df["niche_high"].reindex(adata.obs_names)
    if labels.isna().any():
        logging.warning(
            "%d spots missing niche_high labels; filling with False.",
            int(labels.isna().sum()),
        )
        labels = labels.fillna(False)

    adata.obs["niche_high"] = labels.astype("category")

    # 准备表达矩阵（优先使用 log1p 层）
    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"].copy()
    else:
        logging.warning("Layer 'log1p' missing; normalizing current X.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Wilcoxon 差异表达
    cats = list(adata.obs["niche_high"].cat.categories)
    if len(cats) < 2:
        logging.warning("Only one niche_high category found; skipping DE.")
        return

    # True / "True" / 1 均视为阳性
    target = None
    for cand in (True, "True", "true", 1, "1"):
        if cand in cats:
            target = cand
            break
    if target is None:
        target = cats[-1]

    sc.tl.rank_genes_groups(adata, groupby="niche_high", method="wilcoxon")
    df_de = sc.get.rank_genes_groups_df(adata, group=str(target))
    genes = df_de["names"].head(top_n).tolist()

    chc23_sig_out.parent.mkdir(parents=True, exist_ok=True)
    chc23_sig_out.write_text("\n".join(genes) + "\n", encoding="utf-8")
    logging.info(
        "Phase 3: CHC23 DE signature genes (%d) saved: %s", len(genes), chc23_sig_out
    )


# ============================================================
# 最终报告：并排展示 CHC20 / CHC23 Top-50 基因重叠
# ============================================================

def _plot_gene_overlap_summary(
    genes_chc20: list[str],
    genes_chc23: list[str],
    path: Path,
) -> None:
    """
    绘制 CHC20 和 CHC23 Top-50 signature gene 的并排汇总表格图。

    左列：CHC20 基因（共有基因标红）
    右列：CHC23 基因（共有基因标红）
    用于论文中直观展示两切片发现的核心基因及其重叠情况。

    参数
    ----
    genes_chc20 : CHC20 Top-50 基因列表（按 log2FC 降序）
    genes_chc23 : CHC23 Top-50 基因列表（按 log2FC 降序）
    path        : 输出图片路径
    """
    common = set(genes_chc20) & set(genes_chc23)
    n_rows = max(len(genes_chc20), len(genes_chc23), 1)
    fig_h = min(max(6, n_rows * 0.28 + 1.5), 20)

    fig, axes = plt.subplots(1, 2, figsize=(10, fig_h))

    for ax, genes, title, color in zip(
        axes,
        [genes_chc20, genes_chc23],
        ["CHC20 Top Signature Genes", "CHC23 Top Signature Genes"],
        ["#5b9bd5", "#ed7d31"],
    ):
        ax.set_title(title, fontsize=11, fontweight="bold", color=color)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, n_rows + 1)
        ax.axis("off")
        for i, gene in enumerate(genes):
            row_y = n_rows - i
            is_common = gene in common
            ax.text(
                0.5, row_y, gene,
                ha="center", va="center",
                fontsize=8,
                color="#d62728" if is_common else "#333333",
                fontweight="bold" if is_common else "normal",
            )

    # 添加图例
    handles = [
        mpatches.Patch(color="#d62728", label=f"Shared genes (n={len(common)})"),
        mpatches.Patch(color="#333333", label="Slice-specific genes"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle(
        f"Cross-Slice Signature Gene Comparison\n"
        f"(Red = shared by both slices, n={len(common)} / {len(set(genes_chc20) | set(genes_chc23))} total)",
        fontsize=11, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Gene overlap summary figure saved: %s", path)


# ============================================================
# 主流程
# ============================================================

def main() -> int:
    """
    主函数：依次执行 CHC23 验证的三个阶段。

    Phase 1: 调用 run_spatial_niche_analysis.py 对 CHC23 重跑 Niche 分析
    Phase 2: 跨切片一致性量化（Jaccard、Spearman、niche_high 占比对比）
    Phase 3: CHC23 DE 验证（Wilcoxon 检验，输出 signature gene 列表）

    返回 0 表示正常完成。
    """
    _setup_logging()
    args = _parse_args()
    np.random.seed(args.seed)

    root = Path(__file__).resolve().parents[1]
    results = root / "results"

    # ── Phase 1: CHC23 Niche 分析 ────────────────────────────────────────────
    phase1_rc = run_chc23_niche_analysis(
        chc23_adata=args.chc23_adata,
        chc23_niche_dir=args.chc23_niche_dir,
        chc23_sig_out=args.chc23_sig_out,
        n_neighbors=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
        niche_high_quantile=args.niche_high_quantile,
        top_niche_genes=args.top_niche_genes,
        seed=args.seed,
    )
    if phase1_rc != 0:
        logging.warning(
            "Phase 1 returned non-zero exit code %d; proceeding to Phase 2 anyway.",
            phase1_rc,
        )

    # ── Phase 2: 跨切片一致性量化 ───────────────────────────────────────────
    consistency = run_cross_slice_consistency(
        chc20_niche_dir=args.chc20_niche_dir,
        chc23_niche_dir=args.chc23_niche_dir,
        results_dir=results,
        validation_dir=args.validation_dir,
        top_n=args.top_n_overlap,
        jaccard_threshold=args.jaccard_threshold,
        chc20_adata_path=results / "adata_vis_post.h5ad",
        chc23_adata_path=args.chc23_adata,
    )

    # 最终并排基因图
    genes_chc20 = _load_ranked_genes(
        args.chc20_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv",
        args.top_n_overlap,
    )
    genes_chc23 = _load_ranked_genes(
        args.chc23_niche_dir / "immunosuppressive_niche_signature_genes_ranked.csv",
        args.top_n_overlap,
    )
    if genes_chc20 or genes_chc23:
        _plot_gene_overlap_summary(
            genes_chc20, genes_chc23,
            args.validation_dir / "gene_overlap_parallel_list.png",
        )

    # ── Phase 3: CHC23 DE 验证 ───────────────────────────────────────────────
    run_chc23_de_validation(
        chc23_adata_path=args.chc23_adata,
        chc23_scores_csv=args.chc23_niche_dir / "spatial_niche_scores.csv",
        chc23_sig_out=args.chc23_sig_out,
        top_n=args.top_niche_genes,
    )

    # ── 最终日志 ─────────────────────────────────────────────────────────────
    jaccard = consistency.get("jaccard", float("nan"))
    logging.info(
        "\n===== CHC23 Validation Complete =====\n"
        "  Phase 1 (Niche analysis):     %s\n"
        "  Phase 2 (Jaccard):            %.3f → %s\n"
        "  Phase 3 (DE signature genes): %s\n"
        "  All outputs → %s\n"
        "=====================================",
        "OK" if phase1_rc == 0 else f"WARN (rc={phase1_rc})",
        jaccard,
        "PASS ✓" if consistency.get("jaccard_pass") else "FAIL ✗",
        args.chc23_sig_out,
        args.validation_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
