"""
================================================================================
脚本名称: run_spatial_niche_analysis.py
功能概述: 空间免疫抑制生态位（Niche）识别与特征基因提取 —— Step 2
================================================================================

【整体任务说明】
    本脚本基于 Cell2location 反卷积结果，采用"邻域组成聚类（无监督）+ 功能评分
    语义注释"的两阶段策略识别肝癌空间免疫抑制生态位（Schürch et al., 2020），
    完整流程共 15 个步骤：

      Step 1   加载反卷积空间 AnnData，提取细胞类型丰度矩阵；
      Step 2   构建 k-NN 空间邻接图（Squidpy），统计每个 spot 的邻域细胞组成向量；
      Step 3   对邻域组成向量进行 Leiden 无监督聚类，识别 cellular neighborhoods；
      Step 4   计算邻域统计量，汇总每个聚类簇的细胞类型均值及 spot 数；
      Step 5   构建半径邻域矩阵（radius_multiplier 控制范围），
               计算基于半径的邻域细胞组成向量（为敏感性分析提供第二种邻域定义）；
      Step 6   计算各 spot 到最近 hep_high spot 的距离（niche 空间边界辅助指标）；
      Step 7   计算免疫抑制功能基因模块评分（immunosuppressive_gene_score），
               结合 FOXP3 / IL10 / CTLA4 / TGFB1 / IDO1 等先验基因集；
      Step 8   计算 L-R 配体-受体通讯强度评分（基于 LR_PAIRS 预设对）；
      Step 9   多层次评分计算：
                 - Treg_like_score（Treg 比例 + 免疫抑制基因评分 Z-score 加总）；
                 - immune_stroma_score（Treg + T/NK + Myeloid + Fibroblast Z-score）；
                 - immunosuppressive_niche_score（综合五维 Z-score，用于分位数截断）；
                 - niche_high 标签（score ≥ niche_high_quantile 分位数阈值）；
      Step 10  辅助空间区域标注（规则化：tumor_core / stroma_immune / tumor_edge / other），
               仅用于可视化，不参与 niche 发现；
      Step 11  敏感性分析（v2 连续指标）：
                 对 kNN 邻居数（k=10/15/20）和 radius_multiplier（1.0/1.25/1.5）
                 两个维度分别做参数扰动，采用连续性指标衡量 niche 评分稳健性：
                   · Spearman ρ：衡量全局评分排序的一致性，不受硬截断影响；
                   · 加权 Jaccard（Continuous Jaccard）：对高分区域天然加权，
                     ∑min(a,b) / ∑max(a,b)，基于 Min-Max 归一化连续向量计算；
                 以 k=15 为参考配置，输出 2×2 矩阵图展示两指标在两个参数维度的稳定性；
      Step 12  保存完整评分表（spatial_niche_scores.csv）和分析参数元数据；
      Step 13  可视化（共 15+ 张图）：
                 细胞类型空间分布、综合 niche 评分空间图、Leiden 聚类图、
                 语义 niche 标签图、辅助区域标注图、距离依赖折线图、
                 区域评分箱线图、Hepatocyte-Treg 散点图、细胞类型相关性热图、
                 L-R 通讯热图、敏感性分析稳定性图、niche_high 二值分布图、
                 聚类法 vs 评分法并排对比图；
      Step 14  特征基因提取（三层筛选策略 + 检出率双维度评分）：
                 【问题背景】Visium spot 覆盖 5-50 个细胞，FOXP3 等 Treg 标志基因
                 因细胞稀释效应在大多数 spot 中表达为 0，导致均值 log2FC 接近 0，
                 纯 Wilcoxon+FC 策略无法将目标基因选入 Top-N；
                 【解决方案】引入检出率（Fraction of spots）作为补充维度：
                   · frac_high / frac_low：各组中检出该基因（表达>0）的 spot 比例；
                   · delta_frac = frac_high - frac_low：检出率差值；
                   · composite_score = 0.5×norm(log2FC) + 0.5×norm(delta_frac)：综合排序；
                 Layer 1 - Wilcoxon + BH-FDR + 综合排序：
                   筛选条件：FDR<0.05 且（log2FC>0.5 OR delta_frac>0.10）；
                   输出新列：frac_high / frac_low / delta_frac / composite_score；
                 Layer 2 - 先验功能基因集 AUC 检验（Treg/TAM/CAF 三组）：
                   增加 delta_frac 列；AUC>0.6 且 delta_frac>0.05 的先验基因强制合并；
                  Layer 3 - Gini Index 特异性评分（局灶性高表达稀有免疫基因检测）：
                    Gini>0.3 且 log2FC>0（阈值由 0.5 降至 0.3，覆盖 FOXP3 等稀有基因）；
                 绘制火山图（FC vs FDR）+ 检出率散点图（FC vs delta_frac）；
      Step 15  参数扫描稳定性评估（k × niche_high_quantile 二维网格搜索）：
               【修复说明】原 resolution × quantile 扫描存在根本缺陷——Leiden
               resolution 仅决定聚类粒度，不影响 niche_score 数值，导致所有列
               DEG 结果完全相同，热图无意义。改为 k × quantile 扫描，k（kNN
               邻居数）直接影响邻域组成向量，进而影响 niche_score 和 DEG 结果：
                 · 横轴：niche_high_quantile（0.70 / 0.75 / 0.80 / 0.85）
                 · 纵轴：k（8 / 10 / 15 / 20 / 25）
                 · 颜色：FDR<0.05 且 log2FC>0.5 的 DEG 数量
               选参标准：颜色最深且处于"高原"区域（与相邻格子结果相近）的
               参数组合为推荐；若主流程参数（k=15, q=0.80）位于高原区，则合理。

【niche 识别策略说明（参考 docs/niche修改.md）】
    - niche 发现主体：空间邻接图 → 邻域组成向量 → Leiden 无监督聚类；
    - 阈值规则（hep_high_quantile、niche_high_quantile）仅用于语义注释和签名提取分组；
    - 双层验证：Leiden 聚类法（无监督）与评分阈值法（规则化）并排可视化对比。

【输入文件】
    results/adata_vis_post.h5ad     - CHC20 Cell2location 反卷积后的空间 AnnData
                                      （由 run_preprocessing.py Step 6 生成）

【输出文件】
    results/spatial_niche/
      spatial_niche_scores.csv                           - 完整 spot 级别评分表
      spatial_niche_parameters.csv                       - 分析参数与阈值元数据
      sensitivity_analysis.csv                           - 敏感性分析结果
                                                           （含 spearman_rho / weighted_jaccard 两指标）
      neighborhood_cluster_stats.csv                     - 邻域聚类簇统计信息
      immunosuppressive_niche_signature_genes_ranked.csv - Layer 1 签名基因排名表
                                                           列：gene / mean_high / mean_low /
                                                               log2_fc / frac_high / frac_low /
                                                               delta_frac / composite_score /
                                                               pvalue / fdr
      immunosuppressive_niche_signature_genes.txt        - 签名基因列表（TCGA 投影接口）
      prior_gene_set_auc.csv                             - Layer 2 先验基因集 AUC 检验结果
                                                           （含 frac_high / frac_low / delta_frac 列）
      gini_score_genes.csv                               - Layer 3 Gini Index 特异性基因
                                                           （Gini>0.3 且 log2FC>0，
                                                            降低阈值以覆盖 FOXP3 等稀有免疫基因）
      param_scan_deg_stability.csv                       - 参数扫描 DEG 稳定性数据
                                                           列：k / quantile / n_sig_deg /
                                                               mean_log2fc_topN / top_genes_str
      plots/
        spatial_hepatocyte.png                           - Hepatocyte 空间分布图
        spatial_treg.png                                 - Treg 空间分布图
        spatial_myeloid.png                              - Myeloid 空间分布图
        spatial_fibroblast.png                           - Fibroblast 空间分布图
        spatial_immunosuppressive_niche_score.png        - 综合 niche 评分空间图
        spatial_neighborhood_clusters.png                - Leiden 邻域聚类空间图
        spatial_niche_semantic_labels.png                - 语义 niche 标签空间图
        spatial_region_labels.png                        - 辅助区域标注空间图
        spatial_niche_high_score_spots.png               - niche_high 二值分布图
        spatial_niche_cluster_vs_score_comparison.png   - 聚类法 vs 评分法对比图
        distance_to_hep_high_vs_niche_score.png         - 距离-niche评分折线图
        region_score_boxplots.png                        - 各区域评分箱线图
        hepatocyte_vs_treg_niche_score.png               - Hepatocyte-Treg 散点图
        celltype_niche_correlation.png                   - 细胞类型相关性热图
        lr_communication_heatmap.png                     - L-R 配体受体通讯热图
        sensitivity_niche_stability.png                  - 敏感性分析 2×2 指标图
        niche_signature_volcano.png                      - 全基因火山图（log2FC vs -log10 FDR）
        niche_fraction_scatter.png                       - 检出率差值散点图
                                                           （delta_frac vs log2FC；突出 FOXP3 等
                                                            稀有免疫基因的 niche 富集，补充火山图）
        param_scan_deg_stability_heatmap.png             - 参数扫描稳定性热图
    results/spatial_signature_genes.txt                  - 签名基因（TCGA 投影用途）

【依赖关系】
    上游：run_preprocessing.py（生成 adata_vis_post.h5ad）
    下游：run_chc23_validation.py（读取 spatial_niche/ 进行跨切片验证）
          run_de_analysis.py（读取 spatial_niche_scores.csv 进行 DE 分析）
          tcga_survival_analysis.R（读取 spatial_signature_genes.txt 进行 ssGSEA 预后分析）

【参考文献】
    - Schürch et al., Cell, 2020 (cellular neighborhoods 核心方法)
    - Palla et al., Nature Methods, 2022 (Squidpy 空间分析框架)
    - Keren et al., Cell, 2018 (肿瘤-免疫空间结构分型)
    - Spearman, 1904 (秩相关系数)
    - Jaccard, 1912 (Jaccard 相似度；本脚本采用连续加权版本)
================================================================================
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import cast

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.sparse as sp
import scipy.stats as scipy_stats
from scipy.spatial import KDTree
from statsmodels.stats.multitest import multipletests


# ============================================================
# 全局常量：免疫抑制相关基因（标志基因集）
#
# 选择依据（均有文献支持）：
#   FOXP3 / IL2RA (CD25)  — Treg 核心标志（Sakaguchi et al., 2020）
#   CTLA4 / TIGIT / LAG3 / PDCD1 (PD-1) / HAVCR2 (TIM-3)
#                          — 抑制性免疫检查点分子
#   TGFB1 / IL10           — 免疫抑制细胞因子
#   CXCL12 / CCL22         — 趋化因子（L-R 通讯关键分子）
#   TNFRSF18 (GITR) / TNFRSF4 (OX40) / IKZF2 (Helios)
#                          — Treg 激活与稳定相关因子
# ============================================================
IMMUNOSUPPRESSIVE_GENES: tuple[str, ...] = (
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "TIGIT",
    "LAG3",
    "PDCD1",
    "HAVCR2",
    "TGFB1",
    "IL10",
    "CXCL12",
    "CCL22",
    "TNFRSF18",
    "TNFRSF4",
    "IKZF2",
)

# ============================================================
# 全局常量：配体-受体通讯对（Treg-TAM-CAF 轴关键信号）
#
# 格式：(ligand, receptor, display_label)
# 参考文献：Efremova et al. (CellPhoneDB, 2020); Cang & Nie (COMMOT, 2023)
#
# 补充说明：
#   原有通路（CCL22-CCR4、IL10-IL10RA）因 Visium spot 级别的细胞稀释效应
#   在空间数据中表达极低，导致 product score 接近 0。
#   新增 TAM 相关高表达通路（SPP1-CD44、MIF-CD74 等），这类基因在 bulk spot
#   中信号更强，适合空间转录组的通讯分析。
# ============================================================
LR_PAIRS: list[tuple[str, str, str]] = [
    # 原有关键免疫轴（在 Visium 数据中可检测）
    ("CXCL12", "CXCR4",   "CXCL12–CXCR4"),
    ("TGFB1",  "TGFBR1",  "TGFB1–TGFBR1"),
    ("PDCD1",  "CD274",   "PD-1–PD-L1"),
    ("TIGIT",  "NECTIN2", "TIGIT–NECTIN2"),
    ("LAG3",   "HLA-DRA", "LAG3–MHC-II"),
    # 稀疏通路（保留，fallback 机制将尝试替代基因）
    ("CCL22",  "CCR4",    "CCL22–CCR4"),
    ("IL10",   "IL10RA",  "IL10–IL10RA"),
    # 新增：TAM/CAF 相关高表达通路（Visium 级别可检测）
    ("SPP1",   "CD44",    "SPP1–CD44"),       # 骨桥蛋白-CD44，TAM 分泌
    ("MIF",    "CD74",    "MIF–CD74"),         # 巨噬细胞迁移抑制因子
    ("VEGFA",  "KDR",     "VEGFA–KDR"),        # 血管生成，与 CAF 相关
    ("LGALS9", "HAVCR2",  "Galectin9–TIM-3"), # 抑制性检查点，可检测
    ("CCL2",   "CCR2",    "CCL2–CCR2"),        # 巨噬细胞招募通路
]

# ============================================================
# 全局常量：基因 Fallback 字典
#
# 当主基因在数据集中缺失时，自动尝试同家族替代基因。
# 适用于 L-R 通讯分析及先验基因集检验。
# ============================================================
GENE_FALLBACKS: dict[str, list[str]] = {
    "CCR4":   ["CCR2", "CCR5"],       # 趋化因子受体家族
    "IL10RA": ["IL10RB"],             # IL-10 受体亚基
    "NECTIN2": ["PVRL2", "CD112"],    # Nectin 家族别名
    "HAVCR2": ["TIM3", "TIMD4"],      # TIM-3 别名
    "HLA-DRA": ["HLA-DRB1", "CD74"], # MHC-II 相关
    "CXCR4":  ["CXCR7", "ACKR3"],    # CXCR4 替代受体
    "TGFBR1": ["TGFBR2"],            # TGF-β 受体
    "KDR":    ["FLT1", "FLT4"],      # VEGF 受体家族
}

# ============================================================
# 全局常量：先验功能基因集（用于 Mann-Whitney AUC 检验）
#
# 按细胞类型分三组，直接针对免疫抑制 niche 的核心基因进行
# 统计检验，弥补纯 log2FC 筛选无法发现低表达稀有基因的缺陷。
# ============================================================
PRIOR_GENE_SETS: dict[str, list[str]] = {
    "Treg_markers": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
    "TAM_features": ["CD163", "MRC1", "TGFB1", "IL10", "CXCL12"],
    "CAF_activation": ["FAP", "ACTA2", "POSTN", "COL1A1", "CCL22"],
}


# ============================================================
# 参数解析
# ============================================================

def _parse_args() -> argparse.Namespace:
    """
    解析命令行参数，为整个分析流程提供所有配置选项。

    所有参数均有合理默认值，可直接运行而无需传参。
    """
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    parser = argparse.ArgumentParser(
        description=(
            "Immunosuppressive spatial niche analysis using neighborhood "
            "composition clustering (Schürch et al., Cell 2020)."
        )
    )
    # --- 输入/输出路径 ---
    parser.add_argument(
        "--adata", type=Path,
        default=results / "adata_vis_post.h5ad",
        help="输入 AnnData .h5ad 文件路径；默认使用 CHC20 主分析结果",
    )
    parser.add_argument(
        "--out-dir", type=Path,
        default=results / "spatial_niche",
        help="结果输出目录",
    )
    parser.add_argument(
        "--signature-out", type=Path,
        default=results / "spatial_signature_genes.txt",
        help="特征基因列表输出路径（供 TCGA 分析使用）",
    )
    # --- 细胞类型列名 ---
    parser.add_argument("--abundance-key",  default="means_cell_abundance_w_sf")
    parser.add_argument("--hepatocyte-col", default="Hepatocyte")
    parser.add_argument("--treg-col",       default="Treg")
    parser.add_argument("--myeloid-col",    default="Myeloid")
    parser.add_argument("--fibroblast-col", default="Fibroblast")
    parser.add_argument("--tnk-col",        default="T/NK")
    # --- 邻域聚类参数 ---
    parser.add_argument(
        "--n-neighbors", type=int, default=15,
        help="构建邻域组成 k-NN 图时的邻居数（主分析 k=15）",
    )
    parser.add_argument(
        "--leiden-resolution", type=float, default=0.5,
        help="Leiden 聚类分辨率（越大 cluster 越多）",
    )
    # --- 阈值参数（仅用于辅助注释和签名提取，不用于 niche 发现主体） ---
    parser.add_argument(
        "--hep-high-quantile", type=float, default=0.75,
        help="高肝细胞区域分位数阈值（辅助空间区域标注用）",
    )
    parser.add_argument(
        "--niche-high-quantile", type=float, default=0.80,
        help="高免疫抑制 niche 分位数阈值（用于签名提取分组）",
    )
    parser.add_argument(
        "--neighbor-radius-multiplier", type=float, default=1.25,
        help="半径邻域搜索倍增系数（主分析参数，补充敏感性分析）",
    )
    # --- 其他参数 ---
    parser.add_argument(
        "--top-niche-genes", type=int, default=80,
        help="输出 niche 特征基因数量",
    )
    parser.add_argument("--seed", type=int, default=1234, help="随机数种子")
    return parser.parse_args()


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
# 数据读取与预处理工具函数
# ============================================================

def _clean_abundance_columns(columns: pd.Index) -> list[str]:
    """
    去除 Cell2location 后验矩阵列名中的统计量前缀，只保留细胞类型名称。

    示例：means_cell_abundance_w_sf_Hepatocyte → Hepatocyte
    """
    import re
    result = []
    for col in columns:
        v = str(col)
        v = re.sub(r"^(means|q\d+|median)_?cell_abundance_w_sf_", "", v)
        v = re.sub(r"^means_per_cluster_mu_fg_", "", v)
        result.append(v)
    return result


def _get_abundance(adata: ad.AnnData, key: str) -> pd.DataFrame:
    """
    从 AnnData 对象中提取细胞类型丰度矩阵，并清洗列名。

    参数
    ----
    adata : AnnData 空间转录组数据对象
    key   : adata.obsm 中细胞丰度矩阵的键名

    返回
    ----
    pd.DataFrame，行=spot，列=细胞类型名称。
    """
    if key not in adata.obsm:
        raise KeyError(
            f"Missing adata.obsm['{key}']; available: {list(adata.obsm.keys())}"
        )
    raw = adata.obsm[key]
    if isinstance(raw, pd.DataFrame):
        df = raw.copy()
    else:
        factors = adata.uns.get("mod", {}).get("factor_names")
        df = pd.DataFrame(raw, index=adata.obs_names, columns=factors)
    df.index = adata.obs_names
    df.columns = _clean_abundance_columns(df.columns)
    return df


def _zscore(series: pd.Series) -> pd.Series:
    """
    Z-score 标准化：z = (x - μ) / σ。

    标准差为零或无效时返回全零序列，保证数值稳定性。
    """
    std = float(series.std())
    if not np.isfinite(std) or std == 0:
        return pd.Series(0.0, index=series.index)
    return (series - float(series.mean())) / std


def _expression_frame(adata: ad.AnnData) -> pd.DataFrame:
    """
    将 AnnData 原始 count 矩阵转换为 CP10K + log1p 标准化的 DataFrame。

    流程：
      1. CP10K：count / total_count × 10000（消除测序深度差异）
      2. log1p：log(x + 1)（压缩数值范围，趋近正态分布）

    自动处理稀疏矩阵（scipy.sparse）和稠密矩阵（numpy.ndarray）。
    """
    x = adata.X
    if x is None:
        raise ValueError("adata.X is empty; cannot normalize expression matrix.")
    if sp.issparse(x):
        x = cast(sp.spmatrix, x).copy().astype(float).tocsr()
        totals = np.asarray(x.sum(axis=1)).ravel()
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        x = sp.diags(scale) @ x
        x.data = np.log1p(x.data)
        x = x.toarray()
    else:
        x = np.array(x, dtype=float, copy=True)
        totals = x.sum(axis=1)
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        x = np.log1p(x * scale[:, None])
    return pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names)


def _module_score(
    adata: ad.AnnData,
    genes: tuple[str, ...],
) -> tuple[pd.Series, list[str]]:
    """
    计算每个 spot 的免疫抑制基因模块评分（目标基因集的平均标准化表达量）。

    返回
    ----
    (score, available_genes)：
      score           - pd.Series，每个 spot 的基因模块评分；
      available_genes - 实际使用的基因列表（过滤不在数据集中的基因后）。
    """
    available = [g for g in genes if g in adata.var_names]
    if not available:
        logging.warning("No immunosuppressive marker genes found in adata.var_names.")
        return pd.Series(0.0, index=adata.obs_names, name="immunosuppressive_gene_score"), []
    score = _expression_frame(adata)[available].mean(axis=1)
    score.name = "immunosuppressive_gene_score"
    return score, available


# ============================================================
# 空间邻域工具函数
# ============================================================

def _spatial_neighbors(
    coords: np.ndarray,
    radius_multiplier: float,
) -> tuple[list[np.ndarray], float]:
    """
    基于固定半径构建每个 spot 的空间邻居列表。

    搜索半径 = 所有 spot 最近邻距离的中位数 × radius_multiplier。
    该自适应策略确保每个 spot 至少有一个邻居，同时避免纳入过远的无关 spot
    （参考 Schürch et al., 2020 的空间邻域分析方法）。

    参数
    ----
    coords            : (n_spots, 2) 空间坐标数组
    radius_multiplier : 搜索半径倍增系数（主分析 1.25）

    返回
    ----
    (邻居索引列表, 搜索半径数值)
    """
    tree = KDTree(coords)
    nearest = tree.query(coords, k=2)[0][:, 1]
    radius = float(np.median(nearest) * radius_multiplier)
    raw_lists = tree.query_ball_point(coords, r=radius)
    cleaned = [
        np.array([j for j in idxs if j != i], dtype=int)
        for i, idxs in enumerate(raw_lists)
    ]
    return cleaned, radius


def _build_knn_neighbors(coords: np.ndarray, k: int) -> list[np.ndarray]:
    """
    构建 k 最近邻（kNN）邻居列表，用于邻域组成聚类。

    排除自身，返回每个 spot 的 k 个最近邻索引。

    参数
    ----
    coords : (n_spots, 2) 空间坐标数组
    k      : 每个 spot 的邻居数量
    """
    tree = KDTree(coords)
    _, indices = tree.query(coords, k=k + 1)   # 第 0 列为自身，需排除
    return [indices[i, 1:] for i in range(len(coords))]


def _neighbor_mean(values: pd.Series, neighbors: list[np.ndarray]) -> pd.Series:
    """
    计算每个 spot 在其空间邻居范围内指定值的均值。

    无邻居的孤立 spot 输出 NaN。
    """
    arr = values.to_numpy()
    out = np.full(len(values), np.nan, dtype=float)
    for i, idx in enumerate(neighbors):
        if len(idx):
            out[i] = float(np.mean(arr[idx]))
    return pd.Series(out, index=values.index)


def _distance_to_mask(coords: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """
    计算每个 spot 到目标区域（mask=True 的 spot 集合）的最短欧式距离。

    目标区域为空时返回全 NaN 数组。
    """
    if not bool(mask.any()):
        return np.full(coords.shape[0], np.nan)
    tree = KDTree(coords[mask])
    return tree.query(coords, k=1)[0]


# ============================================================
# 核心：邻域组成聚类（Cellular Neighborhood Discovery）
# ============================================================

def _compute_neighborhood_composition(
    proportions: pd.DataFrame,
    neighbors: list[np.ndarray],
) -> pd.DataFrame:
    """
    计算每个 spot 的邻域细胞组成向量（局部邻居细胞类型比例的均值）。

    这是 Schürch et al. (2020) "Cellular Neighborhoods" 方法的核心步骤：
    不依赖先验阈值，而是让局部细胞组成自然形成 neighborhood 类型。

    参数
    ----
    proportions : (n_spots, n_celltypes) 细胞类型比例矩阵
    neighbors   : 每个 spot 的邻居索引列表（kNN 或半径邻域）

    返回
    ----
    pd.DataFrame，与 proportions 同形状，值为邻域均值。
    孤立 spot（无邻居）退回到自身组成。
    """
    arr = proportions.to_numpy()
    out = np.zeros_like(arr, dtype=float)
    for i, idx in enumerate(neighbors):
        out[i] = arr[idx].mean(axis=0) if len(idx) else arr[i]
    return pd.DataFrame(out, index=proportions.index, columns=proportions.columns)


def _leiden_cluster_neighborhood(
    neighborhood_comp: pd.DataFrame,
    n_neighbors: int,
    resolution: float,
    seed: int,
) -> pd.Series:
    """
    对邻域组成向量进行 Leiden 无监督聚类，识别 cellular neighborhoods。

    步骤：
      1. 以邻域组成向量（细胞类型特征空间）构建 k-NN 图；
      2. 运行 Leiden 社区检测；
      3. 返回每个 spot 的 cluster 标签。

    若 scanpy 不可用，自动退回到 KMeans（6 个聚类）。

    参数
    ----
    neighborhood_comp : (n_spots, n_celltypes) 邻域组成特征矩阵
    n_neighbors       : 特征 k-NN 图的邻居数
    resolution        : Leiden 分辨率参数
    seed              : 随机数种子

    返回
    ----
    pd.Series，cluster 标签字符串，索引与 neighborhood_comp 一致。
    """
    try:
        import scanpy as sc
    except ImportError:
        logging.warning("scanpy not available; falling back to KMeans clustering.")
        return _kmeans_cluster_neighborhood(neighborhood_comp, n_clusters=6, seed=seed)

    adata_nc = ad.AnnData(neighborhood_comp.to_numpy().astype(float))
    adata_nc.obs_names = neighborhood_comp.index.tolist()
    adata_nc.var_names = neighborhood_comp.columns.tolist()
    # 特征维度通常 < 20（细胞类型数），无需 PCA 降维，直接在原始特征空间建图
    sc.pp.neighbors(adata_nc, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata_nc, resolution=resolution, key_added="nc_cluster",
                 random_state=seed)
    return pd.Series(
        adata_nc.obs["nc_cluster"].tolist(),
        index=neighborhood_comp.index,
        name="neighborhood_cluster",
    )


def _kmeans_cluster_neighborhood(
    neighborhood_comp: pd.DataFrame,
    n_clusters: int = 6,
    seed: int = 1234,
) -> pd.Series:
    """
    备用聚类方案：KMeans（当 scanpy 不可用时使用）。
    """
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    labels = km.fit_predict(neighborhood_comp.to_numpy())
    return pd.Series(
        [str(x) for x in labels],
        index=neighborhood_comp.index,
        name="neighborhood_cluster",
    )


def _annotate_neighborhood_clusters(
    df: pd.DataFrame,
    cluster_col: str,
    treg_col: str,
    myeloid_col: str,
    fibroblast_col: str,
    gene_score_col: str = "immunosuppressive_gene_score",
) -> tuple[pd.Series, pd.DataFrame]:
    """
    基于功能评分对邻域 cluster 进行语义注释，识别免疫抑制 niche cluster。

    注释逻辑（"组成为主、功能为辅"）：
      1. 计算每个 cluster 四个维度的均值：
         Treg 比例、Myeloid 比例、Fibroblast 比例、免疫抑制基因评分；
      2. 对每个维度进行降序排名（rank=1 表示该维度最高）；
      3. 综合排名 = 四个维度 rank 之和（越小说明综合免疫抑制特征越强）；
      4. 综合排名最低（最小）的 cluster 注释为 "immunosuppressive_niche"；
      5. 其余 cluster 保留 "cluster_X" 格式标签。

    参数
    ----
    df            : 包含细胞比例、评分和聚类标签的 DataFrame
    cluster_col   : 邻域 cluster 标签列名
    treg_col      : Treg 比例列名
    myeloid_col   : Myeloid 比例列名
    fibroblast_col: Fibroblast 比例列名
    gene_score_col: 免疫抑制基因评分列名

    返回
    ----
    (niche_labels, cluster_stats)：
      niche_labels  - pd.Series，每个 spot 的语义标签；
      cluster_stats - pd.DataFrame，每个 cluster 的功能评分统计。
    """
    stats = (
        df.groupby(cluster_col)[[treg_col, myeloid_col, fibroblast_col, gene_score_col]]
        .mean()
        .rename(columns={
            treg_col:       "mean_treg",
            myeloid_col:    "mean_myeloid",
            fibroblast_col: "mean_fibroblast",
            gene_score_col: "mean_immune_score",
        })
    )
    for col in stats.columns:
        stats[f"rank_{col}"] = stats[col].rank(ascending=False)
    stats["combined_rank"] = (
        stats["rank_mean_treg"]
        + stats["rank_mean_myeloid"]
        + stats["rank_mean_fibroblast"]
        + stats["rank_mean_immune_score"]
    )

    niche_cluster_id = stats["combined_rank"].idxmin()
    logging.info(
        "Annotated immunosuppressive niche → cluster %s  "
        "(combined_rank=%.2f, mean_Treg=%.4f, mean_Myeloid=%.4f, "
        "mean_Fibroblast=%.4f, mean_ImmScore=%.4f)",
        niche_cluster_id,
        stats.loc[niche_cluster_id, "combined_rank"],
        stats.loc[niche_cluster_id, "mean_treg"],
        stats.loc[niche_cluster_id, "mean_myeloid"],
        stats.loc[niche_cluster_id, "mean_fibroblast"],
        stats.loc[niche_cluster_id, "mean_immune_score"],
    )

    niche_labels = df[cluster_col].apply(
        lambda c: "immunosuppressive_niche" if str(c) == str(niche_cluster_id)
        else f"cluster_{c}"
    )
    return niche_labels, stats


# ============================================================
# 敏感性分析
# ============================================================

def _sensitivity_analysis(
    proportions: pd.DataFrame,
    coords: np.ndarray,
    gene_score: pd.Series,
    treg_col: str,
    myeloid_col: str,
    fibroblast_col: str,
    n_neighbors_list: tuple[int, ...] = (10, 15, 20),
    radius_multiplier_list: tuple[float, ...] = (1.0, 1.25, 1.5),
    niche_high_quantile: float = 0.80,
    ref_k: int = 15,
) -> pd.DataFrame:
    """
    对邻域参数进行敏感性分析，评估 niche 评分分布的连续稳定性。

    【改进说明 v2 — 2026-06-26】
    v1 版本改用 Jaccard 相似度（硬截断集合），存在三个根本缺陷：
      1. Magnitude Loss：只判断是否过 80% 阈值，丢失评分绝对量级信息；
      2. 固定边缘概率：|A|=|B|=n×q 恒成立，分母被死死锁住，Jaccard 变动
         范围极窄，不反映真实差异；
      3. 边界极度敏感：阈值边缘处极微小的分数扰动即可导致 Jaccard 大幅跳变，
         而该变动并非空间模式的实质性变化（硬截断噪声）。

    v2 版本改用两个互补的连续型指标，彻底避免硬截断：

      指标一：Spearman 秩相关系数（spearman_rho）
        - 直接计算参考评分序列与实验评分序列之间所有 spot 排名的相关性；
        - 保留全部 spot 信息，不受任何截断影响；
        - 衡量「参数变化是否系统性地改变了 spot 的相对排序」；
        - 值域 [−1, 1]，越接近 1.0 说明排序越稳定。

      指标二：连续型 Weighted Jaccard（weighted_jaccard）
        - 先将每个评分向量归一化到 [0, 1]（min-max）；
        - 计算 Σ min(a_i, b_i) / Σ max(a_i, b_i)；
        - 对高分 spot 天然赋予更高权重（重点关注高 niche 区域），
          弥补 Spearman 对顶部区域与底部区域等权重的不足；
        - 值域 [0, 1]，越接近 1.0 说明评分分布高度相似。

    两指标联合解读：
      - spearman_rho > 0.95 且 weighted_jaccard > 0.90：结果极稳健；
      - spearman_rho > 0.90 且 weighted_jaccard > 0.80：结果稳健，参数不敏感；
      - spearman_rho < 0.85 或 weighted_jaccard < 0.70：参数变化有实质影响，
        建议重新审视参数选择。

    参数
    ----
    proportions             : 细胞类型比例矩阵
    coords                  : 空间坐标
    gene_score              : 免疫抑制基因模块评分
    treg_col                : Treg 比例列名
    myeloid_col             : Myeloid 比例列名
    fibroblast_col          : Fibroblast 比例列名
    n_neighbors_list        : k-NN 邻居数扫描列表
    radius_multiplier_list  : 半径倍增系数扫描列表
    niche_high_quantile     : 高 niche 区域分位数阈值（仅用于记录 n_niche_high，
                              不再作为相似度计算的截断点）
    ref_k                   : 参考参数（主分析 kNN k 值，作为基准评分向量）

    返回
    ----
    pd.DataFrame，每行为一个参数组合的统计结果，包含：
      param_mode / param_label / is_reference / n_niche_high / niche_pct /
      spearman_rho / weighted_jaccard / score_std
    """
    from scipy.stats import spearmanr as _spearmanr

    gs = gene_score.reindex(proportions.index).fillna(0.0)

    def _compute_score(
        nbrs: list[np.ndarray],
    ) -> pd.Series:
        """
        计算给定邻居列表下的综合 niche 连续评分（不截断）。
        返回 score Series，索引与 proportions 一致。
        """
        nc = _compute_neighborhood_composition(proportions, nbrs)
        s = (
            _zscore(nc[treg_col])
            + _zscore(nc[myeloid_col])
            + _zscore(nc[fibroblast_col])
            + _zscore(gs)
        )
        return s

    def _weighted_jaccard(a: np.ndarray, b: np.ndarray) -> float:
        """
        计算两个非负连续向量的 Weighted Jaccard 相似度。

        先做 min-max 归一化将值域映射到 [0, 1]，再计算
          WJ = Σ min(a_i, b_i) / Σ max(a_i, b_i)
        对高分元素天然赋予更高权重，适合「关注高分区域一致性」的场景。
        若分母为 0（所有元素均为 0），返回 1.0（视为完全一致）。
        """
        def _minmax(v: np.ndarray) -> np.ndarray:
            lo, hi = v.min(), v.max()
            return (v - lo) / (hi - lo) if hi > lo else np.zeros_like(v)

        a_n = _minmax(a)
        b_n = _minmax(b)
        denom = np.maximum(a_n, b_n).sum()
        if denom == 0:
            return 1.0
        return float(np.minimum(a_n, b_n).sum() / denom)

    # ── 计算参考评分向量（主分析 k=ref_k）────────────────────────────────────
    ref_nbrs = _build_knn_neighbors(coords, k=ref_k)
    try:
        ref_score = _compute_score(ref_nbrs)
        ref_arr = ref_score.to_numpy(dtype=float)
    except Exception as exc:
        logging.warning("Could not compute reference score for k=%d: %s", ref_k, exc)
        ref_arr = None

    rows = []
    all_params: list[tuple[str, int | None, float | None]] = (
        [("knn", k, None) for k in n_neighbors_list]
        + [("radius", None, rm) for rm in radius_multiplier_list]
    )

    for mode, k, rm in all_params:
        try:
            if mode == "knn":
                nbrs = _build_knn_neighbors(coords, k=k)
                label = f"kNN k={k}"
                is_ref = (k == ref_k)
            else:
                nbrs, _ = _spatial_neighbors(coords, radius_multiplier=rm)
                label = f"radius×{rm}"
                is_ref = False

            score = _compute_score(nbrs)
            curr_arr = score.to_numpy(dtype=float)

            # 连续评分相似度（无硬截断）
            if ref_arr is not None and len(ref_arr) == len(curr_arr):
                # 指标一：Spearman 秩相关（全局排序稳定性）
                rho, _ = _spearmanr(ref_arr, curr_arr)
                spearman_rho = round(float(rho), 4)

                # 指标二：Weighted Jaccard（高分区域加权重叠度）
                wj = _weighted_jaccard(ref_arr, curr_arr)
                weighted_jac = round(wj, 4)
            else:
                spearman_rho = 1.0 if is_ref else float("nan")
                weighted_jac = 1.0 if is_ref else float("nan")

            # 保留 n_niche_high 用于信息记录（不参与相似度计算）
            thr = float(score.quantile(niche_high_quantile))
            n_high = int((score >= thr).sum())
            pct = n_high / len(score) * 100

            rows.append({
                "param_mode":      mode,
                "param_label":     label,
                "is_reference":    is_ref,
                "n_niche_high":    n_high,
                "niche_pct":       pct,
                "spearman_rho":    spearman_rho,
                "weighted_jaccard": weighted_jac,
                "score_std":       float(score.std()),
            })
            logging.info(
                "Sensitivity [%s]: n_high=%d (%.1f%%), Spearman_rho=%.4f, "
                "Weighted_Jaccard=%.4f%s",
                label, n_high, pct, spearman_rho, weighted_jac,
                " [REF]" if is_ref else "",
            )
        except Exception as exc:
            logging.warning("Sensitivity analysis skipped for %s: %s", label, exc)

    return pd.DataFrame(rows)


# ============================================================
# 可视化函数
# ============================================================

def _spatial_scatter(
    df: pd.DataFrame,
    value: str,
    path: Path,
    title: str,
    cmap: str = "viridis",
    categorical: bool = False,
) -> None:
    """
    绘制空间散点图（每个点代表一个 spot，颜色编码特征值或类别标签）。

    连续值模式（categorical=False）：颜色映射 + 颜色条；
    分类变量模式（categorical=True）：固定颜色 + 图例。
    """
    fig, ax = plt.subplots(figsize=(7, 6))
    if categorical:
        palette = {
            "tumor_core":              "#d62728",
            "tumor_edge":              "#ff7f0e",
            "stroma_immune":           "#1f77b4",
            "immunosuppressive_niche": "#9467bd",
            "other":                   "#bdbdbd",
        }
        colors = df[value].map(palette).fillna("#bdbdbd")
        ax.scatter(df["spatial_x"], df["spatial_y"], c=colors, s=18, linewidths=0)
        present = df[value].unique()
        handles = [
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=palette.get(lab, "#bdbdbd"),
                   label=lab, markersize=8)
            for lab in present if lab in palette
        ]
        ax.legend(handles=handles, frameon=False, loc="best", fontsize=8)
    else:
        sc_obj = ax.scatter(
            df["spatial_x"], df["spatial_y"],
            c=df[value], s=18, cmap=cmap, linewidths=0,
        )
        fig.colorbar(sc_obj, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title)
    ax.set_xlabel("spatial x")
    ax.set_ylabel("spatial y")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_neighborhood_clusters(df: pd.DataFrame, path: Path) -> None:
    """
    绘制邻域组成聚类（Cellular Neighborhoods）的空间分布图。

    用 tab20 颜色板自动为每个 cluster 分配颜色，图例标注 cluster 编号。
    """
    clusters = sorted(df["neighborhood_cluster"].unique(), key=lambda x: str(x))
    cmap_nc = plt.get_cmap("tab20", len(clusters))
    palette = {c: mcolors.to_hex(cmap_nc(i)) for i, c in enumerate(clusters)}
    colors = df["neighborhood_cluster"].map(palette).fillna("#bdbdbd")

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df["spatial_x"], df["spatial_y"], c=colors, s=18, linewidths=0)
    handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=palette[c], label=f"NC-{c}", markersize=7)
        for c in clusters
    ]
    ax.legend(handles=handles, frameon=False, loc="best",
              fontsize=7, ncol=2, title="Neighborhood\nCluster")
    ax.set_title("Neighborhood Composition Clusters\n(Leiden, Schürch et al. 2020)")
    ax.set_xlabel("spatial x")
    ax.set_ylabel("spatial y")
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_distance(df: pd.DataFrame, path: Path) -> None:
    """
    绘制"到高肝细胞区域的距离"与"免疫抑制 niche 评分"的折线图。

    将距离按分位数分成5个等频区间，
    计算各区间内 niche 评分的均值 ± 标准误差，
    用于验证免疫抑制信号的距离依赖性。
    """
    cols = ["distance_to_hep_high", "immunosuppressive_niche_score"]
    use = df.dropna(subset=cols).copy()
    if use.empty:
        return
    use["distance_bin"] = pd.qcut(use["distance_to_hep_high"], q=5, duplicates="drop")
    summary = (
        use.groupby("distance_bin", observed=True)["immunosuppressive_niche_score"]
        .agg(["mean", "sem", "count"])
        .reset_index()
    )
    summary["bin"] = np.arange(1, len(summary) + 1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(summary["bin"], summary["mean"], yerr=summary["sem"],
                marker="o", capsize=3, color="#5b9bd5")
    ax.set_xlabel("Distance bin from Hepatocyte-high area")
    ax.set_ylabel("Mean immunosuppressive niche score")
    ax.set_xticks(summary["bin"])
    ax.set_xticklabels([str(v) for v in summary["distance_bin"]], rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_region_box(df: pd.DataFrame, path: Path) -> None:
    """
    绘制不同空间区域评分的箱线图（Treg_like_score 和 immunosuppressive_niche_score）。

    区域顺序：tumor_core → tumor_edge → stroma_immune → other。
    """
    order = ["tumor_core", "tumor_edge", "stroma_immune", "other"]
    use = df[df["spatial_region"].isin(order)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    sns.boxplot(data=use, x="spatial_region", y="Treg_like_score",
                order=order, ax=axes[0])
    sns.boxplot(data=use, x="spatial_region", y="immunosuppressive_niche_score",
                order=order, ax=axes[1])
    for ax in axes:
        ax.tick_params(axis="x", rotation=30)
        ax.set_xlabel("")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_correlation(df: pd.DataFrame, columns: list[str], path: Path) -> None:
    """
    绘制指定列之间的 Pearson 相关系数热图（红-白-蓝 vlag 配色）。

    用于探索细胞类型比例与 niche 评分之间的共定位关系。
    """
    corr = df[columns].corr()
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(corr, cmap="vlag", center=0, annot=True,
                fmt=".2f", square=True, ax=ax)
    ax.set_title("Cell abundance and niche score correlation")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_hep_treg(df: pd.DataFrame, path: Path) -> None:
    """
    绘制 Hepatocyte 比例 vs Treg 比例散点图，颜色编码 niche 评分（magma 配色）。

    虚线标注 0.75 分位数，划分四象限，直观展示肝细胞-Treg 共定位关系。
    """
    fig, ax = plt.subplots(figsize=(5.5, 5))
    sc_obj = ax.scatter(
        df["Hepatocyte"], df["Treg"],
        c=df["immunosuppressive_niche_score"],
        s=20, cmap="magma", linewidths=0,
    )
    ax.axvline(df["Hepatocyte"].quantile(0.75), color="black", linestyle="--", linewidth=1)
    ax.axhline(df["Treg"].quantile(0.75),       color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Hepatocyte proportion")
    ax.set_ylabel("Treg proportion")
    fig.colorbar(sc_obj, ax=ax, label="niche score")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _resolve_gene(gene: str, available: set[str]) -> tuple[str, bool]:
    """
    尝试在数据集中找到基因本身或其同家族替代基因。

    先检查主基因是否存在，若不存在则按 GENE_FALLBACKS 字典依次尝试替代基因，
    返回第一个找到的基因名及是否使用了替代基因的标志。

    参数
    ----
    gene      : 主基因名称
    available : 数据集中所有基因名称的集合

    返回
    ----
    (resolved_gene, is_fallback)：
      resolved_gene : 找到的基因名（主基因或替代基因）
      is_fallback   : True 表示使用了替代基因
    若主基因和所有替代基因均不存在，返回 (gene, False) 但后续调用方应检查基因是否实际存在。
    """
    if gene in available:
        return gene, False
    for fb in GENE_FALLBACKS.get(gene, []):
        if fb in available:
            logging.info("L-R fallback: %s → %s", gene, fb)
            return fb, True
    return gene, False


def _plot_lr_communication(
    adata: ad.AnnData,
    df: pd.DataFrame,
    lr_pairs: list[tuple[str, str, str]],
    path: Path,
) -> None:
    """
    绘制配体-受体（L-R）通讯分析热图，比较 niche_high vs niche_low 的信号强度。

    【改进说明】
    原版本使用 product = ligand × receptor 直接计算通讯强度，存在缺陷：
      - 任一基因在 spot 级别表达接近 0（细胞稀释效应），乘积即为 0；
      - CCL22、IL10 等 Treg 特异基因在 Visium spot 中几乎不可检测，
        导致整行全为 0，热图失去对比性。

    改进后采用 rank-normalized product（秩归一化乘积）策略：
      1. 对每个基因的表达向量做秩归一化（rankdata / n_spots），
         将值域映射到 [0, 1]；
      2. 计算配体秩 × 受体秩的均值作为通讯强度代理，
         避免零值乘积问题；
      3. 对缺失基因自动尝试 GENE_FALLBACKS 中的同家族替代基因，
         并在标签中注明替代（用 * 标注）。

    综述依据：CXCL12-CXCR4、CCL22-CCR4 是 Treg 招募的关键通讯轴
    （Efremova et al., 2020; Cang & Nie, 2023）；
    SPP1-CD44 是 TAM 介导的免疫抑制重要通路（Zhu et al., 2022）。
    """
    from scipy.stats import rankdata as _rankdata

    expr = _expression_frame(adata)
    available_genes = set(expr.columns)
    high_mask = df["niche_high"].astype(bool).reindex(
        pd.Index(range(len(df)))
    ).fillna(False)

    # 对齐索引：df 的行索引可能是 spot_id，需按位置提取 mask
    if df.index.equals(pd.RangeIndex(len(df))):
        high_arr = df["niche_high"].astype(bool).to_numpy()
    else:
        high_arr = df["niche_high"].astype(bool).to_numpy()
    low_arr = ~high_arr

    records = []
    for ligand, receptor, label in lr_pairs:
        # 尝试 fallback 基因解析
        actual_lig, lig_fb = _resolve_gene(ligand, available_genes)
        actual_rec, rec_fb = _resolve_gene(receptor, available_genes)

        missing = [g for g, actual in [(ligand, actual_lig), (receptor, actual_rec)]
                   if actual not in available_genes]
        if missing:
            logging.warning("L-R [%s]: genes %s not found (no fallback); skipping.", label, missing)
            continue

        # Rank-normalized product：对稀疏高表达基因更鲁棒
        lig_vals = expr[actual_lig].to_numpy(dtype=float)
        rec_vals = expr[actual_rec].to_numpy(dtype=float)
        n = len(lig_vals)
        lig_ranked = _rankdata(lig_vals) / n
        rec_ranked = _rankdata(rec_vals) / n
        product = lig_ranked * rec_ranked   # 秩归一化乘积，值域 [0, 1]

        mean_high = float(product[high_arr].mean()) if high_arr.any() else 0.0
        mean_low  = float(product[low_arr].mean())  if low_arr.any()  else 0.0
        log2fc    = float(np.log2((mean_high + 1e-6) / (mean_low + 1e-6)))

        # 在标签中注明是否使用了替代基因
        display_label = label
        if lig_fb or rec_fb:
            display_label = f"{label}*"
            logging.info("L-R [%s] → %s–%s (fallback used)", label, actual_lig, actual_rec)

        records.append({
            "lr_pair":   display_label,
            "mean_high": mean_high,
            "mean_low":  mean_low,
            "log2fc":    log2fc,
        })

    if not records:
        logging.warning("No valid L-R pairs; skipping L-R communication heatmap.")
        return

    res_df = pd.DataFrame(records).set_index("lr_pair")
    hm_data = res_df[["mean_high", "mean_low"]].rename(
        columns={"mean_high": "Niche-High", "mean_low": "Niche-Low"}
    )
    # 行最大值归一化（0–1），便于颜色对比
    row_max = hm_data.max(axis=1).replace(0, 1.0)
    hm_norm = hm_data.div(row_max, axis=0)

    n_rows = len(records)
    fig, axes = plt.subplots(
        1, 2,
        figsize=(10, max(4, n_rows * 0.65 + 1.5)),
        gridspec_kw={"width_ratios": [3, 1]},
    )
    # 左图：信号强度热图
    sns.heatmap(
        hm_norm, cmap="YlOrRd", vmin=0, vmax=1,
        annot=hm_data.round(4), fmt="g",
        linewidths=0.5, ax=axes[0],
    )
    axes[0].set_title(
        "L-R Communication Signal Intensity\n(Niche-High vs Niche-Low)",
        fontsize=10,
    )
    axes[0].set_xlabel("Niche Group")
    axes[0].set_ylabel("Ligand–Receptor Pair")

    # 右图：log2FC 条形图
    lfc_vals = res_df["log2fc"]
    bar_colors = ["#d73027" if v >= 0 else "#4575b4" for v in lfc_vals]
    axes[1].barh(lfc_vals.index, lfc_vals.values, color=bar_colors, edgecolor="none")
    axes[1].axvline(0, color="black", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("log2FC\n(High / Low)")
    axes[1].set_title("Enrichment\nin Niche-High", fontsize=9)
    axes[1].set_yticks([])

    fig.suptitle(
        "Treg–TAM–CAF Ligand-Receptor Communication",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_sensitivity(sensitivity_df: pd.DataFrame, path: Path) -> None:
    """
    绘制敏感性分析结果图（v2）：双指标展示不同参数设置下 niche 评分的连续稳定性。

    【改进说明 v2 — 2026-06-26】
    v1 版本改用了 Jaccard 相似度（硬截断集合指标），但仍存在三个缺陷：
      1. 只判断是否过阈值，丢失评分量级信息（Magnitude Loss）；
      2. 固定边缘概率使 Jaccard 变动范围极窄；
      3. 边界微小扰动引入硬截断噪声。

    v2 版本改用两个连续型指标，彻底消除硬截断：

      左图组（kNN）/ 右图组（radius）各含两个子图：
        上子图：Spearman 秩相关系数（全局排序稳定性，权重均匀）
          - 稳健区间 > 0.90，优秀区间 > 0.95
        下子图：Weighted Jaccard（高分区域加权重叠度，高分 spot 权重更高）
          - 稳健区间 > 0.80，优秀区间 > 0.90

    红色柱 = 参考参数（k=15 kNN）；蓝色柱 = 其他参数。
    橙色虚线 = 各指标的稳健性阈值；红色虚线 = 参考值（=1.0）。
    """
    if sensitivity_df.empty:
        return

    has_spearman = "spearman_rho" in sensitivity_df.columns
    has_wj       = "weighted_jaccard" in sensitivity_df.columns

    # 若列名为旧版（jaccard_vs_ref），退化为单指标兼容模式
    if not has_spearman and not has_wj:
        legacy_col = "jaccard_vs_ref" if "jaccard_vs_ref" in sensitivity_df.columns else "niche_pct"
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for ax, mode in zip(axes, ["knn", "radius"]):
            sub = sensitivity_df[sensitivity_df["param_mode"] == mode].copy()
            if sub.empty:
                ax.set_visible(False)
                continue
            colors = ["#d62728" if bool(r.get("is_reference", False)) else "#5b9bd5"
                      for _, r in sub.iterrows()]
            ax.bar(sub["param_label"], sub[legacy_col], color=colors, edgecolor="none")
            ax.set_xlabel("Parameter Setting")
            ax.set_ylabel(legacy_col)
            ax.tick_params(axis="x", rotation=20)
        fig.suptitle("Sensitivity Analysis (legacy mode)", fontsize=10)
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return

    # ── 主绘图逻辑：2×2 布局（kNN 左列 / radius 右列；Spearman 上行 / WJ 下行）──
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    mode_list  = ["knn", "radius"]
    mode_title = ["kNN Neighbor Count (k)", "Radius Multiplier"]

    metric_list  = ["spearman_rho", "weighted_jaccard"]
    metric_label = [
        "Spearman ρ (rank correlation vs ref k=15)",
        "Weighted Jaccard (score overlap vs ref k=15)",
    ]
    metric_thresh = [0.90, 0.80]   # 各指标稳健性阈值
    metric_colors_ok = ["#1a9641", "#1a9641"]  # 阈值线颜色

    for col_idx, (mode, mtitle) in enumerate(zip(mode_list, mode_title)):
        sub = sensitivity_df[sensitivity_df["param_mode"] == mode].copy()

        for row_idx, (metric, ylabel, thresh) in enumerate(
            zip(metric_list, metric_label, metric_thresh)
        ):
            ax = axes[row_idx, col_idx]

            if sub.empty or metric not in sub.columns:
                ax.set_visible(False)
                continue

            # 参考参数用红色高亮
            colors = [
                "#d62728" if bool(r.get("is_reference", False)) else "#5b9bd5"
                for _, r in sub.iterrows()
            ]

            bars = ax.bar(
                sub["param_label"], sub[metric],
                color=colors, edgecolor="none", width=0.5,
            )

            # 在柱顶标注数值
            for bar, (_, row) in zip(bars, sub.iterrows()):
                val = row[metric]
                if not np.isnan(val):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.005,
                        f"{val:.3f}",
                        ha="center", va="bottom", fontsize=8,
                    )

            # 稳健性阈值线
            ax.axhline(thresh, color="orange", linestyle="--",
                       linewidth=1.2, label=f"threshold={thresh}")
            # 参考值线（=1.0）
            ax.axhline(1.0, color="#d62728", linestyle="--",
                       linewidth=0.8, alpha=0.5, label="ref (k=15) = 1.0")

            ax.set_ylim(max(0, sub[metric].min() - 0.1), 1.05)
            ax.set_xlabel("Parameter Setting", fontsize=8)
            ax.set_ylabel(ylabel, fontsize=8)
            ax.set_title(f"{mtitle}\n({metric})", fontsize=9)
            ax.tick_params(axis="x", rotation=20, labelsize=8)
            ax.legend(fontsize=7, frameon=False)

    fig.suptitle(
        "Sensitivity Analysis: Niche Score Stability Under Parameter Variation\n"
        "(Continuous metrics — no hard thresholding; red bar = reference k=15)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# 特征基因提取（升级版）
# ============================================================

def _gini_index(values: np.ndarray) -> float:
    """
    计算给定数值数组的 Gini Index（基尼系数），衡量表达不均匀性。

    Gini Index 取值范围 [0, 1]：
      - 0：完全均匀（每个 spot 表达量相同）
      - 1：完全集中（只有一个 spot 有表达量）

    对于空间上局灶性高表达的稀有基因（如 FOXP3），Gini 值偏高，
    能弥补 log2FC 因均值稀释而失效的不足。

    参数
    ----
    values : 数值数组（如某基因在一组 spot 中的表达量）

    返回
    ----
    float，Gini 系数。
    """
    v = np.sort(np.abs(values).ravel())
    n = len(v)
    if n == 0 or v.sum() == 0:
        return 0.0
    k = np.arange(1, n + 1, dtype=float)
    return float(1.0 - (2.0 * (k * v).sum()) / (n * v.sum()) + 1.0 / n)


def _rank_niche_genes(
    adata: ad.AnnData,
    labels: pd.Series,
    top_n: int,
) -> pd.DataFrame:
    """
    对 niche 高分 spot 进行多维度差异基因分析，筛选特征性高表达基因。

    【设计背景：为什么 FOXP3 等目标基因不会被 Wilcoxon+FC 选中？】
    Visium 每个 spot 覆盖 5-50 个细胞，FOXP3、CTLA4 等 Treg 标志基因在绝大多数
    spot 中表达值为 0（即使该 spot 内确实含有少量 Treg），导致：
      1. 均值（mean_high / mean_low）均接近 0，log2FC 被稀释至接近 0；
      2. 两组均为大量 0 值时，Wilcoxon 检验产生大量秩并列，p 值虽小，但排名
         被高表达的管家/结构基因挤占，目标基因无缘 Top-N；
      3. 全局 Wilcoxon 将所有 spot 作为对象，无法识别"在少数 spot 中偶发高表达"
         的稀有免疫细胞标志物。

    【解决策略：引入检出率（Fraction of spots）作为主要筛选维度】
    检出率（frac）= 表达量 > 0 的 spot 占各组总数的比例。
    对于 FOXP3 这类基因：可能在 niche_high 组中有 30% 的 spot 能检出，
    而在 niche_low 中仅 3% 能检出。这个 delta_frac=0.27 的差异远比均值
    差异更具生物学意义，且不受稀释效应影响。

    【三层筛选策略（并行运行）】
      Layer 1（Wilcoxon + FDR + 综合排序）：
        - 对每个基因做 Mann-Whitney U 检验（等价于 Wilcoxon 秩和检验）；
        - 用 Benjamini-Hochberg 方法对 p 值进行 FDR 校正；
        - 新增 frac_high / frac_low / delta_frac（检出率及其差值）；
        - 综合排序分 = 0.5 × norm(log2FC) + 0.5 × norm(delta_frac)，
          使低均值但高富集度的基因（如 FOXP3）也能进入 Top-N；
        - 主筛选条件：FDR < 0.05 且（log2FC > 0.5 OR delta_frac > 0.10）；
        - 补充条件：若严格条件不足 Top-N，放宽至 FDR < 0.2 且
                     （log2FC > 0.3 OR delta_frac > 0.05）；
        - 产出：immunosuppressive_niche_signature_genes_ranked.csv（含新列）
                + niche_signature_volcano.png（FC vs -log10 FDR 火山图）
                + niche_fraction_scatter.png（delta_frac vs log2FC 散点图）

      Layer 2（先验功能基因集 AUC + Fraction 检验）：
        - 对 PRIOR_GENE_SETS 中的三组目标基因（Treg/TAM/CAF）分别做检验；
        - 报告每个基因的 AUC、p 值（FDR 校正）、均值表达量和检出率差值；
        - AUC > 0.6 且 FDR < 0.05（或 delta_frac > 0.10）视为有意义的
          niche 标志基因；
        - 满足条件的先验基因强制合并进 Layer 1 结果（补充入 Top-N）；
        - 产出：prior_gene_set_auc.csv

      Layer 3（Gini Index 特异性评分）：
        - 对 niche_high spot 中每个基因计算 Gini 系数；
        - Gini 高说明该基因在少数 spot 中高度富集（局灶性表达）；
        - 与 log2FC > 0 联合筛选，识别高度特异的稀有免疫基因；
        - 产出：gini_score_genes.csv

    函数本身返回 Layer 1 + Layer 2 合并后的结果，其余层结果通过额外属性附加
    到返回 DataFrame 上（df.attrs 字典），供调用方保存。

    参数
    ----
    adata  : AnnData 数据对象（提供基因表达矩阵）
    labels : 布尔 Series，True = niche_high
    top_n  : 返回前 top_n 个基因

    返回
    ----
    pd.DataFrame，列：
      gene / mean_high / mean_low / log2_fc / frac_high / frac_low /
      delta_frac / composite_score / pvalue / fdr
    按 composite_score 降序（FDR < 0.05 或 delta_frac > 0.10 的基因优先）。
    DataFrame.attrs 额外包含：
      "prior_auc_df" : 先验基因集 AUC 检验结果 DataFrame
      "gini_df"      : Gini Index 评分 DataFrame（niche_high 子集）
    """
    expr = _expression_frame(adata)
    high = labels.astype(bool).reindex(expr.index).fillna(False)

    if int(high.sum()) < 3 or int((~high).sum()) < 3:
        logging.warning("Too few spots in niche_high/niche_low; skipping signature ranking.")
        empty = pd.DataFrame(columns=[
            "gene", "mean_high", "mean_low", "log2_fc",
            "frac_high", "frac_low", "delta_frac",
            "composite_score", "pvalue", "fdr",
        ])
        empty.attrs["prior_auc_df"] = pd.DataFrame()
        empty.attrs["gini_df"] = pd.DataFrame()
        return empty

    high_expr = expr.loc[high]
    low_expr  = expr.loc[~high]
    mean_high = high_expr.mean(axis=0)
    mean_low  = low_expr.mean(axis=0)

    # ── 检出率（Fraction of spots with expression > 0） ──────────────────────
    # 这是解决稀释效应的核心指标：不依赖均值大小，只看有多少比例的 spot 能"检出"该基因
    frac_high = (high_expr > 0).mean(axis=0)  # niche_high 组中检出该基因的 spot 比例
    frac_low  = (low_expr  > 0).mean(axis=0)  # niche_low  组中检出该基因的 spot 比例
    delta_frac = frac_high - frac_low          # 检出率差值（正值表示 niche_high 中更多 spot 能检出）

    # ── Layer 1：log2FC + Wilcoxon + FDR ──────────────────────────────────────
    log2fc = np.log2((mean_high + 1.0) / (mean_low + 1.0))
    pvalues = np.full(len(expr.columns), 1.0)

    # 候选基因：log2FC > 0.1 OR delta_frac > 0.05（兼顾高表达基因和高富集基因）
    candidate_mask = (log2fc.to_numpy() > 0.1) | (delta_frac.to_numpy() > 0.05)
    candidate_genes = expr.columns[candidate_mask].tolist()
    logging.info(
        "Running Wilcoxon test on %d candidate genes (log2FC>0.1 OR delta_frac>0.05)...",
        len(candidate_genes),
    )
    for i, gene in enumerate(candidate_genes):
        gene_idx = list(expr.columns).index(gene)
        h_vals = high_expr[gene].to_numpy()
        l_vals = low_expr[gene].to_numpy()
        if h_vals.std() == 0 and l_vals.std() == 0:
            continue
        try:
            _, p = scipy_stats.mannwhitneyu(h_vals, l_vals, alternative="greater")
            pvalues[gene_idx] = p
        except Exception:
            pass

    _, fdr, _, _ = multipletests(pvalues, method="fdr_bh")

    # ── 综合排序分（composite_score）─────────────────────────────────────────
    # 同时考虑 log2FC（表达量倍数差异）和 delta_frac（检出率差异），
    # 使 FOXP3 等低表达但高富集的基因也能获得合理排名。
    # Min-Max 归一化到 [0, 1] 再各赋 0.5 权重。
    log2fc_arr     = log2fc.to_numpy()
    delta_frac_arr = delta_frac.to_numpy()

    def _minmax_norm(v: np.ndarray) -> np.ndarray:
        lo, hi = v.min(), v.max()
        return (v - lo) / (hi - lo + 1e-12)

    composite = 0.5 * _minmax_norm(log2fc_arr) + 0.5 * _minmax_norm(delta_frac_arr)

    ranked = pd.DataFrame({
        "gene":            expr.columns.tolist(),
        "mean_high":       mean_high.to_numpy(),
        "mean_low":        mean_low.to_numpy(),
        "log2_fc":         log2fc_arr,
        "frac_high":       frac_high.to_numpy(),
        "frac_low":        frac_low.to_numpy(),
        "delta_frac":      delta_frac_arr,
        "composite_score": composite,
        "pvalue":          pvalues,
        "fdr":             fdr,
    })

    # 主签名基因：FDR < 0.05 且（log2FC > 0.5 OR delta_frac > 0.10）
    sig_strict = ranked[
        (ranked["fdr"] < 0.05) &
        ((ranked["log2_fc"] > 0.5) | (ranked["delta_frac"] > 0.10))
    ]
    sig_strict = sig_strict.sort_values("composite_score", ascending=False)

    # 若严格条件下不足 top_n，放宽至 FDR < 0.2 且（log2FC > 0.3 OR delta_frac > 0.05）
    if len(sig_strict) < top_n:
        loose = ranked[
            (ranked["fdr"] < 0.2) &
            ((ranked["log2_fc"] > 0.3) | (ranked["delta_frac"] > 0.05))
        ]
        loose = loose.sort_values("composite_score", ascending=False)
        sig_strict = pd.concat([sig_strict, loose]).drop_duplicates("gene")

    result = sig_strict.head(top_n).reset_index(drop=True)

    # ── Layer 2：先验功能基因集 AUC 检验 ──────────────────────────────────────
    prior_rows = []
    for group_name, gene_list in PRIOR_GENE_SETS.items():
        for gene in gene_list:
            # fallback：主基因不存在时尝试替代基因
            actual_gene = gene
            if gene not in expr.columns:
                for fb in GENE_FALLBACKS.get(gene, []):
                    if fb in expr.columns:
                        actual_gene = fb
                        logging.info("Prior gene set: %s → fallback to %s", gene, fb)
                        break
            if actual_gene not in expr.columns:
                logging.debug("Prior gene set [%s]: %s not found in data; skipped.", group_name, gene)
                continue

            h_vals = high_expr[actual_gene].to_numpy()
            l_vals = low_expr[actual_gene].to_numpy()
            # AUC（用 Mann-Whitney U 统计量计算）
            try:
                stat, p = scipy_stats.mannwhitneyu(h_vals, l_vals, alternative="greater")
                auc = stat / (len(h_vals) * len(l_vals))
            except Exception:
                stat, p, auc = 0.0, 1.0, 0.5

            # 检出率：该先验基因在两组 spot 中的表达检出率及差值
            frac_h = float((h_vals > 0).mean())
            frac_l = float((l_vals > 0).mean())
            d_frac = frac_h - frac_l

            prior_rows.append({
                "gene_set":     group_name,
                "gene":         gene,
                "actual_gene":  actual_gene,
                "mean_high":    float(h_vals.mean()),
                "mean_low":     float(l_vals.mean()),
                "log2_fc":      float(np.log2((h_vals.mean() + 1.0) / (l_vals.mean() + 1.0))),
                "frac_high":    frac_h,
                "frac_low":     frac_l,
                "delta_frac":   d_frac,
                "pvalue":       float(p),
                "auc":          float(auc),
            })

    prior_df = pd.DataFrame(prior_rows)
    if not prior_df.empty and "pvalue" in prior_df.columns:
        _, prior_fdr, _, _ = multipletests(prior_df["pvalue"].to_numpy(), method="fdr_bh")
        prior_df["fdr"] = prior_fdr

    n_sig_prior = 0
    if not prior_df.empty:
        n_sig_prior = int(
            ((prior_df.get("auc", pd.Series()) > 0.6) &
             ((prior_df.get("fdr", pd.Series(1.0)) < 0.05) |
              (prior_df.get("delta_frac", pd.Series(0.0)) > 0.10))).sum()
        )
    logging.info(
        "Prior gene set results: %d genes tested, %d with AUC>0.6 & (FDR<0.05 or delta_frac>0.10)",
        len(prior_df),
        n_sig_prior,
    )

    # ── Layer 2 → 强制合并入 result：将满足条件的先验基因补充入签名列表 ─────
    # 解决问题：FOXP3 等目标基因可能因 log2FC 过低而无法通过 Layer 1 筛选，
    # 但其 delta_frac（检出率差值）或 AUC 体现了真实的生物学富集。
    # 在此强制将 AUC > 0.6 且 delta_frac > 0.05 的先验基因纳入最终列表。
    if not prior_df.empty:
        prior_sig_genes = prior_df[
            (prior_df.get("auc", pd.Series(0.0)) > 0.6) &
            (prior_df.get("delta_frac", pd.Series(0.0)) > 0.05)
        ]["actual_gene"].dropna().unique().tolist()

        if prior_sig_genes:
            logging.info(
                "Layer 2 forced-merge: %d prior genes will be added to signature "
                "(AUC>0.6 & delta_frac>0.05): %s",
                len(prior_sig_genes),
                ", ".join(prior_sig_genes[:10]),
            )
            # 从 ranked 中取出这些基因的全部统计量
            forced_rows = ranked[ranked["gene"].isin(prior_sig_genes)].copy()
            # 合并：已在 result 中的不重复添加
            already = set(result["gene"].tolist())
            new_rows = forced_rows[~forced_rows["gene"].isin(already)]
            if not new_rows.empty:
                result = pd.concat([result, new_rows], ignore_index=True)
                result = result.sort_values("composite_score", ascending=False).reset_index(drop=True)

    # ── Layer 3：Gini Index 特异性评分 ─────────────────────────────────────────
    gini_rows = []
    high_arr = high_expr.to_numpy()
    for gene in expr.columns:
        gene_idx = list(expr.columns).index(gene)
        g_vals = high_arr[:, gene_idx]
        l2fc = float(log2fc.iloc[gene_idx])
        if l2fc <= 0:
            continue  # 仅对在 niche_high 中更高表达的基因计算 Gini
        gini = _gini_index(g_vals)
        if gini > 0.3:  # 保留 Gini 较高（局灶性表达）的基因；0.3 覆盖稀有免疫基因
            gini_rows.append({
                "gene":    gene,
                "gini":    round(gini, 4),
                "log2_fc": round(l2fc, 4),
                "mean_high": float(mean_high.iloc[gene_idx]),
            })
    gini_df = pd.DataFrame(gini_rows).sort_values("gini", ascending=False) if gini_rows \
        else pd.DataFrame(columns=["gene", "gini", "log2_fc", "mean_high"])
    logging.info(
        "Gini index: %d genes with Gini>0.5 and log2FC>0", len(gini_df)
    )

    result.attrs["prior_auc_df"] = prior_df
    result.attrs["gini_df"] = gini_df
    return result


def _plot_volcano(
    ranked: pd.DataFrame,
    path: Path,
    fc_threshold: float = 0.5,
    fdr_threshold: float = 0.05,
    highlight_genes: list[str] | None = None,
) -> None:
    """
    绘制差异基因火山图（Volcano Plot）。

    横轴：log2FC（niche_high vs niche_low）
    纵轴：-log10(FDR)
    颜色编码：
      红色 = 上调显著（log2FC > fc_threshold 且 FDR < fdr_threshold）
      蓝色 = 下调显著
      灰色 = 不显著

    参数
    ----
    ranked          : _rank_niche_genes 返回的 DataFrame（含 log2_fc、fdr 列）
    path            : 输出图片路径
    fc_threshold    : log2FC 显著性阈值（默认 0.5）
    fdr_threshold   : FDR 显著性阈值（默认 0.05）
    highlight_genes : 需要特别标注名称的基因列表（如免疫抑制目标基因）
    """
    if "fdr" not in ranked.columns or ranked.empty:
        logging.warning("No FDR column found; skipping volcano plot.")
        return

    if highlight_genes is None:
        # 默认高亮先验免疫抑制目标基因
        highlight_genes = [
            g for gs in PRIOR_GENE_SETS.values() for g in gs
        ] + ["FOXP3", "TGFB1", "FAP", "ACTA2", "CCL22", "CXCL12", "SPP1"]

    df = ranked.copy()
    df["neg_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))

    # 分类着色
    conditions = [
        (df["log2_fc"] > fc_threshold) & (df["fdr"] < fdr_threshold),
        (df["log2_fc"] < -fc_threshold) & (df["fdr"] < fdr_threshold),
    ]
    colors_map = {True: "#d62728", False: {True: "#4575b4", False: "#bdbdbd"}}

    c_list = []
    for _, row in df.iterrows():
        if row["log2_fc"] > fc_threshold and row["fdr"] < fdr_threshold:
            c_list.append("#d62728")
        elif row["log2_fc"] < -fc_threshold and row["fdr"] < fdr_threshold:
            c_list.append("#4575b4")
        else:
            c_list.append("#bdbdbd")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df["log2_fc"], df["neg_log10_fdr"],
               c=c_list, s=12, alpha=0.7, linewidths=0)
    ax.axvline(fc_threshold,  color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(fdr_threshold), color="black", linestyle="--", linewidth=0.8)

    # 标注目标基因
    labeled = set()
    for _, row in df.iterrows():
        gene = row.get("gene", "")
        if gene in highlight_genes and gene not in labeled:
            ax.annotate(
                gene, (row["log2_fc"], row["neg_log10_fdr"]),
                textcoords="offset points", xytext=(4, 2),
                fontsize=7, color="#333333",
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.5),
            )
            labeled.add(gene)

    n_up = int(((df["log2_fc"] > fc_threshold) & (df["fdr"] < fdr_threshold)).sum())
    n_dn = int(((df["log2_fc"] < -fc_threshold) & (df["fdr"] < fdr_threshold)).sum())
    ax.set_xlabel(f"log₂FC (niche_high / niche_low)")
    ax.set_ylabel(f"-log₁₀(FDR)")
    ax.set_title(
        f"Niche Signature Genes: Volcano Plot\n"
        f"Up={n_up}  Down={n_dn}  (|FC|>{fc_threshold}, FDR<{fdr_threshold})",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Volcano plot saved: %s", path)


# ============================================================
# 新增：检出率差值散点图（Fraction Scatter）
# ============================================================

def _plot_fraction_scatter(
    ranked: pd.DataFrame,
    path: Path,
    delta_frac_threshold: float = 0.10,
    fc_threshold: float = 0.5,
    highlight_genes: list[str] | None = None,
) -> None:
    """
    绘制"检出率差值 vs log2FC"散点图（Fraction Scatter Plot）。

    【设计目的】
    火山图以均值 log2FC 为横轴，在 Visium spot-level bulk 数据中，FOXP3 等
    稀有细胞标志基因的均值被大量 0 值压低，即使这些基因在 niche_high 中有更高
    的检出率（更多 spot 能检测到表达），也无法在火山图上突显。本图以
    delta_frac（检出率差值）为纵轴，直接展示"niche_high 相对 niche_low 有
    多少更多 spot 能检出该基因"，是火山图的重要补充视角。

    【颜色说明】
      🟠 橙色：仅 delta_frac > 阈值（检出率富集，但 FC 不高；典型稀有免疫基因）
      🔴 红色：同时满足 log2FC > 阈值 且 delta_frac > 阈值（双重显著）
      🔵 蓝色：仅 log2FC > 阈值（表达量倍数高，但检出率差异不大；可能是高表达管家基因）
      ⚫ 灰色：两者均不显著

    参数
    ----
    ranked               : _rank_niche_genes 返回的 DataFrame（含 log2_fc、delta_frac 列）
    path                 : 输出图片路径
    delta_frac_threshold : delta_frac 显著性阈值（默认 0.10，即检出率差 10 个百分点）
    fc_threshold         : log2FC 显著性阈值（默认 0.5）
    highlight_genes      : 需要特别标注名称的基因列表（如免疫抑制目标基因）
    """
    if "delta_frac" not in ranked.columns or ranked.empty:
        logging.warning("No delta_frac column found; skipping fraction scatter plot.")
        return

    if highlight_genes is None:
        highlight_genes = [
            g for gs in PRIOR_GENE_SETS.values() for g in gs
        ] + ["FOXP3", "TGFB1", "FAP", "ACTA2", "CCL22", "CXCL12", "SPP1"]

    df = ranked.copy()
    # 防止 log2_fc 列不存在时报错
    if "log2_fc" not in df.columns:
        logging.warning("No log2_fc column; skipping fraction scatter plot.")
        return

    def _color(row: pd.Series) -> str:
        has_frac = row["delta_frac"] > delta_frac_threshold
        has_fc   = row["log2_fc"]   > fc_threshold
        if has_fc and has_frac:
            return "#d62728"   # 红色：双重显著
        elif has_frac:
            return "#ff7f0e"   # 橙色：检出率富集
        elif has_fc:
            return "#4575b4"   # 蓝色：仅 FC 高
        else:
            return "#bdbdbd"   # 灰色

    c_list = [_color(row) for _, row in df.iterrows()]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df["log2_fc"], df["delta_frac"],
               c=c_list, s=14, alpha=0.7, linewidths=0)

    # 阈值线
    ax.axvline(fc_threshold,         color="black", linestyle="--", linewidth=0.8)
    ax.axhline(delta_frac_threshold, color="black", linestyle="--", linewidth=0.8)

    # 标注目标基因
    labeled = set()
    for _, row in df.iterrows():
        gene = row.get("gene", "")
        if gene in highlight_genes and gene not in labeled:
            ax.annotate(
                gene, (row["log2_fc"], row["delta_frac"]),
                textcoords="offset points", xytext=(4, 2),
                fontsize=7, color="#333333",
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.5),
            )
            labeled.add(gene)

    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d62728", label=f"Both: log2FC>{fc_threshold} & Δfrac>{delta_frac_threshold}"),
        Patch(facecolor="#ff7f0e", label=f"Fraction-only: Δfrac>{delta_frac_threshold}"),
        Patch(facecolor="#4575b4", label=f"FC-only: log2FC>{fc_threshold}"),
        Patch(facecolor="#bdbdbd", label="Not significant"),
    ]
    ax.legend(handles=legend_elements, fontsize=7, frameon=False, loc="upper left")

    n_frac = int((df["delta_frac"] > delta_frac_threshold).sum())
    n_both = int(((df["delta_frac"] > delta_frac_threshold) & (df["log2_fc"] > fc_threshold)).sum())
    ax.set_xlabel("log₂FC (niche_high / niche_low)")
    ax.set_ylabel(f"Δ Fraction (frac_high − frac_low)\n(proportion of spots with expression > 0)")
    ax.set_title(
        f"Niche Signature Genes: Fraction of Spots Enrichment\n"
        f"Δfrac-sig={n_frac}  Both-sig={n_both}  (|FC|>{fc_threshold}, Δfrac>{delta_frac_threshold})",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Fraction scatter plot saved: %s", path)


# ============================================================
# 新增：聚类 vs 评分并排对比图
# ============================================================

def _plot_niche_comparison(
    df: pd.DataFrame,
    path: Path,
    quantile_pct: int = 80,
) -> None:
    """
    绘制聚类方法 vs 评分方法的 niche 空间分布并排对比图。

    左图：Leiden 聚类注释的 immunosuppressive niche（niche_semantic_label 列，
          语义标签为 "immunosuppressive" 的 spot 高亮，其余用灰色）
    右图：评分阈值切割的 niche_high spot（niche_high 布尔列，True 为红色，False 为灰色）

    两图共用同一空间坐标系，用于评估两种方法的一致性：
      - 高度重叠：两种方法得到一致的结论，结果稳健
      - 聚类有但评分无：邻域组成像免疫抑制 niche 但自身评分不够高（边缘区域）
      - 评分有但聚类无：自身评分高但周围邻域不典型（孤立的免疫抑制岛）

    参数
    ----
    df          : 主分析 DataFrame（含 spatial_x/y、niche_semantic_label、niche_high 列）
    path        : 输出图片路径
    quantile_pct: niche_high 切割的分位数百分位（如 80 代表 80th percentile）
    """
    if "spatial_x" not in df.columns or "spatial_y" not in df.columns:
        logging.warning("spatial_x/y not found; skipping niche comparison plot.")
        return

    x = df["spatial_x"].to_numpy()
    y = df["spatial_y"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── 左图：Leiden 聚类方法 ────────────────────────────────────────────────
    ax0 = axes[0]
    if "niche_semantic_label" in df.columns:
        # 高亮含"immunosuppressive"或"Immunosuppressive"的语义标签 spot
        is_niche = df["niche_semantic_label"].astype(str).str.contains(
            "immunosuppressive", case=False, na=False
        )
        c0 = np.where(is_niche, "#d62728", "#bdbdbd")
        ax0.scatter(x[~is_niche], y[~is_niche], c="#bdbdbd", s=6, alpha=0.5, linewidths=0)
        ax0.scatter(x[is_niche],  y[is_niche],  c="#d62728", s=8, alpha=0.8, linewidths=0,
                    label="Immunosuppressive niche")
    else:
        ax0.scatter(x, y, c="#bdbdbd", s=6, alpha=0.5, linewidths=0)
        ax0.text(0.5, 0.5, "niche_semantic_label\nnot available",
                 transform=ax0.transAxes, ha="center", va="center", fontsize=9)

    ax0.set_title("Method 1: Leiden Cluster Annotation\n(immunosuppressive semantic label)",
                  fontsize=9, fontweight="bold")
    ax0.set_aspect("equal")
    ax0.axis("off")
    ax0.legend(loc="lower right", fontsize=7, frameon=False, markerscale=1.5)

    # ── 右图：评分阈值方法 ───────────────────────────────────────────────────
    ax1 = axes[1]
    if "niche_high" in df.columns:
        niche_high = df["niche_high"].astype(bool).to_numpy()
        ax1.scatter(x[~niche_high], y[~niche_high], c="#bdbdbd", s=6, alpha=0.5, linewidths=0)
        ax1.scatter(x[niche_high],  y[niche_high],  c="#d62728", s=8, alpha=0.8, linewidths=0,
                    label=f"niche_high (score >{quantile_pct}th pct)")
    else:
        ax1.scatter(x, y, c="#bdbdbd", s=6, alpha=0.5, linewidths=0)
        ax1.text(0.5, 0.5, "niche_high\nnot available",
                 transform=ax1.transAxes, ha="center", va="center", fontsize=9)

    ax1.set_title(f"Method 2: Score-Based Thresholding\n(niche score > {quantile_pct}th percentile)",
                  fontsize=9, fontweight="bold")
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.legend(loc="lower right", fontsize=7, frameon=False, markerscale=1.5)

    # 计算两方法的 Jaccard 相似度（如果两列都存在）
    jaccard_str = ""
    if "niche_semantic_label" in df.columns and "niche_high" in df.columns:
        is_niche_arr = df["niche_semantic_label"].astype(str).str.contains(
            "immunosuppressive", case=False, na=False
        ).to_numpy()
        niche_high_arr = df["niche_high"].astype(bool).to_numpy()
        inter = int((is_niche_arr & niche_high_arr).sum())
        union = int((is_niche_arr | niche_high_arr).sum())
        if union > 0:
            jaccard = inter / union
            jaccard_str = f"  Jaccard similarity: {jaccard:.3f}"

    fig.suptitle(
        f"Immunosuppressive Niche: Leiden Clustering vs Score Thresholding\n"
        f"(Left = cluster-based; Right = score-based{jaccard_str})",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Niche comparison plot saved: %s", path)


# ============================================================
# 新增：参数扫描（k × niche_high_quantile 网格搜索）
# ============================================================

def _param_scan_deg_stability(
    adata: ad.AnnData,
    proportions: pd.DataFrame,
    coords: np.ndarray,
    gene_score: pd.Series,
    treg_col: str,
    myeloid_col: str,
    fibroblast_col: str,
    n_neighbors: int = 15,
    k_list: tuple[int, ...] = (8, 10, 15, 20, 25),
    quantile_list: tuple[float, ...] = (0.70, 0.75, 0.80, 0.85),
    top_n: int = 50,
    seed: int = 1234,
) -> pd.DataFrame:
    """
    对 kNN 邻居数（k）× niche_high 分位数阈值进行网格扫描，
    以 DEG 稳定性（显著基因数量 + 相邻参数组合 Jaccard 相似度）
    作为目标函数，辅助选择最优参数组合。

    【为什么改为 k × quantile，而不是 resolution × quantile？】
    原来的 resolution × quantile 扫描存在根本性缺陷：
      - niche_score 的计算（Step 9）基于邻域组成向量的 Z-score 加总；
      - 邻域组成向量由 k-NN 邻居数决定，k 不同则每个 spot 的邻域大小不同；
      - Leiden resolution 仅决定聚类粒度（n_clusters），不影响 niche_score 数值；
      - 因此，固定 k 时改变 resolution 对 DEG 结果没有任何影响，热图列完全相同。

    k 才是真正影响 niche 评分的参数：
      - 小 k（如 k=8）：邻域小，捕捉微观局灶性免疫聚集；
      - 大 k（如 k=25）：邻域大，捕捉宏观区域性免疫浸润模式；
      - 不同 k 下 niche_score 向量不同，DEG 结果自然不同，热图有意义差异。

    扫描逻辑（k × quantile 共 5×4 = 20 种参数组合）：
      1. 对每种 k，重新构建 kNN 邻接图，重新计算邻域组成向量和 niche_score；
      2. 对每种 quantile 阈值切割 niche_high；
      3. 对 niche_high vs niche_low 做 Wilcoxon 检验，记录：
           - n_sig_deg      : FDR < 0.05 且 log2FC > 0.5 的基因数量
           - mean_log2fc    : Top-N 基因的平均 log2FC
           - top_genes_str  : 前 20 个基因名（用于跨参数 Jaccard 比较）
      4. 对相邻参数组合（k ±1档 或 quantile ±0.05），
         计算 Top-N 基因列表之间的 Jaccard 相似度（稳健性指标）

    【热图解读】
      - 横轴：niche_high_quantile（阈值越高 = niche_high 越少 = 越严格）
      - 纵轴：k（邻居数越大 = 邻域越大 = 捕捉宏观模式）
      - 颜色：n_sig_deg 越深 = 该参数组合下可重复 DEG 越多
      - 选参标准：颜色最深且处于"高原"区域（与相邻参数结果相近）的组合即为推荐

    参数
    ----
    adata          : AnnData 数据对象
    proportions    : 细胞类型比例矩阵
    coords         : 空间坐标
    gene_score     : 免疫抑制基因模块评分
    treg_col       : Treg 比例列名
    myeloid_col    : Myeloid 比例列名
    fibroblast_col : Fibroblast 比例列名
    n_neighbors    : 主流程使用的 kNN 邻居数（参考值，出现在热图中对应格子）
    k_list         : kNN 邻居数扫描列表（主轴参数，真正影响 niche_score）
    quantile_list  : niche_high 分位数阈值扫描列表
    top_n          : 用于计算 Jaccard 的 Top-N 基因数
    seed           : 随机数种子（供 Leiden 使用，本步骤不重新运行 Leiden）

    返回
    ----
    pd.DataFrame，每行为一种参数组合的结果：
      k / quantile / n_sig_deg / mean_log2fc_topN / top_genes_str
    """
    import itertools

    gs = gene_score.reindex(proportions.index).fillna(0.0)
    expr = _expression_frame(adata)
    expr_columns = expr.columns.tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(expr_columns)}
    rows = []

    for k, quantile in itertools.product(k_list, quantile_list):
        try:
            # 1. 重新构建 k 近邻图（k 变化 → 邻域大小变化 → niche_score 真正不同）
            neighbors_k = _build_knn_neighbors(coords, k=k)

            # 2. 用新的邻域大小重新计算邻域组成向量和 niche_score
            neighborhood_comp = _compute_neighborhood_composition(proportions, neighbors_k)
            niche_score = (
                _zscore(neighborhood_comp[treg_col])
                + _zscore(neighborhood_comp[myeloid_col])
                + _zscore(neighborhood_comp[fibroblast_col])
                + _zscore(gs)
            )

            # 3. 切割 niche_high（quantile 变化 → 高分组大小变化）
            thr = float(niche_score.quantile(quantile))
            niche_high = niche_score >= thr

            if int(niche_high.sum()) < 3 or int((~niche_high).sum()) < 3:
                logging.debug(
                    "Param scan [k=%d, q=%.2f]: skipped (too few spots).", k, quantile
                )
                continue

            # 4. Wilcoxon 检验（快速版：仅对 log2FC > 0.3 的基因检验）
            # niche_score 的索引来自 proportions，expr 的索引来自 adata.obs_names，
            # 两者可能不完全一致，必须对齐后再用布尔索引，否则触发
            # "Unalignable boolean Series" 错误。
            niche_high_aligned = niche_high.reindex(expr.index).fillna(False).astype(bool)
            high_expr = expr.loc[niche_high_aligned]
            low_expr  = expr.loc[~niche_high_aligned]
            mean_h = high_expr.mean()
            mean_l = low_expr.mean()
            lfc = np.log2((mean_h + 1.0) / (mean_l + 1.0))

            cands = expr.columns[lfc.to_numpy() > 0.3].tolist()
            pvals = np.ones(len(expr.columns))
            for gene in cands:
                gidx = gene_to_idx[gene]
                try:
                    _, p = scipy_stats.mannwhitneyu(
                        high_expr[gene].to_numpy(), low_expr[gene].to_numpy(),
                        alternative="greater",
                    )
                    pvals[gidx] = p
                except Exception:
                    pass
            _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
            sig_mask = (fdr < 0.05) & (lfc.to_numpy() > 0.5)
            n_sig = int(sig_mask.sum())

            # Top-N 基因（按 log2FC 降序）
            sig_df = pd.DataFrame({
                "gene":    expr_columns,
                "log2_fc": lfc.to_numpy(),
                "fdr":     fdr,
            })
            top_genes = (
                sig_df[sig_df["fdr"] < 0.05]
                .sort_values("log2_fc", ascending=False)
                .head(top_n)["gene"]
                .tolist()
            )
            if not top_genes:
                top_genes = sig_df.sort_values("log2_fc", ascending=False).head(top_n)["gene"].tolist()

            mean_lfc_top = float(lfc.reindex(top_genes).dropna().mean()) if top_genes else 0.0

            rows.append({
                "k":               k,
                "quantile":        quantile,
                "n_sig_deg":       n_sig,
                "mean_log2fc_topN": round(mean_lfc_top, 4),
                "top_genes_str":   ";".join(top_genes[:20]),
            })
            logging.info(
                "Param scan [k=%d, q=%.2f]: n_sig_deg=%d, mean_log2fc_top%d=%.3f",
                k, quantile, n_sig, top_n, mean_lfc_top,
            )
        except Exception as exc:
            logging.warning(
                "Param scan [k=%d, q=%.2f] failed: %s", k, quantile, exc
            )

    return pd.DataFrame(rows)


def _plot_param_scan_heatmap(
    param_scan_df: pd.DataFrame,
    path: Path,
) -> None:
    """
    绘制参数扫描热图：横轴为 niche_high_quantile，纵轴为 Leiden resolution，
    颜色编码为显著 DEG 数量（n_sig_deg）。

    热图用于可视化最优参数区域：颜色最深（DEG 数量最多）且处于"高原"区域
    （与相邻参数结果相近）的参数组合即为推荐选择。

    同时在每个格子中标注 n_sig_deg 数值，便于直接读取。

    参数
    ----
    param_scan_df : _param_scan_deg_stability 返回的 DataFrame
    path          : 输出图片路径
    """
    if param_scan_df.empty or "n_sig_deg" not in param_scan_df.columns:
        logging.warning("Empty param scan DataFrame; skipping heatmap.")
        return

    # 透视表：行=k（kNN邻居数），列=quantile，值=n_sig_deg
    pivot = param_scan_df.pivot(
        index="k", columns="quantile", values="n_sig_deg"
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(
        pivot,
        cmap="YlOrRd",
        annot=True, fmt="d",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Number of significant DEGs\n(FDR<0.05, log2FC>0.5)"},
    )
    ax.set_xlabel("niche_high quantile threshold")
    ax.set_ylabel("kNN neighbors (k)")
    ax.set_title(
        "Parameter Scan: DEG Stability Heatmap\n"
        "(k × quantile; darker = more stable DEGs; optimal = darkest plateau)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    logging.info("Parameter scan heatmap saved: %s", path)


# ============================================================
# 主流程
# ============================================================

def main() -> int:
    """
    主函数：执行完整的空间免疫抑制生态位分析流程。

    返回 0 表示正常完成（供 Shell 脚本通过 $? 检查）。
    """
    # ── Step 1: 初始化 ──────────────────────────────────────────────────────────
    _setup_logging()
    args = _parse_args()
    np.random.seed(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 2: 数据加载与校验 ───────────────────────────────────────────────────
    logging.info("Loading spatial AnnData: %s", args.adata)
    adata = ad.read_h5ad(args.adata)
    if "spatial" not in adata.obsm:
        raise KeyError("adata.obsm['spatial'] is required for spatial niche analysis.")

    abundance = _get_abundance(adata, args.abundance_key)
    required_cols = [
        args.hepatocyte_col, args.treg_col,
        args.myeloid_col, args.fibroblast_col, args.tnk_col,
    ]
    missing_cols = [c for c in required_cols if c not in abundance.columns]
    if missing_cols:
        raise KeyError(
            f"Missing cell abundance columns: {missing_cols}; "
            f"available: {list(abundance.columns)}"
        )

    # ── Step 3: 归一化细胞丰度为比例 ────────────────────────────────────────────
    proportions = (
        abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0)
        .fillna(0.0)
    )
    coords = np.asarray(adata.obsm["spatial"])

    # ── Step 4: 构建空间半径邻域（用于计算邻域均值特征） ─────────────────────────
    neighbors_radius, radius = _spatial_neighbors(coords, args.neighbor_radius_multiplier)
    logging.info(
        "Spatial neighbor radius (multiplier=%.2f): %.3f",
        args.neighbor_radius_multiplier, radius,
    )

    # ── Step 5: 计算免疫抑制基因模块评分 ────────────────────────────────────────
    gene_score, marker_genes = _module_score(adata, IMMUNOSUPPRESSIVE_GENES)
    logging.info(
        "Immunosuppressive marker genes used (%d): %s",
        len(marker_genes),
        ", ".join(marker_genes) if marker_genes else "none",
    )

    # ── Step 6: 邻域组成聚类（Cellular Neighborhood Discovery）────────────────
    logging.info(
        "Step 6: Building kNN neighbors (k=%d) for neighborhood composition clustering...",
        args.n_neighbors,
    )
    neighbors_knn = _build_knn_neighbors(coords, k=args.n_neighbors)
    neighborhood_comp = _compute_neighborhood_composition(proportions, neighbors_knn)

    logging.info(
        "Step 6: Running Leiden clustering on neighborhood composition "
        "(resolution=%.2f, seed=%d)...",
        args.leiden_resolution, args.seed,
    )
    nc_labels = _leiden_cluster_neighborhood(
        neighborhood_comp,
        n_neighbors=args.n_neighbors,
        resolution=args.leiden_resolution,
        seed=args.seed,
    )
    n_clusters = len(nc_labels.unique())
    logging.info("Leiden clustering identified %d neighborhood clusters.", n_clusters)

    # ── Step 7: 功能评分注释 niche cluster ──────────────────────────────────────
    logging.info("Step 7: Annotating neighborhood clusters with functional scores...")
    tmp_df = proportions.copy()
    tmp_df["immunosuppressive_gene_score"] = (
        gene_score.reindex(tmp_df.index).fillna(0.0).to_numpy()
    )
    tmp_df["neighborhood_cluster"] = nc_labels.reindex(tmp_df.index).to_numpy()

    niche_semantic, cluster_stats = _annotate_neighborhood_clusters(
        df=tmp_df,
        cluster_col="neighborhood_cluster",
        treg_col=args.treg_col,
        myeloid_col=args.myeloid_col,
        fibroblast_col=args.fibroblast_col,
    )
    cluster_stats.to_csv(args.out_dir / "neighborhood_cluster_stats.csv")

    # ── Step 8: 构建综合分析 DataFrame ──────────────────────────────────────────
    df = proportions.copy()
    df.insert(0, "spot_id", adata.obs_names)
    df["spatial_x"]   = coords[:, 0]
    df["spatial_y"]   = coords[:, 1]
    df["neighbor_radius"] = radius
    df["immunosuppressive_gene_score"] = gene_score.reindex(df.index).to_numpy()
    df["neighborhood_cluster"]   = nc_labels.reindex(df.index).to_numpy()
    df["niche_semantic_label"]   = niche_semantic.reindex(df.index).to_numpy()

    # 高肝细胞区域标志（辅助注释，不用于 niche 发现主体）
    hep = df[args.hepatocyte_col]
    hep_thr = float(hep.quantile(args.hep_high_quantile))
    df["hep_high"] = hep >= hep_thr
    df["distance_to_hep_high"] = _distance_to_mask(coords, df["hep_high"].to_numpy())

    # 各细胞类型的空间邻域均值
    for col in [args.hepatocyte_col, "Treg", args.tnk_col,
                args.myeloid_col, args.fibroblast_col]:
        out_col = col.replace("/", "_")
        df[f"neighbor_{out_col}"] = _neighbor_mean(df[col], neighbors_radius)

    # 是否有高肝细胞邻居（用于 tumor_edge 标注）
    hep_high_arr = df["hep_high"].to_numpy()
    df["has_hep_high_neighbor"] = [
        bool(len(idx) and hep_high_arr[idx].any())
        for idx in neighbors_radius
    ]

    # ── Step 9: 多层次评分计算 ──────────────────────────────────────────────────
    # Treg 样评分：Treg 比例 + 免疫抑制基因评分（Z-score 加总）
    df["Treg_like_score"] = (
        _zscore(df["Treg"])
        + _zscore(df["immunosuppressive_gene_score"])
    )
    # 免疫基质评分：Treg + T/NK + Myeloid + Fibroblast 的 Z-score 和
    df["immune_stroma_score"] = (
        _zscore(df["Treg"])
        + _zscore(df[args.tnk_col])
        + _zscore(df[args.myeloid_col])
        + _zscore(df[args.fibroblast_col])
    )
    # 综合免疫抑制 niche 评分（用于签名提取阈值的参考评分，不用于 niche 发现主体）
    df["immunosuppressive_niche_score"] = (
        _zscore(df["Treg"])
        + _zscore(df[args.myeloid_col])
        + _zscore(df[args.fibroblast_col])
        + _zscore(df["immunosuppressive_gene_score"])
        + _zscore(df["neighbor_Hepatocyte"].fillna(df[args.hepatocyte_col]))
    )
    niche_thr = float(
        df["immunosuppressive_niche_score"].quantile(args.niche_high_quantile)
    )
    # niche_high：用于签名基因提取分组（基于综合评分分位数）
    df["niche_high"] = df["immunosuppressive_niche_score"] >= niche_thr

    # ── Step 10: 辅助空间区域标注（基于规则，仅用于可视化） ─────────────────────
    stroma_thr = float(df["immune_stroma_score"].quantile(0.60))
    region = pd.Series("other", index=df.index)
    region[df["hep_high"]] = "tumor_core"
    region[
        (~df["hep_high"]) & (df["immune_stroma_score"] >= stroma_thr)
    ] = "stroma_immune"
    region[
        (~df["hep_high"]) & df["has_hep_high_neighbor"]
        & (df["immune_stroma_score"] >= stroma_thr)
    ] = "tumor_edge"
    df["spatial_region"] = region

    # ── Step 11: 敏感性分析 ──────────────────────────────────────────────────────
    logging.info("Step 11: Running sensitivity analysis...")
    sensitivity_df = _sensitivity_analysis(
        proportions=proportions,
        coords=coords,
        gene_score=gene_score,
        treg_col=args.treg_col,
        myeloid_col=args.myeloid_col,
        fibroblast_col=args.fibroblast_col,
    )
    sensitivity_df.to_csv(args.out_dir / "sensitivity_analysis.csv", index=False)
    logging.info("Sensitivity analysis saved.")

    # ── Step 12: 保存结果 ─────────────────────────────────────────────────────────
    out_csv = args.out_dir / "spatial_niche_scores.csv"
    df.to_csv(out_csv, index=False)
    logging.info("Niche score table saved: %s", out_csv)

    metadata = pd.DataFrame(
        [
            ("abundance_key",         args.abundance_key),
            ("n_neighbors_knn",       args.n_neighbors),
            ("leiden_resolution",     args.leiden_resolution),
            ("n_neighborhood_clusters", n_clusters),
            ("neighbor_radius_multiplier", args.neighbor_radius_multiplier),
            ("neighbor_radius",       radius),
            ("hep_high_quantile",     args.hep_high_quantile),
            ("hep_high_threshold",    hep_thr),
            ("niche_high_quantile",   args.niche_high_quantile),
            ("niche_high_threshold",  niche_thr),
            ("n_niche_high",          int(df["niche_high"].sum())),
            ("marker_genes_used",     ",".join(marker_genes)),
        ],
        columns=["parameter", "value"],
    )
    metadata.to_csv(args.out_dir / "spatial_niche_parameters.csv", index=False)

    # ── Step 13: 可视化 ──────────────────────────────────────────────────────────
    logging.info("Step 13: Generating plots...")
    # 各细胞类型空间分布图
    _spatial_scatter(df, args.hepatocyte_col,
                     plot_dir / "spatial_hepatocyte.png", "Hepatocyte proportion")
    _spatial_scatter(df, "Treg",
                     plot_dir / "spatial_treg.png", "Treg-like abundance")
    _spatial_scatter(df, args.myeloid_col,
                     plot_dir / "spatial_myeloid.png", "Myeloid proportion")
    _spatial_scatter(df, args.fibroblast_col,
                     plot_dir / "spatial_fibroblast.png", "Fibroblast proportion")
    # 免疫抑制 niche 综合评分空间分布图（magma 配色）
    _spatial_scatter(
        df, "immunosuppressive_niche_score",
        plot_dir / "spatial_immunosuppressive_niche_score.png",
        "Immunosuppressive niche score",
        cmap="magma",
    )
    # 邻域组成聚类空间分布图
    _plot_neighborhood_clusters(df, plot_dir / "spatial_neighborhood_clusters.png")
    # 语义 niche 标签空间分布图
    _spatial_scatter(
        df, "niche_semantic_label",
        plot_dir / "spatial_niche_semantic_labels.png",
        "Immunosuppressive Niche (Leiden Clusters)",
        categorical=True,
    )
    # 辅助区域标注图（规则化）
    _spatial_scatter(
        df, "spatial_region",
        plot_dir / "spatial_region_labels.png",
        "Spatial region labels (rule-based annotation)",
        categorical=True,
    )
    # 距离依赖性折线图
    _plot_distance(df, plot_dir / "distance_to_hep_high_vs_niche_score.png")
    # 各区域评分箱线图
    _plot_region_box(df, plot_dir / "region_score_boxplots.png")
    # Hepatocyte vs Treg 散点图
    _plot_hep_treg(df, plot_dir / "hepatocyte_vs_treg_niche_score.png")
    # 细胞类型相关性热图
    corr_cols = [
        args.hepatocyte_col, "Treg", args.tnk_col,
        args.myeloid_col, args.fibroblast_col,
        "Treg_like_score", "immunosuppressive_niche_score",
    ]
    _plot_correlation(df, corr_cols, plot_dir / "celltype_niche_correlation.png")
    # L-R 通讯热图（新增）
    _plot_lr_communication(
        adata, df, LR_PAIRS,
        plot_dir / "lr_communication_heatmap.png",
    )
    # 敏感性分析结果图（新增）
    _plot_sensitivity(sensitivity_df, plot_dir / "sensitivity_niche_stability.png")

    # ── 新增图1：niche_high 分位数阈值切割的 spot 空间二值分布图 ─────────────────
    # 颜色编码：红色 = 评分高于 niche_high_quantile（如 80 百分位）的 spot
    #           灰色 = 其余 spot
    # 与 spatial_niche_semantic_labels.png（Leiden 聚类方法）对比，
    # 展示"评分法"与"聚类法"在空间上的差异（重叠/分歧）。
    logging.info("Step 13+: Generating niche_high score-based binary map...")
    _spatial_scatter(
        df, "niche_high",
        plot_dir / "spatial_niche_high_score_spots.png",
        f"Niche-High Spots (score > {int(args.niche_high_quantile * 100)}th percentile)",
        cmap="RdGy_r",
    )

    # ── 新增图2：聚类方法 vs 评分方法并排对比图 ──────────────────────────────────
    # 左图：Leiden 聚类注释的 immunosuppressive niche（niche_semantic_label）
    # 右图：评分阈值切割的 niche_high spot 分布（niche_high）
    # 两图共用同一空间坐标系，便于直接对比空间位置
    logging.info("Step 13+: Generating clustering vs score comparison plot...")
    _plot_niche_comparison(
        df,
        plot_dir / "spatial_niche_cluster_vs_score_comparison.png",
        quantile_pct=int(args.niche_high_quantile * 100),
    )

    # ── Step 14: 特征基因提取（升级版） ─────────────────────────────────────────
    logging.info("Step 14: Extracting niche signature genes (with Wilcoxon + FDR + Gini)...")
    ranked = _rank_niche_genes(adata, df["niche_high"], args.top_niche_genes)

    # 保存 Layer 1：主签名基因（Wilcoxon + FDR）
    ranked.to_csv(
        args.out_dir / "immunosuppressive_niche_signature_genes_ranked.csv",
        index=False,
    )
    sig_path = args.out_dir / "immunosuppressive_niche_signature_genes.txt"
    sig_text = "\n".join(ranked["gene"].tolist()) + "\n"
    sig_path.write_text(sig_text, encoding="utf-8")

    args.signature_out.parent.mkdir(parents=True, exist_ok=True)
    args.signature_out.write_text(sig_text, encoding="utf-8")

    # 保存 Layer 2：先验功能基因集 AUC 检验结果
    prior_auc_df = ranked.attrs.get("prior_auc_df", pd.DataFrame())
    if not prior_auc_df.empty:
        prior_auc_df.to_csv(
            args.out_dir / "prior_gene_set_auc.csv", index=False
        )
        logging.info(
            "Prior gene set AUC results saved: %s",
            args.out_dir / "prior_gene_set_auc.csv",
        )

    # 保存 Layer 3：Gini Index 特异性评分
    # 无论结果是否为空，始终写出 CSV（保证文件路径可预期）
    gini_df = ranked.attrs.get("gini_df", pd.DataFrame(
        columns=["gene", "gini", "log2_fc", "mean_high"]
    ))
    gini_out = args.out_dir / "gini_score_genes.csv"
    gini_df.to_csv(gini_out, index=False)
    logging.info(
        "Gini index genes saved: %s  (%d genes with Gini>0.3 & log2FC>0)",
        gini_out,
        len(gini_df),
    )

    # 绘制火山图 + 检出率散点图（基于全基因范围的统计量）
    # 注：此处需传入包含全部基因统计量的 DataFrame（ranked 仅含 Top-N），
    # 因此重新计算全基因范围的 log2FC / delta_frac / pvalue / FDR，
    # 供火山图和检出率散点图共享，无需对同一数据集做两次独立计算。
    logging.info("Step 14: Generating volcano plot + fraction scatter for all tested genes...")
    try:
        from scipy.stats import mannwhitneyu as _mwu
        expr_full = _expression_frame(adata)
        high_full = df["niche_high"].astype(bool).reindex(expr_full.index).fillna(False)

        # ── 全基因均值和 log2FC ──────────────────────────────────────────────
        mean_h = expr_full.loc[high_full].mean()
        mean_l = expr_full.loc[~high_full].mean()
        lfc_all = np.log2((mean_h + 1.0) / (mean_l + 1.0))

        # ── 全基因检出率（fraction of spots with expression > 0） ───────────
        # 解决稀疏基因（如 FOXP3）在均值计算中被稀释的问题
        frac_h_all = (expr_full.loc[high_full]  > 0).mean(axis=0)
        frac_l_all = (expr_full.loc[~high_full] > 0).mean(axis=0)
        delta_frac_all = frac_h_all - frac_l_all

        # ── 全基因 Wilcoxon 检验（仅对候选基因，减少计算量） ────────────────
        # 候选条件：log2FC > 0.05 OR delta_frac > 0.03（覆盖稀疏目标基因）
        pvals_all = np.full(len(expr_full.columns), 1.0)
        cand_mask = (lfc_all.to_numpy() > 0.05) | (delta_frac_all.to_numpy() > 0.03)
        cands = expr_full.columns[cand_mask].tolist()
        for gene in cands:
            gidx = list(expr_full.columns).index(gene)
            try:
                _, p = _mwu(
                    expr_full.loc[high_full, gene].to_numpy(),
                    expr_full.loc[~high_full, gene].to_numpy(),
                    alternative="greater",
                )
                pvals_all[gidx] = p
            except Exception:
                pass
        _, fdr_all, _, _ = multipletests(pvals_all, method="fdr_bh")

        # ── 构建全基因 DataFrame（同时包含 FC 和 Fraction 信息）─────────────
        volcano_df = pd.DataFrame({
            "gene":       expr_full.columns.tolist(),
            "log2_fc":    lfc_all.to_numpy(),
            "frac_high":  frac_h_all.to_numpy(),
            "frac_low":   frac_l_all.to_numpy(),
            "delta_frac": delta_frac_all.to_numpy(),
            "fdr":        fdr_all,
        })

        # 火山图（log2FC vs -log10 FDR）
        _plot_volcano(volcano_df, plot_dir / "niche_signature_volcano.png")

        # 检出率散点图（log2FC vs delta_frac，突出 FOXP3 等稀有基因的富集）
        logging.info("Step 14: Generating fraction scatter plot...")
        _plot_fraction_scatter(
            volcano_df,
            plot_dir / "niche_fraction_scatter.png",
        )
    except Exception as exc:
        logging.warning("Volcano/fraction plot failed: %s", exc)

    # ── Step 15: 参数扫描（resolution × niche_high_quantile 网格搜索）────────────
    logging.info("Step 15: Running parameter scan (k × quantile grid)...")
    try:
        param_scan_df = _param_scan_deg_stability(
            adata=adata,
            proportions=proportions,
            coords=coords,
            gene_score=gene_score,
            treg_col=args.treg_col,
            myeloid_col=args.myeloid_col,
            fibroblast_col=args.fibroblast_col,
            n_neighbors=args.n_neighbors,
            k_list=(8, 10, 15, 20, 25),
            quantile_list=(0.70, 0.75, 0.80, 0.85),
            top_n=args.top_niche_genes,
            seed=args.seed,
        )
        param_scan_df.to_csv(
            args.out_dir / "param_scan_deg_stability.csv", index=False
        )
        _plot_param_scan_heatmap(
            param_scan_df,
            plot_dir / "param_scan_deg_stability_heatmap.png",
        )
        logging.info("Parameter scan completed. Results saved.")
    except Exception as exc:
        logging.warning("Parameter scan failed: %s", exc)

    logging.info("Niche signature genes saved: %s", sig_path)
    logging.info("TCGA-ready signature also saved: %s", args.signature_out)
    logging.info(
        "Done. Next: run DE validation with:\n"
        "  python code/run_de_analysis.py --coloc %s --coloc-column niche_high",
        out_csv,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())