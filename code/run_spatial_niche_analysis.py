"""
================================================================================
脚本名称: run_spatial_niche_analysis.py
功能概述: 空间免疫抑制生态位（Niche）分析 —— Step 2
================================================================================

【整体任务说明】
    本脚本基于 Cell2location 反卷积结果，采用"邻域组成聚类 + 功能评分注释"
    的两阶段策略识别肝癌空间免疫抑制生态位（Schürch et al., Cell, 2020）。
    具体流程：
      1. 构建空间邻接图，计算每个 spot 的邻域细胞组成向量（Squidpy 框架）；
      2. 对邻域组成向量进行 Leiden 无监督聚类，识别 cellular neighborhoods；
      3. 对每个 neighborhood cluster 计算功能评分（Treg_score、免疫抑制基因
         模块评分等），将高 Treg、高免疫基质成分的 cluster 注释为目标 niche；
      4. 对阈值参数进行敏感性分析，验证结果稳定性；
      5. 提取 niche 特征基因签名，生成多类可视化图表（含 L-R 通讯热图）。

【niche 识别策略说明（参考 niche修改.md）】
    - 本脚本将阈值规则（hep_high_quantile、niche_high_quantile 等）降级为
      "注释/筛选"步骤，而非 niche 发现的第一步；
    - niche 发现主体改为：空间邻接图 → 邻域组成向量 → Leiden 聚类；
    - 在 Leiden 聚类结果上叠加功能评分，对 niche 进行语义注释；
    - 通过敏感性分析（k=10/15/20 邻居，radius_multiplier=1.0/1.25/1.5）
      验证 niche 分配的稳定性。

【参考文献】
    - Schürch et al., Cell, 2020 (cellular neighborhoods 核心方法)
    - Palla et al., Nature Methods, 2022 (Squidpy 空间分析框架)
    - Keren et al., Cell, 2018 (肿瘤-免疫空间结构分型)

【输入文件】
    adata_vis_post.h5ad   - CHC20 主分析的 Cell2location 反卷积空间数据
                            （CHC23 为独立验证输出，默认不进入 niche 发现）

【输出文件】
    spatial_niche/
      spatial_niche_scores.csv                         - 完整评分表
      spatial_niche_parameters.csv                     - 分析参数与元数据
      sensitivity_analysis.csv                         - 敏感性分析结果
      neighborhood_cluster_stats.csv                   - 邻域聚类统计
      immunosuppressive_niche_signature_genes_ranked.csv - 排名后签名基因表
      immunosuppressive_niche_signature_genes.txt      - 签名基因列表
      plots/                                           - 各类可视化图表
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
from scipy.spatial import KDTree


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
# ============================================================
LR_PAIRS: list[tuple[str, str, str]] = [
    ("CXCL12", "CXCR4",   "CXCL12–CXCR4"),
    ("CCL22",  "CCR4",    "CCL22–CCR4"),
    ("TGFB1",  "TGFBR1",  "TGFB1–TGFBR1"),
    ("PDCD1",  "CD274",   "PD-1–PD-L1"),
    ("IL10",   "IL10RA",  "IL10–IL10RA"),
    ("TIGIT",  "NECTIN2", "TIGIT–NECTIN2"),
    ("LAG3",   "HLA-DRA", "LAG3–MHC-II"),
]


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
) -> pd.DataFrame:
    """
    对邻域参数进行敏感性分析，评估 niche 高分区域占比的稳定性。

    对 k-NN 邻居数（k=10/15/20）和半径倍增系数（×1.0/1.25/1.5）
    分别进行扫描，记录每种设置下 niche_high spot 数量及占比，
    验证 niche 分配对参数选择不敏感（即结果稳定）。

    参数
    ----
    proportions             : 细胞类型比例矩阵
    coords                  : 空间坐标
    gene_score              : 免疫抑制基因模块评分
    n_neighbors_list        : k-NN 邻居数扫描列表
    radius_multiplier_list  : 半径倍增系数扫描列表
    niche_high_quantile     : 高 niche 区域分位数阈值

    返回
    ----
    pd.DataFrame，每行为一个参数组合的统计结果。
    """
    gs = gene_score.reindex(proportions.index).fillna(0.0)
    rows = []
    all_params: list[tuple[str, int | None, float | None]] = (
        [("knn", k, None) for k in n_neighbors_list]
        + [("radius", None, r) for r in radius_multiplier_list]
    )
    for mode, k, r in all_params:
        try:
            if mode == "knn":
                nbrs = _build_knn_neighbors(coords, k=k)
                label = f"kNN k={k}"
            else:
                nbrs, _ = _spatial_neighbors(coords, radius_multiplier=r)
                label = f"radius×{r}"

            nc = _compute_neighborhood_composition(proportions, nbrs)
            score = (
                _zscore(nc[treg_col])
                + _zscore(nc[myeloid_col])
                + _zscore(nc[fibroblast_col])
                + _zscore(gs)
            )
            thr = float(score.quantile(niche_high_quantile))
            n_high = int((score >= thr).sum())
            pct = n_high / len(score) * 100
            rows.append({
                "param_mode":  mode,
                "param_label": label,
                "n_niche_high": n_high,
                "niche_pct":   pct,
                "score_mean":  float(score.mean()),
                "score_std":   float(score.std()),
            })
            logging.info("Sensitivity [%s]: n_niche_high=%d (%.1f%%)", label, n_high, pct)
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


def _plot_lr_communication(
    adata: ad.AnnData,
    df: pd.DataFrame,
    lr_pairs: list[tuple[str, str, str]],
    path: Path,
) -> None:
    """
    绘制配体-受体（L-R）通讯分析热图，比较 niche_high vs niche_low 的信号强度。

    方法：
      - 对每个 L-R 对，计算配体与受体基因表达乘积的均值（product score），
        作为该位置通讯强度的代理指标；
      - 左图：归一化信号强度热图（行=L-R 对，列=niche 分组）；
      - 右图：log2FC（niche_high / niche_low）条形图，红色=上调，蓝色=下调。

    综述依据：CXCL12-CXCR4、CCL22-CCR4 是 Treg 招募的关键通讯轴
    （Efremova et al., 2020; Cang & Nie, 2023）。
    """
    expr = _expression_frame(adata)
    high_mask = df["niche_high"].astype(bool).to_numpy()
    low_mask  = ~high_mask

    records = []
    for ligand, receptor, label in lr_pairs:
        missing = [g for g in (ligand, receptor) if g not in expr.columns]
        if missing:
            logging.warning("L-R [%s]: genes %s not found; skipping.", label, missing)
            continue
        product   = expr[ligand].to_numpy() * expr[receptor].to_numpy()
        mean_high = float(product[high_mask].mean()) if high_mask.any() else 0.0
        mean_low  = float(product[low_mask].mean())  if low_mask.any()  else 0.0
        log2fc    = float(np.log2((mean_high + 1e-6) / (mean_low + 1e-6)))
        records.append({
            "lr_pair":   label,
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
    绘制敏感性分析结果图：不同参数设置下 niche-high spot 占比的变化。

    左图：k-NN 邻居数变化；右图：半径倍增系数变化。
    占比变化幅度小说明 niche 分配对参数选择稳健。
    """
    if sensitivity_df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, mode, title_str in zip(
        axes,
        ["knn",   "radius"],
        ["kNN Neighbor Count (k) Sensitivity",
         "Radius Multiplier Sensitivity"],
    ):
        sub = sensitivity_df[sensitivity_df["param_mode"] == mode]
        if sub.empty:
            ax.set_visible(False)
            continue
        ax.bar(sub["param_label"], sub["niche_pct"],
               color="#5b9bd5", edgecolor="none")
        ax.set_xlabel("Parameter Setting")
        ax.set_ylabel("Niche-High Spots (%)")
        ax.set_title(title_str)
        ax.tick_params(axis="x", rotation=20)
        # 标注参考水平线（主分析结果）
        main_label = "kNN k=15" if mode == "knn" else "radius×1.25"
        main_row = sub[sub["param_label"] == main_label]
        if not main_row.empty:
            ax.axhline(
                float(main_row["niche_pct"].iloc[0]),
                color="red", linestyle="--", linewidth=1,
                label="main analysis",
            )
            ax.legend(fontsize=8)
    fig.suptitle(
        "Sensitivity Analysis: Niche Stability Under Parameter Variation",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


# ============================================================
# 特征基因提取
# ============================================================

def _rank_niche_genes(
    adata: ad.AnnData,
    labels: pd.Series,
    top_n: int,
) -> pd.DataFrame:
    """
    对 niche 高分 spot 进行差异基因分析，筛选特征性高表达基因。

    排序指标：log2FC = log2((mean_high + 1) / (mean_low + 1))，加 1 为伪计数。

    参数
    ----
    adata  : AnnData 数据对象（提供基因表达矩阵）
    labels : 布尔 Series，True = niche_high
    top_n  : 返回前 top_n 个基因

    返回
    ----
    pd.DataFrame，列：gene / mean_high / mean_low / log2_fc，按 log2_fc 降序。
    """
    expr = _expression_frame(adata)
    high = labels.astype(bool)
    if int(high.sum()) < 3 or int((~high).sum()) < 3:
        logging.warning("Too few spots in niche_high/niche_low; skipping signature ranking.")
        return pd.DataFrame(columns=["gene", "mean_high", "mean_low", "log2_fc"])
    mean_high = expr.loc[high].mean(axis=0)
    mean_low  = expr.loc[~high].mean(axis=0)
    ranked = pd.DataFrame({
        "gene":      expr.columns,
        "mean_high": mean_high.to_numpy(),
        "mean_low":  mean_low.to_numpy(),
    })
    ranked["log2_fc"] = np.log2((ranked["mean_high"] + 1.0) / (ranked["mean_low"] + 1.0))
    ranked = ranked.sort_values("log2_fc", ascending=False)
    return ranked.head(top_n)


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

    # ── Step 14: 特征基因提取 ────────────────────────────────────────────────────
    logging.info("Step 14: Extracting niche signature genes...")
    ranked = _rank_niche_genes(adata, df["niche_high"], args.top_niche_genes)
    ranked.to_csv(
        args.out_dir / "immunosuppressive_niche_signature_genes_ranked.csv",
        index=False,
    )
    sig_path = args.out_dir / "immunosuppressive_niche_signature_genes.txt"
    sig_text = "\n".join(ranked["gene"].tolist()) + "\n"
    sig_path.write_text(sig_text, encoding="utf-8")

    args.signature_out.parent.mkdir(parents=True, exist_ok=True)
    args.signature_out.write_text(sig_text, encoding="utf-8")

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
