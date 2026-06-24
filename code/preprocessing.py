"""
================================================================================
模块名称: preprocessing.py
功能概述: 数据预处理与 Cell2location 反卷积核心函数库
================================================================================

【模块说明】
    本模块封装了单细胞转录组（scRNA-seq）和空间转录组（Visium）数据预处理的
    核心函数，以及 Cell2location 两阶段建模的完整流程。
    所有函数均通过 run_preprocessing.py 主脚本调用。

【函数索引】
    数据加载：
      check_cell2location_available()     - 导入并返回 Cell2location 模型类
      load_scrna_h5ad()                   - 加载 scRNA-seq .h5ad，基础质控
      load_visium()                       - 加载 Visium 空间数据，完整质控流程

    预处理：
      align_shared_genes()               - 对齐 scRNA 与 Visium 的共享基因
      subset_t_cells()                   - 提取 T/NK 细胞亚群
      assign_treg_label()                - 将指定 T/NK 簇重注释为 Treg

    Cell2location 建模：
      setup_and_train_regression_model() - 训练 RegressionModel（参考签名学习）
      export_signatures()                - 导出 RegressionModel 后验签名
      extract_cell_state_df()            - 提取细胞类型参考签名矩阵
      sanitize_cell_state_df()           - 清理签名矩阵中的无效值

    后处理：
      _strip_abundance_prefix()          - 去除 cell2location 列名前缀
      compute_spot_cell_proportion()     - 将后验丰度转换为细胞比例表

【输入文件】（由调用方 run_preprocessing.py 传入，本模块不直接读取文件）
    data/scRNA_reference.h5ad          - scRNA-seq 参考数据（load_scrna_h5ad 读取）
    data/chc20_visium.h5ad             - CHC20 Visium 空间数据（load_visium 读取）
    data/chc23_visium.h5ad             - CHC23 Visium 空间数据（load_visium 读取）
    或 data/CHC20_Visium/              - Space Ranger 输出目录（load_visium 直接读取原始格式）
    或 data/CHC23_Visium/              - Space Ranger 输出目录

【输出文件】（由调用方 run_preprocessing.py 接收并保存，本模块以返回值形式输出）
    （内存中）inf_aver DataFrame       - Cell2location 参考签名矩阵（extract_cell_state_df 输出）
    （内存中）adata_vis AnnData        - 含细胞丰度估计的空间数据（setup_and_train_regression_model 输出）
    （内存中）proportion DataFrame     - 每个 spot 的细胞类型比例表（compute_spot_cell_proportion 输出）
================================================================================
"""
from time import perf_counter
from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


# ============================================================
# Cell2location 模型导入
# ============================================================

def check_cell2location_available():
    """
    导入并返回 Cell2location 两阶段模型类。

    返回
    ----
    (RegressionModel, Cell2location) — 参考签名学习模型类和空间建模模型类。

    注意
    ----
    如果 cell2location 包未安装，将抛出 ImportError。
    建议在环境配置阶段运行 `pip install cell2location` 确保安装。
    """
    from cell2location.models import Cell2location, RegressionModel
    return RegressionModel, Cell2location


# ============================================================
# 数据加载函数
# ============================================================

def load_scrna_h5ad(path, t0=None):
    """
    加载 scRNA-seq .h5ad 文件并执行基础质控过滤。

    质控标准（参考 Luecken & Theis, Mol Syst Biol, 2019）：
      - 保留总表达量 > 200 的细胞（过低可能为空液滴或碎片）；
      - 保留总表达量 > 10 的基因（过低为噪声信号）。

    参数
    ----
    path : str 或 Path，.h5ad 文件路径
    t0   : 计时起始时间戳（perf_counter），用于日志记录

    返回
    ----
    ad.AnnData，质控后的 scRNA-seq 数据对象。
    """
    if t0 is None:
        t0 = perf_counter()

    print(f"-> Loading scRNA-seq from h5ad: {path} ...")
    adata_sc = ad.read_h5ad(path)
    print(f"<- Loaded in {perf_counter() - t0:.2f}s | Shape: {adata_sc.shape}")

    # 处理行索引：若 obs 中有 "Cell" 列，将其设为行索引
    obs = cast(pd.DataFrame, adata_sc.obs).copy()
    if "Cell" in obs.columns:
        obs = obs.set_index("Cell")
    adata_sc.obs = obs

    # 将 var_names 备份至 "SYMBOL" 列，便于后续基因名引用
    adata_sc.var["SYMBOL"] = adata_sc.var_names

    # 质控过滤
    gene_filter = np.array(adata_sc.X.sum(axis=0)).flatten() > 10
    adata_sc = adata_sc[:, gene_filter]
    cell_filter = np.array(adata_sc.X.sum(axis=1)).flatten() > 200
    adata_sc = adata_sc[cell_filter, :]
    print(f"After QC filtering: {adata_sc.shape}")

    # 输出细胞类型统计
    celltype_counts = (
        adata_sc.obs["celltype"]
        .astype("string")
        .fillna("NA")
        .value_counts()
        .rename_axis("celltype")
        .reset_index(name="count")
    )
    print("Cell type counts:\n", celltype_counts.to_string(index=False))

    return adata_sc


def load_visium(path, sample_name: str = "slice"):
    """
    加载 10x Genomics Visium 空间转录组数据，执行完整质控、归一化和聚类。

    质控标准（参考 Williams et al., Genome Medicine, 2022）：
      - 保留总 UMI 计数在 5,000–35,000 之间的 spot；
      - 过滤线粒体基因比例超过 20% 的 spot（提示细胞损伤）；
      - 仅保留在至少 10 个 spot 中表达的基因。

    参数
    ----
    path        : str 或 Path，Space Ranger 输出目录路径
    sample_name : 切片名称标签（写入 adata.obs["slice"]）

    返回
    ----
    ad.AnnData，质控并完成归一化/聚类的 Visium 数据对象。
    """
    adata_vis = sq.read.visium(path, load_images=True)
    print(f"Loaded Visium [{sample_name}]: {adata_vis.shape}")

    adata_vis.var_names_make_unique()
    adata_vis.var["mt"] = adata_vis.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata_vis, qc_vars=["mt"], inplace=True)
    print("MT gene count:", adata_vis.var["mt"].sum())

    adata_vis.obs["slice"] = sample_name

    # 质控过滤
    sc.pp.filter_cells(adata_vis, min_counts=5000)
    sc.pp.filter_cells(adata_vis, max_counts=35000)
    adata_vis = adata_vis[adata_vis.obs["pct_counts_mt"] < 20].copy()
    sc.pp.filter_genes(adata_vis, min_cells=10)
    print(f"After QC filtering [{sample_name}]: {adata_vis.shape}")

    # 保存原始 count 层（Cell2location 需要原始 count 作为输入）
    adata_vis.layers["counts"] = adata_vis.X.copy()

    # 归一化与降维（用于质控可视化和探索性分析）
    sc.pp.normalize_total(adata_vis, inplace=True)
    sc.pp.log1p(adata_vis)
    sc.pp.highly_variable_genes(adata_vis, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata_vis)
    sc.pp.neighbors(adata_vis)
    sc.tl.umap(adata_vis)
    sc.tl.leiden(adata_vis, key_added="clusters")

    return adata_vis


# ============================================================
# 预处理函数
# ============================================================

def align_shared_genes(adata_sc, adata_vis):
    """
    对齐 scRNA-seq 与 Visium 数据集的共享基因，返回对齐后的副本。

    Cell2location 要求 scRNA 和 Visium 数据使用相同的基因集。
    取两个数据集 var_names 的交集，确保基因维度一致。

    参数
    ----
    adata_sc  : scRNA-seq AnnData
    adata_vis : Visium AnnData

    返回
    ----
    (shared_genes, adata_sc_aligned, adata_vis_aligned) — 三元组：
      shared_genes      : 共享基因列表（pd.Index）；
      adata_sc_aligned  : 只含共享基因的 scRNA AnnData；
      adata_vis_aligned : 只含共享基因的 Visium AnnData。
    """
    shared_genes = adata_sc.var_names.intersection(adata_vis.var_names)
    print(f"Shared genes: {len(shared_genes)}")
    adata_sc_aligned  = adata_sc[:, shared_genes].copy()
    adata_vis_aligned = adata_vis[:, shared_genes].copy()
    return shared_genes, adata_sc_aligned, adata_vis_aligned


def subset_t_cells(adata_sc, celltype_col: str = "celltype", cluster_col: str = "res.3"):
    """
    从 scRNA-seq 数据中提取 T/NK 细胞亚群，并将聚类标签转为 category 类型。

    参数
    ----
    adata_sc    : scRNA-seq AnnData
    celltype_col: 细胞类型注释列名
    cluster_col : 聚类结果列名（用于 Treg 亚群鉴定）

    返回
    ----
    ad.AnnData，只含 T/NK 细胞的子集。
    """
    adata_t = adata_sc[adata_sc.obs[celltype_col] == "T/NK"].copy()
    # 将聚类标签转为 category 类型，确保 scanpy 绘图函数正确识别分类变量
    adata_t.obs[cluster_col] = (
        adata_t.obs[cluster_col].astype(str).astype("category")
    )
    return adata_t


def assign_treg_label(
    adata_sc,
    celltype_col: str = "celltype",
    cluster_col: str = "res.3",
    treg_clusters: tuple[str, ...] = ("9", "9.0"),
    output_col: str = "final_celltype",
):
    """
    在新注释列中将指定 T/NK 亚聚类重命名为 Treg。

    鉴定依据（参考 Sakaguchi et al., Annual Review of Immunology, 2020）：
      该方法遵循基于标志基因的细胞类型注释经典范式：
      CD3D + CD4 为 T 细胞谱系标志；
      FOXP3 为 Treg 主转录因子；
      IL2RA（CD25）为 Treg 表面标志物。
      通过气泡图确认 res.3 第 9 簇具备上述表达特征后，将其重注释为 Treg。

    参数
    ----
    adata_sc     : scRNA-seq AnnData
    celltype_col : 原始细胞类型注释列名
    cluster_col  : 聚类列名
    treg_clusters: 对应 Treg 的聚类编号（兼容字符串和浮点格式）
    output_col   : 新注释列名

    返回
    ----
    ad.AnnData，新增 output_col 列，Treg 细胞已重新标注。
    """
    adata_sc.obs[output_col] = adata_sc.obs[celltype_col].astype(str)
    treg_mask = (adata_sc.obs[celltype_col] == "T/NK") & (
        adata_sc.obs[cluster_col].astype(str).isin(treg_clusters)
    )
    adata_sc.obs.loc[treg_mask, output_col] = "Treg"
    n_treg = int(treg_mask.sum())
    print(f"Treg cells annotated: {n_treg}")
    return adata_sc


# ============================================================
# Cell2location 建模函数
# ============================================================

def setup_and_train_regression_model(
    adata_sc,
    regression_model_cls,
    max_epochs: int = 250,
    batch_size: int = 1024,
    early_stopping: bool = True,
    early_stopping_patience: int = 30,
    early_stopping_min_delta: float = 1e-4,
):
    """
    训练 Cell2location RegressionModel，从 scRNA-seq 数据中学习细胞类型参考签名。

    模型说明（参考 Kleshchevnikov et al., Nature Biotechnology, 2022）
    ---------------------------------------------------------------
    RegressionModel 是 Cell2location 的第一阶段模型，通过变分推断学习
    每种细胞类型（final_celltype）在各基因上的平均表达分布参数，
    输出用于第二阶段空间建模的参考签名矩阵（cell_state_df）。

    训练参数选择
    -----------
    - max_epochs=250    : 官方推荐 200–300 epochs，确保模型充分收敛；
    - batch_size=1024   : 适合常规 GPU/CPU 内存大小；
    - Early Stopping    : 损失平稳时自动终止，避免无效等待和过拟合；
    - 学习率            : 使用框架默认值 5e-4（经官方调优，不应随意修改）。

    参数
    ----
    adata_sc                 : 已完成 Treg 注释的 scRNA AnnData
    regression_model_cls     : RegressionModel 类
    max_epochs               : 最大训练轮次上限
    batch_size               : 每轮训练批大小
    early_stopping           : 是否启用早停
    early_stopping_patience  : 早停等待轮次
    early_stopping_min_delta : 早停最小改善阈值

    返回
    ----
    训练完成的 RegressionModel 对象。
    """
    regression_model_cls.setup_anndata(
        adata=adata_sc,
        labels_key="final_celltype",
    )
    model = regression_model_cls(adata_sc)

    train_kwargs: dict = dict(max_epochs=max_epochs, batch_size=batch_size)
    if early_stopping:
        train_kwargs.update({
            "early_stopping":         True,
            "early_stopping_patience":  early_stopping_patience,
            "early_stopping_min_delta": early_stopping_min_delta,
            "early_stopping_monitor":   "elbo_train",
        })

    model.train(**train_kwargs)

    # 读取实际训练轮次（兼容不同版本的 scvi-tools history 格式）
    actual_epochs: int | str = "unknown"
    if hasattr(model, "history"):
        for key in ("train_loss_epoch", "elbo_train", "training_loss"):
            try:
                actual_epochs = int(len(model.history[key]))
                break
            except (KeyError, TypeError):
                continue

    print(f"RegressionModel 训练完成（实际轮次: {actual_epochs} / 上限 {max_epochs}）")
    return model


def export_signatures(model, adata_sc):
    """
    导出 RegressionModel 后验签名至 adata_sc.varm。

    参数
    ----
    model    : 训练完成的 RegressionModel 对象
    adata_sc : 训练时使用的 scRNA AnnData

    返回
    ----
    ad.AnnData，后验签名已写入 varm["means_per_cluster_mu_fg"] 的数据对象。
    """
    exported = model.export_posterior(adata_sc)
    return adata_sc if exported is None else exported


def extract_cell_state_df(adata_sc):
    """
    从 RegressionModel 后验中提取细胞类型参考签名矩阵（cell_state_df）。

    cell_state_df 的结构：
      - 行：共享基因（gene × cell_type）
      - 列：细胞类型名称（已去除统计量前缀）
      - 值：每个基因在每种细胞类型中的平均表达强度

    参数
    ----
    adata_sc : 已通过 export_signatures 导出后验的 scRNA AnnData

    返回
    ----
    pd.DataFrame，行=基因，列=细胞类型，值=参考表达签名。
    """
    signatures = adata_sc.varm["means_per_cluster_mu_fg"]
    if isinstance(signatures, pd.DataFrame):
        cell_state_df = signatures.copy()
    else:
        factor_names = adata_sc.uns["mod"]["factor_names"]
        cell_state_df = pd.DataFrame(
            signatures,
            index=adata_sc.var_names,
            columns=factor_names,
        )
    cell_state_df.index   = adata_sc.var_names
    cell_state_df.columns = [
        str(col).replace("means_per_cluster_mu_fg_", "")
        for col in cell_state_df.columns
    ]
    return cell_state_df


def sanitize_cell_state_df(cell_state_df, eps: float = 1e-5):
    """
    清理 cell_state_df 签名矩阵中的无效值（Inf、NaN、负值、全零行）。

    确保矩阵满足 Cell2location 空间建模的数值稳定性要求：
      - 将 Inf 和 NaN 替换为 0；
      - 将负值截断为 0；
      - 对全零行用小正数 eps 填充（避免概率计算中出现零）。

    参数
    ----
    cell_state_df : 原始签名矩阵
    eps           : 全零行的填充值（默认 1e-5）

    返回
    ----
    pd.DataFrame，清理后的签名矩阵。
    """
    df = cell_state_df.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    df = df.clip(lower=0.0)
    if float(df.to_numpy().sum()) <= 0:
        raise ValueError(
            "cell_state_df 全为 0 或无效值，请检查 RegressionModel.export_posterior() "
            "以及 means_per_cluster_mu_fg 的列名是否正确。"
        )
    zero_rows = df.sum(axis=1) <= 0
    if bool(zero_rows.any()):
        df.loc[zero_rows, :] = eps
    return df


# ============================================================
# 后处理函数
# ============================================================

def _strip_abundance_prefix(col: str) -> str:
    """
    去除 Cell2location 后验矩阵列名中的统计量前缀，只保留细胞类型名称。

    处理以下前缀格式：
      means_cell_abundance_w_sf_Hepatocyte → Hepatocyte
      q05_cell_abundance_w_sf_T/NK         → T/NK
      q95_cell_abundance_w_sf_Fibroblast   → Fibroblast
      means_per_cluster_mu_fg_Treg         → Treg

    参数
    ----
    col : 原始列名字符串

    返回
    ----
    str，清洗后的细胞类型名称。
    """
    import re
    cleaned = re.sub(r"^(means|q\d+|median)_?cell_abundance_w_sf_", "", col)
    cleaned = re.sub(r"^means_per_cluster_mu_fg_", "", cleaned)
    return cleaned


def compute_spot_cell_proportion(
    adata_sp,
    abundance_key: str = "means_cell_abundance_w_sf",
    spot_id_col: str = "spot_id",
):
    """
    将 Cell2location 后验细胞丰度转换为每个 spot 的细胞类型比例表。

    处理流程：
      1. 从 adata.obsm 中提取细胞丰度矩阵；
      2. 自动去除列名中的统计量前缀，保留干净的细胞类型名称；
      3. 每个 spot 各细胞类型丰度除以该 spot 总丰度，得到 0–1 的比例值；
      4. 在首列插入 spot_id。

    参数
    ----
    adata_sp      : Cell2location 建模后的空间 AnnData
    abundance_key : 后验丰度矩阵在 adata.obsm 中的键名
    spot_id_col   : 输出 DataFrame 中 spot ID 列的列名

    返回
    ----
    pd.DataFrame，行=spot，列=细胞类型（首列为 spot_id），值=归一化比例。
    """
    if "cell_abundance" in adata_sp.obsm:
        abundance = adata_sp.obsm["cell_abundance"]
    else:
        abundance = adata_sp.obsm[abundance_key]

    if not isinstance(abundance, pd.DataFrame):
        cell_types = adata_sp.uns.get("mod", {}).get("factor_names")
        abundance = pd.DataFrame(
            abundance, index=adata_sp.obs_names, columns=cell_types
        )
    else:
        abundance = abundance.copy()

    # 清理列名：去掉统计量前缀
    abundance.columns = [_strip_abundance_prefix(str(c)) for c in abundance.columns]

    # 归一化为比例
    denom = abundance.sum(axis=1).replace(0, np.nan)
    proportions = abundance.div(denom, axis=0).fillna(0.0)
    proportions.insert(0, spot_id_col, proportions.index)
    return proportions
