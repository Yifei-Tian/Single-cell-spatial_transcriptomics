from time import perf_counter
from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq


def check_cell2location_available():
    """导入 cell2location 模型类。"""
    from cell2location.models import Cell2location, RegressionModel

    return RegressionModel, Cell2location


def load_scrna_h5ad(path, t0=None):
    """读取 scRNA 的 h5ad 并进行基础质控过滤。"""
    if t0 is None:
        t0 = perf_counter()

    print(f"-> Loading directly from h5ad: {path} ...")
    adata_sc = ad.read_h5ad(path)

    print(f"<- Loaded AnnData in {perf_counter() - t0:.2f}s")
    print(f"AnnData shape: {adata_sc.shape}")
    # 以上读入并输出必要信息

    obs = cast(pd.DataFrame, adata_sc.obs).copy()
    if "Cell" in obs.columns:
        obs = obs.set_index("Cell") # 将 Cell 列设置为行索引
    adata_sc.obs = obs

    adata_sc.var["SYMBOL"] = adata_sc.var_names # 将 var_names 复制到 SYMBOL 列，方便后续使用
    # gene_id 列不是Ensembl ID命名方式

    gene_filter = np.array(adata_sc.X.sum(axis=0)).flatten() > 10 # 按列求和，只保留总表达量>10的基因
    adata_sc = adata_sc[:, gene_filter]
    cell_filter = np.array(adata_sc.X.sum(axis=1)).flatten() > 200 # 按行求和，只保留总表达量>200的细胞
    adata_sc = adata_sc[cell_filter, :]

    print("After filtering:")
    print(adata_sc.shape) # 输出筛选后的 AnnData 维度

    celltype_counts = (
        adata_sc.obs["celltype"]
        .astype("string")
        .fillna("NA")
        .value_counts()
        .rename_axis("celltype")
        .reset_index(name="count")
    )
    print("celltype 统计结果:")
    print(celltype_counts)

    return adata_sc


def load_visium(path, sample_name="slice"):
    """读取 Visium 空间数据并完成质控、归一化和聚类。"""
    adata_vis = sq.read.visium(path, load_images=True)
    print(adata_vis)

    adata_vis.var_names_make_unique() # 确保基因名唯一，避免后续分析出错
    adata_vis.var["mt"] = adata_vis.var_names.str.startswith("MT-") # 标记线粒体基因
    sc.pp.calculate_qc_metrics(adata_vis, qc_vars=["mt"], inplace=True) # 在 data_vis 新增三列：total_counts、n_genes_by_counts、pct_counts_mt
    print(adata_vis.var["mt"].value_counts()) # 输出线粒体基因数量统计

    adata_vis.obs["slice"] = sample_name
    
    # 筛选空间数据：总表达量在 5000-35000 之间，线粒体基因比例 < 20%，基因至少在 10 个位置表达
    sc.pp.filter_cells(adata_vis, min_counts=5000)
    sc.pp.filter_cells(adata_vis, max_counts=35000)
    adata_vis = adata_vis[adata_vis.obs["pct_counts_mt"] < 20].copy()
    sc.pp.filter_genes(adata_vis, min_cells=10)

    adata_vis.layers["counts"] = adata_vis.X.copy() # 将原始表达矩阵存在一个单独的count矩阵层中

    sc.pp.normalize_total(adata_vis, inplace=True) # 归一化
    sc.pp.log1p(adata_vis) # 对数转换
    sc.pp.highly_variable_genes(adata_vis, flavor="seurat", n_top_genes=2000)

    sc.pp.pca(adata_vis)
    sc.pp.neighbors(adata_vis)
    sc.tl.umap(adata_vis)
    sc.tl.leiden(adata_vis, key_added="clusters")

    return adata_vis


def align_shared_genes(adata_sc, adata_vis):
    """对齐 scRNA 与空间数据的基因，并返回对齐后的拷贝。"""
    shared_genes = adata_sc.var_names.intersection(adata_vis.var_names)
    adata_sc = adata_sc[:, shared_genes].copy()
    adata_vis_sh = adata_vis[:, shared_genes].copy()
    return shared_genes, adata_sc, adata_vis_sh


def subset_t_cells(adata_sc, celltype_col="celltype", cluster_col="res.3"):
    """筛选 T/NK 细胞并将聚类标签转为类别型便于绘图。"""
    adata_t = adata_sc[adata_sc.obs[celltype_col] == "T/NK"].copy()
    adata_t.obs[cluster_col] = adata_t.obs[cluster_col].astype(str).astype("category") # 将 adata_t 的 res.3 列转换为字符串类型再转换为类别型
    return adata_t


def assign_treg_label(
    adata_sc,
    celltype_col="celltype",
    cluster_col="res.3",
    treg_clusters=("9", "9.0"),
    output_col="final_celltype",
):
    """在新注释列中将指定 T/NK 簇重命名为 Treg。"""
    adata_sc.obs[output_col] = adata_sc.obs[celltype_col].astype(str) # 把celltype列复制到final_celltype列
    treg_mask = (adata_sc.obs[celltype_col] == "T/NK") & (
        adata_sc.obs[cluster_col].astype(str).isin(treg_clusters)
    )
    adata_sc.obs.loc[treg_mask, output_col] = "Treg" # 将final_celltype列选择出来的细胞标记为 Treg
    return adata_sc


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
    从 scRNA 数据中估算参考细胞类型签名（RegressionModel）。

    参数
    ----
    adata_sc                  : 已完成注释（含 final_celltype 列）的 scRNA AnnData
    regression_model_cls      : cell2location.models.RegressionModel 类
    max_epochs                : 最大训练轮次上限，默认 250（官方推荐 200–300）
    batch_size                : 每轮训练使用的细胞数，默认 1024
    early_stopping            : 是否启用 Early Stopping，默认启用
    early_stopping_patience   : 连续多少轮损失无改善则提前停止，默认 30 轮
    early_stopping_min_delta  : 损失改善的最小阈值，低于此值视为"无改善"

    说明
    ----
    RegressionModel 学习每种细胞类型（final_celltype）在各基因上的平均表达特征，
    作为后续 Cell2location 空间建模的参考签名矩阵（cell_state_df）。
    训练轮次过少（< 200）容易导致签名矩阵未充分收敛，影响反卷积精度；
    Early Stopping 可在损失已经平稳时自动终止训练，避免无效等待。

    关于学习率
    ----------
    scvi-tools / cell2location 的 RegressionModel 默认学习率为 5e-4，
    该值已经过官方调优，不应随意修改；学习率过大（如 > 1e-3）会导致
    -ELBO 在训练后期剧烈震荡而无法收敛（如本项目实测图所示）。
    """
    regression_model_cls.setup_anndata(
        adata=adata_sc,
        labels_key="final_celltype",  # 指定细胞类型注释列
    )

    model = regression_model_cls(adata_sc)

    # cell2location / scvi-tools 的 train() 接口说明：
    #   - max_epochs, batch_size：直接传给 train()
    #   - early_stopping*：scvi-tools >= 0.20 支持直接作为 train() 关键字参数
    #   - 学习率：通过 plan_kwargs={"lr": ...} 传递，默认 5e-4 无需修改
    #   - 不额外传 plan_kwargs，保留框架默认学习率（5e-4），避免震荡
    train_kwargs: dict = dict(
        max_epochs=max_epochs,
        batch_size=batch_size,
    )
    if early_stopping:
        train_kwargs["early_stopping"] = True
        train_kwargs["early_stopping_patience"] = early_stopping_patience
        train_kwargs["early_stopping_min_delta"] = early_stopping_min_delta
        # 监控训练集损失（RegressionModel 无验证集，只有 train_loss_epoch）
        train_kwargs["early_stopping_monitor"] = "train_loss_epoch"

    model.train(**train_kwargs)

    # 读取实际训练轮次：scvi-tools 将损失记录在 model.history 中，
    # 键名为 "train_loss_epoch"（注意不是 "elbo_train"）
    actual_epochs: int | str = "unknown"
    if hasattr(model, "history"):
        history = model.history
        # history 可能是 dict-like 或 scvi History 对象，统一用 [] 或 .get()
        for key in ("train_loss_epoch", "elbo_train", "training_loss"):
            try:
                records = history[key]
                # records 可能是 list、ndarray 或 pd.Series
                actual_epochs = int(len(records))
                break
            except (KeyError, TypeError):
                continue

    print(f"RegressionModel 训练完成（实际训练轮次: {actual_epochs} / 上限 {max_epochs}）")
    return model  # 训练得到了每个细胞类型的参考表达谱


def export_signatures(model, adata_sc):
    """导出 RegressionModel 后验中的细胞类型参考签名。"""
    exported = model.export_posterior(adata_sc)
    return adata_sc if exported is None else exported


def extract_cell_state_df(adata_sc):
    """从参考模型后验中提取细胞类型 signature 矩阵。"""
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

    cell_state_df.index = adata_sc.var_names
    cell_state_df.columns = [
        str(col).replace("means_per_cluster_mu_fg_", "") for col in cell_state_df.columns
    ]
    return cell_state_df # 内容是每个基因在每类细胞中的参考表达签名


def sanitize_cell_state_df(cell_state_df, eps=1e-5):
    """清理 signature 矩阵中的无效值。"""
    df = cell_state_df.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    df = df.clip(lower=0.0)
    if float(df.to_numpy().sum()) <= 0:
        raise ValueError(
            "cell_state_df 全部为 0 或无效值，请检查 RegressionModel.export_posterior() "
            "以及 means_per_cluster_mu_fg 的列名是否正确。"
        )
    zero_rows = df.sum(axis=1) <= 0
    if bool(zero_rows.any()):
        df.loc[zero_rows, :] = eps
    return df


def _strip_abundance_prefix(col: str) -> str:
    """
    去除 cell2location 后验矩阵列名中的统计量前缀，只保留细胞类型名称。

    cell2location 导出的 obsm 列名格式为 "<统计量前缀>_<细胞类型>"，例如：
      means_cell_abundance_w_sf_Hepatocyte  →  Hepatocyte
      q05_cell_abundance_w_sf_T/NK          →  T/NK
      q95_cell_abundance_w_sf_Fibroblast    →  Fibroblast
      means_per_cluster_mu_fg_Treg          →  Treg

    采用正则匹配所有已知前缀模式，保证无论使用哪种后验统计量导出，
    最终列名都是干净的细胞类型名称。
    """
    import re
    # 匹配所有已知的 cell2location / scvi 后验统计量前缀
    # 格式：(means|q05|q95|median|q25|q75)[_]cell_abundance_w_sf_ 或 means_per_cluster_mu_fg_
    cleaned = re.sub(
        r"^(means|q\d+|median)_?cell_abundance_w_sf_",
        "",
        col,
    )
    cleaned = re.sub(r"^means_per_cluster_mu_fg_", "", cleaned)
    return cleaned


def compute_spot_cell_proportion(
    adata_sp,
    abundance_key="means_cell_abundance_w_sf",
    spot_id_col="spot_id",
):
    """
    将 Cell2location 后验细胞丰度转换为每个 spot 的细胞类型比例表。

    列名处理：自动去除 cell2location 添加的统计量前缀（如 means_cell_abundance_w_sf_），
    输出列名直接为细胞类型名称（如 Hepatocyte、Treg、T/NK 等）。

    返回 DataFrame：行=spot，列=细胞类型（首列为 spot_id）。
    """
    if "cell_abundance" in adata_sp.obsm:
        abundance = adata_sp.obsm["cell_abundance"]
    else:
        abundance = adata_sp.obsm[abundance_key]

    if not isinstance(abundance, pd.DataFrame):
        cell_types = adata_sp.uns.get("mod", {}).get("factor_names")
        abundance = pd.DataFrame(abundance, index=adata_sp.obs_names, columns=cell_types)
    else:
        abundance = abundance.copy()

    # 清理列名：去掉所有统计量前缀，只保留细胞类型名称
    abundance.columns = [_strip_abundance_prefix(str(col)) for col in abundance.columns]

    denom = abundance.sum(axis=1).replace(0, np.nan)
    proportions = abundance.div(denom, axis=0).fillna(0.0)
    proportions.insert(0, spot_id_col, proportions.index)
    return proportions
# 每个 spot 的细胞类型比例表（行=spot，列=cell_type，首列为 spot_id）
