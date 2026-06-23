from pathlib import Path

import anndata as ad
import cell2location_pre
import pandas as pd
import scanpy as sc


# 1) 指定项目根目录（保持你服务器上的原始目录）
ROOT_DIR = Path("/home/gyw/R/Python/Single_cell_train")
DATA_DIR = ROOT_DIR / "data"

# 2) 输入/输出路径
PATH_COUNTS = DATA_DIR / "GSE149614_HCC.scRNAseq.S71915.count.txt"
PATH_META = DATA_DIR / "GSE149614_HCC.metadata.updated.txt"
PATH_SCRNA_H5AD = DATA_DIR / "scRNA_reference.h5ad"
PATH_CHC20 = DATA_DIR / "CHC20_Visium"
PATH_CHC23 = DATA_DIR / "CHC23_Visium"


def _normalize_gene_ids(index_obj) -> pd.Index:
    """标准化 gene_id 字符串，去掉 Ensembl 版本号后缀（如 ENSG... .15）。"""
    vals = pd.Series(index_obj.astype(str)).str.strip()
    # index_obj.astype(str) 全转字符串；.str 向量化方法；strip() 去掉首尾空格。
    vals = vals.str.replace(r"^(ENSG\d+)\.\d+$", r"\1", regex=True)
    vals = vals.str.replace(r"^(ENSMUSG\d+)\.\d+$", r"\1", regex=True)
    vals = vals.str.replace(r"^(ENST\d+)\.\d+$", r"\1", regex=True)
    return pd.Index(vals)


def _force_gene_id_index(adata: sc.AnnData, name: str, gene_id_column: str | None = None) -> sc.AnnData:
    """
    改变基因索引
    强制使用 gene_id 作为 var_names；若不满足则直接报错。
    gene_id_column: 提供的用于替换的列名
    """
    if gene_id_column is not None:
        if gene_id_column not in adata.var.columns: # 先检查该列是否存在于基因注释表 adata.var。
            raise KeyError(f"{name} 缺少 `{gene_id_column}` 列，无法切换到 gene_id 索引。")
        if "SYMBOL" not in adata.var.columns:
            adata.var.set_index(gene_id_column, drop=True, inplace=True)
            adata.var["SYMBOL"] = adata.var_names.astype(str) # 把当前 var_names 备份到 SYMBOL，用于可读展示和回溯。

        adata.var.set_index(gene_id_column, drop=True, inplace=True) # 把基因索引改成你指定的 gene_id 列。
        # drop=True 表示把原来的索引列丢弃，不保留为普通列；inplace=True 表示直接修改原 DataFrame，不返回新对象。

    adata.var_names = _normalize_gene_ids(adata.var_names)
    adata.var_names_make_unique() # 保证索引唯一，避免后续合并/建模报重复基因名错误。

    if _guess_id_style(adata.var_names) != "ensembl_like":
        raise ValueError(
            f"{name} 当前 var_names 不是 gene_id/Ensembl 风格，"
            "请先提供 gene_id 列或将输入矩阵第一列改为 gene_id。"
        )
    return adata


def _align_shared_gene_ids(adata_sc: sc.AnnData, adata_vis: sc.AnnData) -> tuple[sc.AnnData, sc.AnnData]:
    """
    寻找单细胞与空间组的交集
    仅保留两者共有的 gene_id，并统一顺序。
    adata_sc 单细胞参考
    adata_vis 空间转录组数据
    返回处理后的 adata_sc 和 adata_vis，保证两者 var_names 完全一致且都是 gene_id 风格。
    """
    shared_gene_ids = adata_sc.var_names.intersection(adata_vis.var_names) # 求两个索引的交集，得到共有的 gene_id 列表
    if len(shared_gene_ids) == 0:
        raise ValueError("scRNA 与 spatial 在 gene_id 层面没有交集，无法继续建模。")

    adata_sc = adata_sc[:, shared_gene_ids].copy() # 新建一个单细胞对象，只保留交集基因，保留所有细胞（行不变）
    adata_vis = adata_vis[:, shared_gene_ids].copy() # 在空间对象中也做相同的事情
    print(f"已统一到 gene_id，交集基因数: {len(shared_gene_ids)}")
    return adata_sc, adata_vis


def _guess_id_style(index_obj) -> str:
    """
        粗略判断索引更像 SYMBOL 还是 Ensembl/gene_ids。
        index_obj: pandas Index 对象，通常是 adata.obs_names 或 adata.var_names；
        .obs_names: AnnData的行索引，通常是细胞条形码；.var_names: AnnData的列索引，通常是基因ID。
    """
    if len(index_obj) == 0:
        return "empty"
    sample_vals = pd.Series(index_obj.astype(str)).head(200)
    # .astype(str) 将索引元素统一转成字符串；pd.Series 转成 Series ，便于使用字符串向量化方法；head(200) 取前200个元素进行判断
    ensembl_ratio = sample_vals.str.startswith(("ENSG", "ENSMUSG", "ENST")).mean()
    # 判断是否以这些前缀开头，得到布尔序列，取平均得到比例
    if ensembl_ratio > 0.6:
        return "ensembl_like"
    return "symbol_or_other"


def _print_adata_snapshot(adata: sc.AnnData, name: str) -> None:
    '''检查细胞和基因名字'''
    print(f"\n=== {name} 概览 ===")
    print(f"shape (行obs x 列var): {adata.shape}")
    print(f"obs_names(前5): {list(adata.obs_names[:5])}")
    print(f"var_names(前5): {list(adata.var_names[:5])}")
    print(f"obs_names是否唯一: {adata.obs_names.is_unique}")
    print(f"var_names是否唯一: {adata.var_names.is_unique}")
    print(f"var_names风格判断: {_guess_id_style(adata.var_names)}")
    print(f"obs列(前10): {list(adata.obs.columns[:10])}")
    print(f"var列(前10): {list(adata.var.columns[:10])}")

    if adata.n_obs > 0:
        print("obs内容示例(前3行):")
        print(adata.obs.head(3))
    if adata.n_vars > 0:
        print("var内容示例(前3行):")
        print(adata.var.head(3))


def _print_overlap(idx_a, idx_b, name_a: str, name_b: str) -> pd.Index:
    '''
    比较两个索引集合，共同和差异的基因
    idx_a, idx_b 通常是 adata_sc.var_names 和 adata_vis.var_names 这样的 pandas Index 对象
    name_a, name_b 是它们的描述字符串，用于打印输出
    函数返回一个 Index 对象，表示 idx_a 和 idx_b 的交集
    '''
    overlap = idx_a.intersection(idx_b)
    only_a = idx_a.difference(idx_b)
    only_b = idx_b.difference(idx_a)

    print(f"\n[{name_a}] vs [{name_b}] 交集统计")
    print(f"{name_a} 总数: {len(idx_a)}")
    print(f"{name_b} 总数: {len(idx_b)}")
    print(f"交集数量: {len(overlap)}")
    print(f"仅{ name_a }有: {len(only_a)}")
    print(f"仅{ name_b }有: {len(only_b)}")
    print(f"交集示例(前5): {list(overlap[:5])}")
    print(f"仅{name_a}示例(前5): {list(only_a[:5])}")
    print(f"仅{name_b}示例(前5): {list(only_b[:5])}")
    return overlap


def diagnose_alignment(adata_sc: sc.AnnData, adata_vis: sc.AnnData) -> None:
    """打印索引对齐诊断信息，帮助判断是否需要切换 index。"""
    print("\n================ 索引对齐诊断开始 ================")
    _print_adata_snapshot(adata_sc, "scRNA") # 单细胞
    _print_adata_snapshot(adata_vis, "Spatial(merged)") # 空间组

    overlap_genes = _print_overlap( # 检查单细胞与空间组中重叠与差异基因
        adata_sc.var_names, adata_vis.var_names, "scRNA.var_names", "Spatial.var_names"
    )
    _print_overlap( # 检查单细胞与空间组中重叠与差异细胞/spot
        adata_sc.obs_names, adata_vis.obs_names, "scRNA.obs_names", "Spatial.obs_names"
    )

    if "SYMBOL" in adata_vis.var.columns:
        overlap_symbol = adata_sc.var_names.intersection(adata_vis.var["SYMBOL"].astype(str))
        print("\nscRNA.var_names 与 Spatial.var['SYMBOL'] 交集数量:", len(overlap_symbol))
        print("交集示例(前5):", list(overlap_symbol[:5]))

    if len(overlap_genes) == 0:
        print("\n[提示] 基因索引交集为 0：通常说明两边命名体系不同（如 SYMBOL vs ENSG）。")
        print("[提示] 你需要把 scRNA 与 spatial 的 var_names 统一到同一体系后再建模。")
    else:
        print("\n[提示] 已存在基因交集，可继续评估交集规模是否足够。")

    if "sample" in adata_vis.obs.columns:
        print("Spatial每个切片spot数量:")
        print(adata_vis.obs["sample"].value_counts())
    print("================ 索引对齐诊断结束 ================\n")


def _require_exists(path: Path, desc: str) -> None:
    '''
    前置校验函数，检查某个路径是否存在
    path: 要检查的文件/目录路径
    desc: 路径的描述信息，用于错误提示
    None: 如果路径存在，不返回；如果路径不存在，抛出 FileNotFoundError 异常
    '''
    if not path.exists():
        raise FileNotFoundError(f"{desc} 不存在: {path}")


def build_scrna_reference() -> sc.AnnData: # ->sc.AnnData 返回值注释
    """将 txt 计数矩阵 + metadata 转成 h5ad，并做索引对齐。"""
    _require_exists(PATH_COUNTS, "counts 文件")
    _require_exists(PATH_META, "metadata 文件")

    print("正在读取表达矩阵，7万细胞可能需要1-2分钟...")
    counts = pd.read_csv(PATH_COUNTS, sep="\t", index_col=0)
    adata_sc = ad.AnnData(counts.T) # 将读入的count转化为anndata对象

    print("正在读取并对齐 Metadata...")
    meta = pd.read_csv(PATH_META, sep="\t", index_col=0)
    common_cells = adata_sc.obs_names.intersection(meta.index) # obs_names AnnData的观测名，本质是一个pandas Index
    # intersection 求两个索引的交集
    if len(common_cells) == 0:
        raise ValueError("counts 与 metadata 没有共同细胞条形码，请检查索引。")

    adata_sc = adata_sc[common_cells].copy() # 只保留共同细胞
    adata_sc.obs = meta.loc[common_cells].copy() # 对齐 Metadata 到 adata_sc.obs
    adata_sc = _force_gene_id_index(adata_sc, "scRNA")
    adata_sc.var_names_make_unique()
    adata_sc.write(PATH_SCRNA_H5AD)
    print(f"scRNA 参考数据已保存: {PATH_SCRNA_H5AD}")
    return adata_sc


def _load_visium_slice(path_visium: Path, sample_name: str) -> sc.AnnData:
    _require_exists(path_visium, f"{sample_name} Visium 目录")

    adata_slice = sc.read_visium(path_visium)
    adata_slice.var_names_make_unique()
    adata_slice.obs["sample"] = sample_name
    adata_slice.var["SYMBOL"] = adata_slice.var_names

    adata_slice = _force_gene_id_index(adata_slice, f"{sample_name} spatial", gene_id_column="gene_ids")

    return adata_slice


def load_and_merge_visium() -> sc.AnnData:
    """读取 CHC20/CHC23 并合并为 adata_vis。"""
    adata_chc20 = _load_visium_slice(PATH_CHC20, "CHC20")
    adata_chc23 = _load_visium_slice(PATH_CHC23, "CHC23")

    adata_vis = ad.concat(
        [adata_chc20, adata_chc23],
        join="inner",
        label="sample",
        keys=["CHC20", "CHC23"],
        index_unique="-",
        merge="same",
    )
    adata_vis.layers["counts"] = adata_vis.X.copy()

    print("合并后的空间数据:")
    print(adata_vis)
    return adata_vis


if __name__ == "__main__":
    adata_sc = build_scrna_reference()
    adata_vis = load_and_merge_visium()
    adata_sc, adata_vis = _align_shared_gene_ids(adata_sc, adata_vis)
    diagnose_alignment(adata_sc, adata_vis)

    if "cell_type" not in adata_sc.obs.columns:
        raise KeyError("metadata 中缺少 `cell_type` 列，无法训练 RegressionModel。")

    print("开始训练 reference model...")
    cell2location.models.RegressionModel.setup_anndata(adata_sc, labels_key="cell_type")
    model = cell2location.models.RegressionModel(adata_sc)
    model.train(max_epochs=250)

    print("reference model 训练完成。")
    print(f"可用于后续 Cell2location 空间建模的对象: adata_vis (n_obs={adata_vis.n_obs})")

