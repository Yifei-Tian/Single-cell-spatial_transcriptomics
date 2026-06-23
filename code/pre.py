"""
pre.py — 原始数据整理脚本
============================
职责：将各种原始格式的数据统一转换为 AnnData h5ad 格式，供后续分析流程直接读取。
本脚本 **不做** 任何模型训练、归一化或下游分析，只做数据格式转换与结构检查。

支持的转换任务：
  1. scRNA-seq：GSE149614 txt count 矩阵 + metadata → scRNA_reference.h5ad
  2. Visium 空间转录组：CHC20 / CHC23 Space Ranger 输出目录 → chc20_visium.h5ad / chc23_visium.h5ad
  3. 可选：将两张 Visium 切片合并保存为 merged_visium.h5ad

运行方式：
  python pre.py                          # 转换所有数据
  python pre.py --skip-scrna             # 跳过 scRNA 转换（已有 h5ad 时）
  python pre.py --skip-visium            # 跳过 Visium 转换
  python pre.py --no-merge               # 不生成合并的 Visium h5ad
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


# ---------------------------------------------------------------------------
# 路径配置（基于脚本自身位置动态推断，无需硬编码绝对路径）
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = _SCRIPT_DIR.parent          # 项目根目录
DATA_DIR = ROOT_DIR / "data"           # 原始数据目录

PATH_COUNTS  = DATA_DIR / "GSE149614_HCC.scRNAseq.S71915.count.txt"
PATH_META    = DATA_DIR / "GSE149614_HCC.metadata.updated.txt"
PATH_CHC20   = DATA_DIR / "CHC20_Visium"
PATH_CHC23   = DATA_DIR / "CHC23_Visium"

# 输出 h5ad 路径
OUT_SCRNA   = DATA_DIR / "scRNA_reference.h5ad"
OUT_CHC20   = DATA_DIR / "chc20_visium.h5ad"
OUT_CHC23   = DATA_DIR / "chc23_visium.h5ad"
OUT_MERGED  = DATA_DIR / "merged_visium.h5ad"


# ---------------------------------------------------------------------------
# 通用工具
# ---------------------------------------------------------------------------

def _require_exists(path: Path, desc: str) -> None:
    """检查路径是否存在，不存在则抛出 FileNotFoundError。"""
    if not path.exists():
        raise FileNotFoundError(f"{desc} 不存在: {path}")


def _normalize_gene_ids(index_obj: pd.Index) -> pd.Index:
    """去掉 Ensembl ID 版本号后缀（如 ENSG00000001234.15 → ENSG00000001234）。"""
    vals = pd.Series(index_obj.astype(str)).str.strip()
    vals = vals.str.replace(r"^(ENSG\d+)\.\d+$",    r"\1", regex=True)
    vals = vals.str.replace(r"^(ENSMUSG\d+)\.\d+$", r"\1", regex=True)
    vals = vals.str.replace(r"^(ENST\d+)\.\d+$",    r"\1", regex=True)
    return pd.Index(vals)


def _guess_id_style(index_obj: pd.Index) -> str:
    """粗略判断基因索引是 Ensembl 风格还是 Symbol 风格。"""
    if len(index_obj) == 0:
        return "empty"
    sample = pd.Series(index_obj.astype(str)).head(200)
    ratio = sample.str.startswith(("ENSG", "ENSMUSG", "ENST")).mean()
    return "ensembl_like" if ratio > 0.6 else "symbol_or_other"


def _set_ensembl_index(adata: sc.AnnData, name: str,
                       gene_id_column: str | None = None) -> sc.AnnData:
    """
    将 AnnData 的基因索引（var_names）切换为 Ensembl ID 并去除版本号。

    参数
    ----
    adata          : 待处理的 AnnData 对象
    name           : 数据集名称，用于报错信息
    gene_id_column : adata.var 中存放 Ensembl ID 的列名；
                     若为 None，则认为 var_names 本身已是 Ensembl ID
    """
    if gene_id_column is not None:
        if gene_id_column not in adata.var.columns:
            raise KeyError(f"{name} 缺少 `{gene_id_column}` 列，无法切换到 Ensembl 索引。")
        # 先把当前 Symbol 索引备份到 SYMBOL 列，方便人工回溯
        adata.var["SYMBOL"] = adata.var_names.astype(str)
        adata.var.set_index(gene_id_column, drop=True, inplace=True)

    # 去除版本号、确保唯一
    adata.var_names = _normalize_gene_ids(adata.var_names)
    adata.var_names_make_unique()

    if _guess_id_style(adata.var_names) != "ensembl_like":
        raise ValueError(
            f"{name} 的 var_names 不是 Ensembl 风格，"
            "请确认 gene_id_column 参数或输入文件是否正确。"
        )
    return adata


def print_summary(adata: sc.AnnData, name: str) -> None:
    """打印 AnnData 的关键结构信息，用于转换后的快速核验。"""
    print(f"\n{'='*50}")
    print(f"  {name}")
    print(f"{'='*50}")
    print(f"  cells / spots : {adata.n_obs}")
    print(f"  genes         : {adata.n_vars}")
    print(f"  obs 列        : {list(adata.obs.columns)}")
    print(f"  var 列        : {list(adata.var.columns)}")
    print(f"  var_names 风格: {_guess_id_style(adata.var_names)}")
    print(f"  var_names 示例: {list(adata.var_names[:5])}")
    if "sample" in adata.obs.columns:
        print(f"  样本分布      :\n{adata.obs['sample'].value_counts().to_string()}")
    print()


# ---------------------------------------------------------------------------
# 转换任务 1：scRNA-seq txt → h5ad
# ---------------------------------------------------------------------------

def convert_scrna(
    path_counts: Path = PATH_COUNTS,
    path_meta: Path   = PATH_META,
    out_path: Path    = OUT_SCRNA,
) -> sc.AnnData:
    """
    将 GSE149614 txt 格式的 count 矩阵 + metadata 转换并保存为 h5ad。

    转换内容
    --------
    - 读取 tab 分隔的 count 矩阵（基因×细胞），转置为 AnnData（细胞×基因）
    - 对齐 metadata，将 celltype 等注释写入 adata.obs
    - 将 var_names 标准化为 Ensembl ID（如原始数据已是 Symbol，则跳过此步，
      保留 Symbol 索引并在 var['SYMBOL'] 中备份）
    - 保存为 h5ad

    注意
    ----
    GSE149614 的 count 矩阵 var_names 本身为 Gene Symbol（如 GAPDH），
    不含 Ensembl ID 列。若后续需要与 Visium（Ensembl 索引）对齐，
    需要额外的 Symbol→Ensembl 映射文件，此转换脚本不处理该映射，
    由 run_preprocessing.py 中的 align_shared_genes() 负责。
    """
    _require_exists(path_counts, "scRNA count 矩阵")
    _require_exists(path_meta,   "scRNA metadata")

    print("正在读取 scRNA count 矩阵（7 万细胞级别，约需 1–2 分钟）...")
    counts = pd.read_csv(path_counts, sep="\t", index_col=0)
    # count 矩阵原始排列为 基因×细胞，AnnData 需要 细胞×基因，故转置
    adata = ad.AnnData(counts.T)

    print("正在读取并对齐 metadata...")
    meta = pd.read_csv(path_meta, sep="\t", index_col=0)
    common_cells = adata.obs_names.intersection(meta.index)
    if len(common_cells) == 0:
        raise ValueError("count 矩阵与 metadata 没有共同细胞条形码，请检查两份文件的行/列索引。")

    adata = adata[common_cells].copy()
    adata.obs = meta.loc[common_cells].copy()

    # 将当前 Symbol 索引备份，保留可读性
    adata.var["SYMBOL"] = adata.var_names.astype(str)
    adata.var_names_make_unique()

    print(f"转换完成：{adata.n_obs} 个细胞，{adata.n_vars} 个基因")
    if "celltype" in adata.obs.columns:
        print("细胞类型分布：")
        print(adata.obs["celltype"].value_counts().to_string())

    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    print(f"已保存: {out_path}")
    print_summary(adata, "scRNA_reference")
    return adata


# ---------------------------------------------------------------------------
# 转换任务 2：Visium Space Ranger 目录 → h5ad
# ---------------------------------------------------------------------------

def convert_visium_slice(
    path_visium: Path,
    sample_name: str,
    out_path: Path,
) -> sc.AnnData:
    """
    将单张 Visium Space Ranger 输出目录转换并保存为 h5ad。

    转换内容
    --------
    - 用 scanpy 读取 Visium 目录（filtered_feature_bc_matrix + spatial）
    - 在 adata.obs 中添加 sample 列标记样本来源
    - 将 var_names 备份为 SYMBOL，切换为 Ensembl ID 索引（来自 gene_ids 列）
    - 将原始 count 矩阵额外存入 layers["counts"]，避免后续归一化操作覆盖原始数据
    - 保存为 h5ad
    """
    _require_exists(path_visium, f"{sample_name} Visium 目录")

    print(f"正在读取 Visium 切片: {sample_name} ...")
    adata = sc.read_visium(path_visium)
    adata.var_names_make_unique()

    # 记录样本来源，多切片合并后可通过此列区分
    adata.obs["sample"] = sample_name

    # Space Ranger 的 features.tsv.gz：第二列（Gene Symbol）默认为 var_names，
    # 第一列（Ensembl ID）存放在 adata.var["gene_ids"]
    adata = _set_ensembl_index(adata, sample_name, gene_id_column="gene_ids")

    # 保存原始 count（整数矩阵），供 Cell2location 等需要原始 count 的工具使用
    adata.layers["counts"] = adata.X.copy()

    print(f"转换完成：{adata.n_obs} 个 spot，{adata.n_vars} 个基因")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)
    print(f"已保存: {out_path}")
    print_summary(adata, sample_name)
    return adata


# ---------------------------------------------------------------------------
# 转换任务 3（可选）：合并多张 Visium 切片 → merged_visium.h5ad
# ---------------------------------------------------------------------------

def merge_visium_slices(
    slices: list[tuple[Path, str, Path]],
    out_path: Path = OUT_MERGED,
) -> sc.AnnData:
    """
    将多张已转换的 Visium h5ad 文件合并为单一 AnnData，保存为 merged_visium.h5ad。

    参数
    ----
    slices   : 列表，每个元素为 (h5ad路径, 样本名, 原始Space Ranger目录)；
               若 h5ad 文件不存在，会先调用 convert_visium_slice() 生成
    out_path : 合并结果的输出路径

    合并策略
    --------
    - join="inner"：只保留所有切片共有的基因（Ensembl ID 交集）
    - index_unique="-"：为 spot barcode 加后缀（如 AAACAATCTACTAGTT-1-CHC20），
      避免不同切片间 barcode 重名导致索引冲突
    """
    adatas = []
    keys = []
    for h5ad_path, name, raw_path in slices:
        if h5ad_path.exists():
            print(f"从缓存加载: {h5ad_path}")
            adata = sc.read_h5ad(h5ad_path)
        else:
            adata = convert_visium_slice(raw_path, name, h5ad_path)
        adatas.append(adata)
        keys.append(name)

    print(f"\n正在合并 {len(adatas)} 张切片（取基因交集）...")
    merged = ad.concat(
        adatas,
        join="inner",        # 只保留所有切片共有的基因
        label="sample",
        keys=keys,
        index_unique="-",    # spot barcode 加后缀防重名
        merge="same",
    )
    # 合并后重新写入 counts layer（concat 会保留各切片的 layer，但需确认）
    if "counts" not in merged.layers:
        merged.layers["counts"] = merged.X.copy()

    gene_counts = [a.n_vars for a in adatas]
    print(f"合并后：{merged.n_obs} 个 spot，{merged.n_vars} 个基因")
    print(f"（各切片基因数：{dict(zip(keys, gene_counts))}，交集后保留 {merged.n_vars} 个）")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out_path)
    print(f"已保存: {out_path}")
    print_summary(merged, "merged_visium")
    return merged


# ---------------------------------------------------------------------------
# 命令行入口
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="原始数据格式转换：txt / Visium 目录 → h5ad（不做训练或归一化）"
    )
    parser.add_argument("--skip-scrna",   action="store_true", help="跳过 scRNA txt→h5ad 转换")
    parser.add_argument("--skip-visium",  action="store_true", help="跳过 Visium 目录→h5ad 转换")
    parser.add_argument("--no-merge",     action="store_true", help="不生成合并的 merged_visium.h5ad")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    errors: list[str] = []

    # --- 任务 1：scRNA txt → h5ad ---
    if not args.skip_scrna:
        print("\n[任务 1/3] scRNA-seq txt → h5ad")
        try:
            convert_scrna()
        except Exception as e:
            print(f"[警告] scRNA 转换失败，已跳过：{e}", file=sys.stderr)
            errors.append(f"scRNA: {e}")
    else:
        print("\n[任务 1/3] 已跳过 scRNA 转换（--skip-scrna）")

    # --- 任务 2：Visium 目录 → 各自的 h5ad ---
    if not args.skip_visium:
        print("\n[任务 2/3] Visium Space Ranger 目录 → h5ad")
        for raw_path, name, out_path in [
            (PATH_CHC20, "CHC20", OUT_CHC20),
            (PATH_CHC23, "CHC23", OUT_CHC23),
        ]:
            try:
                convert_visium_slice(raw_path, name, out_path)
            except Exception as e:
                print(f"[警告] {name} Visium 转换失败，已跳过：{e}", file=sys.stderr)
                errors.append(f"{name}: {e}")
    else:
        print("\n[任务 2/3] 已跳过 Visium 转换（--skip-visium）")

    # --- 任务 3（可选）：合并多张 Visium 切片 ---
    if not args.no_merge and not args.skip_visium:
        print("\n[任务 3/3] 合并多张 Visium 切片 → merged_visium.h5ad")
        slices = [
            (OUT_CHC20, "CHC20", PATH_CHC20),
            (OUT_CHC23, "CHC23", PATH_CHC23),
        ]
        try:
            merge_visium_slices(slices)
        except Exception as e:
            print(f"[警告] Visium 合并失败，已跳过：{e}", file=sys.stderr)
            errors.append(f"merge: {e}")
    else:
        print("\n[任务 3/3] 已跳过 Visium 合并（--no-merge 或 --skip-visium）")

    # --- 汇总 ---
    print("\n" + "="*50)
    if errors:
        print(f"完成（{len(errors)} 项任务失败，详见上方警告）：")
        for err in errors:
            print(f"  • {err}")
        return 1
    else:
        print("所有转换任务完成，h5ad 文件已保存至 data/ 目录。")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
