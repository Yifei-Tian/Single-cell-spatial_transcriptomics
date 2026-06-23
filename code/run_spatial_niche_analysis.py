from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import cast

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.sparse as sp
from scipy.spatial import KDTree


# ============================================================
# 免疫抑制相关基因列表（全局常量）
# 这些基因在免疫抑制微环境中高表达，是识别免疫抑制区域的核心标志基因：
#   FOXP3       - Treg 细胞的主转录因子，Treg 最核心的标志基因
#   IL2RA       - IL-2 受体 α 链（CD25），Treg 细胞表面标志物
#   CTLA4       - 抑制性免疫检查点，阻断 T 细胞激活
#   TIGIT       - T 细胞抑制性受体，介导免疫耗竭
#   LAG3        - 免疫检查点分子，抑制 T 细胞和 NK 细胞功能
#   PDCD1       - 即 PD-1，免疫检查点靶点
#   HAVCR2      - 即 TIM-3，促进 T 细胞耗竭
#   TGFB1       - 转化生长因子β1，强免疫抑制细胞因子
#   IL10        - 白介素10，抑制炎症反应的细胞因子
#   CXCL12      - 趋化因子，招募 Treg 和 MDSC 进入肿瘤
#   CCL22       - 趋化因子，特异性招募 Treg 细胞
#   TNFRSF18    - 即 GITR，Treg 激活共刺激分子
#   TNFRSF4     - 即 OX40，调节 Treg 功能
#   IKZF2       - 即 Helios，稳定 Treg 表型的转录因子
# ============================================================
IMMUNOSUPPRESSIVE_GENES = (
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


def _parse_args() -> argparse.Namespace:
    """
    解析命令行参数，为整个分析流程提供所有配置选项。

    【作用】
    使用 argparse 库定义脚本支持的所有参数，并自动处理用户在终端输入的参数。
    每个参数都有默认值，因此不传参数也可直接运行。

    【参数说明】
    --adata                   : 输入的空间转录组 .h5ad 文件路径
    --out-dir                 : 结果输出目录路径
    --signature-out           : 特征基因列表的输出路径（供 TCGA 分析使用）
    --abundance-key           : 在 adata.obsm 中存储细胞丰度矩阵的键名
    --hepatocyte-col          : 肝细胞（Hepatocyte）对应的列名
    --treg-col                : 调节性 T 细胞（Treg）对应的列名
    --myeloid-col             : 髓系细胞（Myeloid）对应的列名
    --fibroblast-col          : 成纤维细胞（Fibroblast）对应的列名
    --tnk-col                 : T/NK 细胞对应的列名
    --hep-high-quantile       : 定义"高肝细胞"区域的分位数阈值（默认 0.75，即前25%）
    --niche-high-quantile     : 定义"高免疫抑制 niche"的分位数阈值（默认 0.80，即前20%）
    --neighbor-radius-multiplier : 空间邻居搜索半径的倍增系数（默认 1.25）
    --top-niche-genes         : 输出 niche 特征基因的数量（默认 80 个）
    --seed                    : 随机数种子，确保结果可复现（默认 1234）

    【返回值】
    argparse.Namespace 对象，可以通过 args.xxx 访问每个参数的值。
    """
    # 找到项目根目录（当前文件向上两级）
    root = Path(__file__).resolve().parents[1]
    results = root / "results"
    parser = argparse.ArgumentParser(
        description="Define an immunosuppressive spatial niche from cell2location output."
    )
    parser.add_argument("--adata", type=Path, default=results / "adata_vis_post.h5ad") 
    parser.add_argument("--out-dir", type=Path, default=results / "spatial_niche")
    parser.add_argument("--signature-out", type=Path, default=results / "spatial_signature_genes.txt")
    parser.add_argument("--abundance-key", default="means_cell_abundance_w_sf")
    parser.add_argument("--hepatocyte-col", default="Hepatocyte")
    parser.add_argument("--treg-col", default="Treg")
    parser.add_argument("--myeloid-col", default="Myeloid")
    parser.add_argument("--fibroblast-col", default="Fibroblast")
    parser.add_argument("--tnk-col", default="T/NK")
    parser.add_argument("--hep-high-quantile", type=float, default=0.75)
    parser.add_argument("--niche-high-quantile", type=float, default=0.80)
    parser.add_argument("--neighbor-radius-multiplier", type=float, default=1.25)
    parser.add_argument("--top-niche-genes", type=int, default=80)
    parser.add_argument("--seed", type=int, default=1234)
    return parser.parse_args()


def _setup_logging() -> None:
    """
    初始化全局日志系统配置。

    【作用】
    在脚本启动时调用一次，配置日志的输出格式、级别和目标流，
    使后续所有 logging.info() / logging.warning() 的调用都按统一格式输出。

    【配置说明】
    - level=logging.INFO  : 只显示 INFO 及以上级别（INFO/WARNING/ERROR/CRITICAL），
                            DEBUG 级别的调试信息会被过滤掉
    - format=...          : 日志格式，每条日志包含三部分：
                              %(asctime)s   → 时间戳，如 2026-06-06 10:23:45,123
                              %(levelname)s → 日志级别，如 INFO / WARNING
                              %(message)s   → 日志正文内容
                            示例输出：2026-06-06 10:23:45,123 | INFO | Loading spatial AnnData...
    - handlers            : 日志处理器列表，这里使用 StreamHandler(sys.stdout)，
                            将日志输出到标准输出（终端屏幕），而非默认的标准错误（stderr）
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def _clean_abundance_columns(columns: pd.Index) -> list[str]:
    """
    去除 cell2location 后验矩阵列名中的统计量前缀，只保留细胞类型名称。

    【背景知识】
    cell2location 是一种空间转录组去卷积工具，它的输出矩阵列名格式为：
      "<统计量前缀>_<细胞类型名>"
    例如：
      means_cell_abundance_w_sf_Hepatocyte  →  Hepatocyte（肝细胞）
      q05_cell_abundance_w_sf_T/NK          →  T/NK（T/NK 细胞）
      means_per_cluster_mu_fg_Treg          →  Treg（调节性T细胞）

    【处理步骤】
    使用正则表达式（re.sub）匹配并删除两种常见前缀：
      1. "means_cell_abundance_w_sf_"  （或 q05_、median_ 等变体）
      2. "means_per_cluster_mu_fg_"    （RegressionModel 后验分布）

    【参数】
    columns : pandas 列索引，原始的带前缀列名

    【返回值】
    list[str]，清洗后的纯细胞类型名称列表
    """
    import re
    cleaned = []
    for col in columns:
        value = str(col)
        # 匹配所有统计量前缀：means / q<数字> / median，后接可选下划线和 cell_abundance_w_sf_
        value = re.sub(r"^(means|q\d+|median)_?cell_abundance_w_sf_", "", value)
        # 匹配 RegressionModel 后验前缀
        value = re.sub(r"^means_per_cluster_mu_fg_", "", value)
        cleaned.append(value)
    return cleaned


def _get_abundance(adata: ad.AnnData, key: str) -> pd.DataFrame:
    """
    从 AnnData 对象中提取细胞类型丰度矩阵，并清洗列名。

    【背景知识】
    AnnData 是单细胞/空间转录组领域通用的数据格式：
      adata.obsm  : 存储每个 spot/细胞的附加矩阵（如细胞丰度、降维坐标等）
      adata.obs_names : 每个 spot 的唯一 ID 列表
      adata.uns   : 非结构化附加信息（如模型因子名称）

    【处理流程】
    1. 检查 key 是否存在于 adata.obsm 中，不存在则报错并列出可用键
    2. 提取丰度矩阵（可能是 DataFrame 或 numpy 数组两种格式）
    3. 如果是 numpy 数组，则从 adata.uns["mod"]["factor_names"] 获取列名
    4. 统一设置行索引为 spot ID，并调用 _clean_abundance_columns 清洗列名

    【参数】
    adata : AnnData 空间转录组数据对象
    key   : adata.obsm 中细胞丰度矩阵的键名（如 "means_cell_abundance_w_sf"）

    【返回值】
    pd.DataFrame，行为 spot，列为细胞类型名称，值为对应丰度数值
    """
    if key not in adata.obsm:
        raise KeyError(f"Missing adata.obsm['{key}']; available keys: {list(adata.obsm.keys())}")
    abundance = adata.obsm[key]
    if isinstance(abundance, pd.DataFrame):
        # 已经是 DataFrame，直接复制使用
        df = abundance.copy()
    else:
        # 是 numpy 数组，需要从 adata.uns 中获取列名（细胞类型名称）
        factors = adata.uns.get("mod", {}).get("factor_names")
        df = pd.DataFrame(abundance, index=adata.obs_names, columns=factors)
    # 统一行索引为 spot ID
    df.index = adata.obs_names
    # 清洗列名，去掉 cell2location 的统计量前缀
    df.columns = _clean_abundance_columns(df.columns)
    return df


def _zscore(series: pd.Series) -> pd.Series:
    """
    对一列数据进行 Z-score 标准化（零均值、单位方差）。

    【数学公式】
    z = (x - mean) / std
    其中：
      x    : 原始数值
      mean : 该列所有值的均值
      std  : 该列所有值的标准差

    【作用】
    Z-score 标准化使不同量纲、不同数值范围的指标可以直接相加比较。
    例如，Treg 比例（0~1）和基因评分（可能0~10）经过 Z-score 后都变成
    以0为中心、标准差为1的无量纲数值，相加时贡献相等。

    【特殊情况处理】
    如果标准差为 0（所有值相同）或无效值（NaN/Inf），
    则返回全 0 的 Series，避免除以零的报错。

    【参数】
    series : 待标准化的数值列（pd.Series）

    【返回值】
    pd.Series，标准化后的数值，与输入保持相同的索引
    """
    std = float(series.std())
    # 如果标准差无效或为0，说明该列没有区分度，直接返回全0
    if not np.isfinite(std) or std == 0:
        return pd.Series(0.0, index=series.index)
    return (series - float(series.mean())) / std


def _expression_frame(adata: ad.AnnData) -> pd.DataFrame:
    """
    将 AnnData 中的原始基因表达矩阵转换为标准化后的表达量 DataFrame。

    【标准化流程（CP10K + log1p）】
    这是单细胞/空间转录组分析中最常用的表达量归一化方式：
    Step 1 - CP10K 归一化：
      每个 spot 的原始计数除以该 spot 的总计数，再乘以 10000。
      目的：消除不同 spot 之间测序深度的差异（测序越深，计数越多，不可直接比较）
      公式：normalized = (count / total_count) × 10000
    Step 2 - log1p 对数变换：
      对 CP10K 值做 log(x+1) 变换，加1是为了避免对0取对数。
      目的：压缩高表达基因的数值范围，使表达值分布更接近正态分布

    【稀疏矩阵处理】
    单细胞数据通常以稀疏矩阵（scipy.sparse）存储（因为大多数基因不表达，值为0）。
    代码会自动检测并分支处理稀疏和稠密两种格式。

    【参数】
    adata : AnnData 数据对象，.X 属性存储原始计数矩阵

    【返回值】
    pd.DataFrame，行为 spot ID，列为基因名，值为 log1p(CP10K) 标准化表达量
    """
    x = adata.X
    if x is None:
        raise ValueError("adata.X is empty, cannot normalize expression matrix.")

    if sp.issparse(x): # 如果是稀疏矩阵
        # ---- 稀疏矩阵处理分支 ----
        x = cast(sp.spmatrix, x).copy().astype(float).tocsr()  # 转为 CSR 格式，便于行操作
        totals = np.asarray(x.sum(axis=1)).ravel()  # 每个 spot 的总计数
        # 计算缩放因子 10000/total，对总计数为0的 spot 缩放因子设为0（避免除零）
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        # 用对角矩阵乘法实现逐行缩放（等价于每行乘以对应的 scale 值）
        x = sp.diags(scale) @ x
        # 只对非零元素做 log1p（稀疏矩阵的0值不需要变换，log1p(0)=0 仍为0）
        x.data = np.log1p(x.data)
        x = x.toarray()  # 最终转为稠密数组
    else:
        # ---- 稠密矩阵处理分支 ----
        x = np.array(x, dtype=float, copy=True)
        totals = x.sum(axis=1)  # 每个 spot 的总计数
        scale = np.divide(1e4, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        # scale[:, None] 将1D数组扩展为列向量，实现逐行广播乘法
        x = np.log1p(x * scale[:, None])
    return pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names)


def _module_score(adata: ad.AnnData, genes: tuple[str, ...]) -> tuple[pd.Series, list[str]]:
    """
    计算每个 spot 的免疫抑制基因模块评分。

    【什么是基因模块评分】
    基因模块评分（Module Score）是一种衡量某个基因集在每个细胞/spot 中
    整体表达水平的方法。这里简单地取目标基因集的平均表达量作为评分。
    评分越高，说明该 spot 的免疫抑制基因整体表达越强。

    【处理流程】
    1. 过滤：从 IMMUNOSUPPRESSIVE_GENES 中只保留在当前数据集中实际存在的基因
       （不同数据集可能缺少某些基因，需要跳过以避免报错）
    2. 如果一个基因都没有，发出警告并返回全0评分
    3. 调用 _expression_frame 获取标准化表达矩阵
    4. 对所有可用基因的表达量按行取均值，得到每个 spot 的评分

    【参数】
    adata : AnnData 数据对象
    genes : 要计算评分的目标基因名称元组（使用全局常量 IMMUNOSUPPRESSIVE_GENES）

    【返回值】
    tuple，包含：
      - pd.Series : 每个 spot 的免疫抑制基因模块评分，名称为 "immunosuppressive_gene_score"
      - list[str] : 实际被使用的基因名称列表（用于记录到元数据）
    """
    # 只保留在数据集中实际存在的基因
    available = [gene for gene in genes if gene in adata.var_names]
    if not available:
        # 没有任何匹配基因，发出警告，返回全0评分（不影响后续流程）
        logging.warning("No immunosuppressive marker genes found in adata.var_names.")
        return pd.Series(0.0, index=adata.obs_names, name="immunosuppressive_gene_score"), []
    # 获取标准化表达矩阵，只取目标基因列，按行求均值
    expr = _expression_frame(adata)
    score = expr[available].mean(axis=1)
    score.name = "immunosuppressive_gene_score"
    return score, available


def _spatial_neighbors(coords: np.ndarray, radius_multiplier: float) -> tuple[list[np.ndarray], float]:
    """
    基于空间坐标计算每个 spot 的空间邻居列表，并自动确定搜索半径。

    【什么是空间邻居】
    在空间转录组中，每个 spot 有 x/y 坐标。空间邻居是指在物理空间上
    距离足够近的其他 spot，它们共同构成该 spot 的局部微环境。

    【搜索半径的自动计算】
    1. 用 cKDTree（k-d 树，一种空间索引数据结构）构建所有 spot 的空间索引
    2. 查询每个 spot 到"最近邻居"（k=2，第1个是自身，取第2个）的距离
    3. 取所有最近邻距离的中位数，乘以 radius_multiplier 作为搜索半径
       - 中位数比均值更稳健，不受少数异常远离 spot 的影响
       - radius_multiplier > 1 确保能搜索到至少一个邻居

    【邻居清洗】
    query_ball_point 返回的结果包含 spot 自身（距离为0），
    需要从邻居列表中排除自身索引 i。

    【参数】
    coords            : 形状为 (n_spots, 2) 的二维坐标数组，每行为一个 spot 的 (x, y)
    radius_multiplier : 搜索半径倍增系数，默认 1.25

    【返回值】
    tuple，包含：
      - list[np.ndarray] : 长度为 n_spots 的列表，每个元素是该 spot 的邻居索引数组
      - float            : 使用的搜索半径值（记录到元数据供参考）
    """
    # 构建 k-d 树空间索引，实现高效的距离查询
    tree = KDTree(coords)
    # 查询每个 spot 的最近2个点（第0个是自身，第1个是最近邻居），取距离列
    nearest = tree.query(coords, k=2)[0][:, 1]
    # 以最近邻距离中位数 × 倍增系数 作为搜索半径
    radius = float(np.median(nearest) * radius_multiplier)
    # 找出每个 spot 在半径内的所有 spot（包括自身）
    neighbor_lists = tree.query_ball_point(coords, r=radius)
    # 从每个邻居列表中排除 spot 自身
    cleaned = []
    for i, indices in enumerate(neighbor_lists):
        arr = np.array([idx for idx in indices if idx != i], dtype=int)
        cleaned.append(arr)
    return cleaned, radius


def _neighbor_mean(values: pd.Series, neighbors: list[np.ndarray]) -> pd.Series:
    """
    计算每个 spot 在其空间邻居范围内的平均值。通过周围邻居spot细胞类型的平均比例判断该spot位置和周围spot类型

    【作用】
    对于某个细胞类型比例（如 Treg 比例），计算每个 spot 周围邻居 spot 的平均比例，
    用于捕捉局部微环境的特征。例如，一个 spot 自身 Treg 比例低，
    但周围 spot Treg 比例高，说明该 spot 处于 Treg 富集区域边缘。

    【处理细节】
    - 没有邻居的 spot（孤立点）输出 NaN
    - 使用 numpy 数组操作代替 Python 循环内的 pandas 操作，提升计算效率

    【参数】
    values    : 待平均的数值列（pd.Series），如某细胞类型的比例
    neighbors : 每个 spot 的邻居索引数组列表（由 _spatial_neighbors 生成）

    【返回值】
    pd.Series，每个 spot 对应其邻居的平均值，无邻居的 spot 为 NaN，
    与输入 values 保持相同的索引
    """
    arr = values.to_numpy()  # 转为 numpy 数组；value（所有spot的treg比例）；arr为numpy数组，后面提取索引更快
    out = np.full(len(values), np.nan, dtype=float)  
    # np.full() 创建一个填充了指定值的数组；参数为长度、填充值和数据类型
    # len 创建数组的长度
    for i, idx in enumerate(neighbors): # neighbors 储存着对应spot所有邻居的索引编号，这里同时获取序号i和值idx
        if len(idx):  # 邻居数组的长度
            out[i] = float(np.mean(arr[idx]))
    return pd.Series(out, index=values.index)


def _distance_to_mask(coords: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """
    计算每个 spot 到"目标区域"（mask 为 True 的 spot 集合）的最短距离。

    【作用】
    用于计算每个 spot 距离"高肝细胞区域"（肿瘤核心）的最短空间距离。
    这个距离可以反映 spot 与肿瘤核心的位置关系：
      - 距离接近0 → 在肿瘤核心内部或紧邻
      - 距离较大  → 远离肿瘤核心，处于基质区

    【参数】
    coords : 形状为 (n_spots, 2) 的空间坐标数组
    mask   : 布尔数组，True 表示该 spot 属于目标区域（高肝细胞 spot）

    【返回值】
    np.ndarray，每个 spot 到目标区域最近 spot 的距离。
    如果目标区域为空（没有高肝细胞 spot），返回全 NaN 数组。
    """
    if not bool(mask.any()):
        # 目标区域为空，无法计算距离，返回全 NaN
        return np.full(coords.shape[0], np.nan)
    # 只用目标区域的 spot 坐标构建 k-d 树
    tree = KDTree(coords[mask])
    # 查询每个 spot 到目标区域最近点的距离（k=1 表示只要最近的1个）
    return tree.query(coords, k=1)[0]


def _spatial_scatter(
    df: pd.DataFrame,
    value: str,
    path: Path,
    title: str,
    cmap: str = "viridis",
    categorical: bool = False,
) -> None:
    """
    绘制空间散点图，将每个 spot 按其坐标位置绘制，颜色表示某一特征的数值。

    【图表说明】
    这是空间转录组分析中最基本的可视化方式：
    - x 轴：spot 的空间 x 坐标
    - y 轴：spot 的空间 y 坐标（y 轴翻转以匹配组织切片的通常方向）
    - 颜色：该 spot 某特征的值（连续值用颜色渐变，分类值用固定颜色）

    【两种模式】
    1. 连续值模式（categorical=False，默认）：
       使用颜色映射（cmap）将数值映射为颜色，右侧显示颜色条
       适用于：细胞比例、niche 评分等连续数值
    2. 分类值模式（categorical=True）：
       每个类别使用固定颜色，右侧显示图例
       适用于：空间区域标签（tumor_core / tumor_edge / stroma_immune / other）
       颜色方案：红色=肿瘤核心, 橙色=肿瘤边缘, 蓝色=免疫基质, 灰色=其他

    【参数】
    df          : 包含 spatial_x、spatial_y 及目标值列的 DataFrame
    value       : 用于着色的列名
    path        : 图片保存路径（.png 格式）
    title       : 图表标题
    cmap        : 连续值模式下使用的颜色映射名称（默认 "viridis"，绿-黄-蓝渐变）
    categorical : 是否为分类变量模式，默认 False
    """
    fig, ax = plt.subplots(figsize=(7, 6))
    if categorical:
        # 分类变量：定义固定颜色映射
        palette = {
            "tumor_core": "#d62728",    # 红色 - 肿瘤核心
            "tumor_edge": "#ff7f0e",    # 橙色 - 肿瘤边缘
            "stroma_immune": "#1f77b4", # 蓝色 - 免疫基质区
            "other": "#bdbdbd",         # 灰色 - 其他区域
        }
        # 将每个 spot 的分类标签映射为对应颜色，未知类别填充灰色
        colors = df[value].map(palette).fillna("#bdbdbd")
        ax.scatter(df["spatial_x"], df["spatial_y"], c=colors, s=18, linewidths=0)
        # 手动创建图例
        handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=color, label=label, markersize=8)
            for label, color in palette.items()
        ]
        ax.legend(handles=handles, frameon=False, loc="best")
    else:
        # 连续变量：使用颜色映射渐变着色
        sc = ax.scatter(
            df["spatial_x"],
            df["spatial_y"],
            c=df[value],   # 连续数值自动映射到颜色
            s=18,          # 每个点的大小（像素面积）
            cmap=cmap,     # 颜色映射方案
            linewidths=0,  # 不画点的边框线，更美观
        )
        # 添加颜色条（色标尺），显示数值范围
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title)
    ax.set_xlabel("spatial x")
    ax.set_ylabel("spatial y")
    ax.invert_yaxis()  # 翻转 y 轴，使图像方向与显微镜切片一致
    ax.set_aspect("equal", adjustable="box")  # 保持 x/y 轴比例相等，避免图像变形
    fig.tight_layout()  # 自动调整布局，防止标题/坐标轴被裁剪
    fig.savefig(path, dpi=180)  # 以 180 DPI 保存高分辨率图片
    plt.close(fig)  # 关闭图形对象，释放内存


def _plot_distance(df: pd.DataFrame, path: Path) -> None:
    """
    绘制"距离高肝细胞区域的距离"与"免疫抑制 niche 评分"的关系折线图。

    【图表说明】
    将每个 spot 到高肝细胞区域的距离分成5个区间（分位数分箱），
    计算每个区间内 niche 评分的均值 ± 标准误差（SEM），绘制折线+误差棒图。

    【生物学意义】
    该图用于验证：距离肿瘤核心越近的 spot，其免疫抑制 niche 评分是否越高。
    如果存在明显的距离依赖性趋势，说明免疫抑制微环境主要分布在肿瘤核心周围。

    【参数】
    df   : 包含 distance_to_hep_high 和 immunosuppressive_niche_score 列的 DataFrame
    path : 图片保存路径
    """
    # 去掉距离或评分缺失的行
    use = df.dropna(subset=["distance_to_hep_high", "immunosuppressive_niche_score"]).copy() # 数据清理，从原始数据框df这两行中删除含有NaN值的行
    if use.empty:
        return  # 数据为空则跳过绘图
    # 将距离按分位数分成5个等频率区间（每组 spot 数量大致相同）
    use["distance_bin"] = pd.qcut(use["distance_to_hep_high"], q=5, duplicates="drop")
    # 数据分箱（离散化）。使用 pandas 的 qcut 按分位数将连续的距离数据分成 5 个区间（即每个区间包含的数据点数量大致相等）。
    # duplicates="drop" 用于丢弃由于数据中有大量重复值可能导致的边界相同的区间。
    # 按距离区间分组，计算均值、标准误差、样本量
    summary = (
        use.groupby("distance_bin", observed=True)["immunosuppressive_niche_score"]
        .agg(["mean", "sem", "count"])
        .reset_index() # 将分组索引重新恢复为普通列
    ) # 按照distance_bin 列进行分组，计算组内immnosuppressive_niche_score列的均值、标准误差和样本量
    summary["bin"] = np.arange(1, len(summary) + 1)  # 用整数序号作为 x 轴
    fig, ax = plt.subplots(figsize=(7, 4.5)) # 初始化 matplotlib 的图形对象 (fig) 和坐标轴对象 (ax)，画布尺寸为 7 x 4.5 英寸
    # 绘制折线 + 误差棒（error bar）
    ax.errorbar(summary["bin"], summary["mean"], yerr=summary["sem"], marker="o", capsize=3)
    ax.set_xlabel("Distance bin from Hepatocyte-high area") # x 轴标签
    ax.set_ylabel("Mean immunosuppressive niche score") # y 轴标签
    ax.set_xticks(summary["bin"]) # 指定 x 轴的刻度位置必须要对齐刚才设置的整数坐标
    # x 轴刻度标签使用实际的距离区间范围，旋转30度防止重叠
    ax.set_xticklabels([str(v) for v in summary["distance_bin"]], rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_region_box(df: pd.DataFrame, path: Path) -> None:
    """
    绘制不同空间区域的评分箱线图，并排展示两个评分指标。

    【图表说明】
    并排绘制两个箱线图：
    - 左图：各空间区域的 Treg_like_score（Treg 样评分）分布
    - 右图：各空间区域的 immunosuppressive_niche_score（免疫抑制 niche 评分）分布

    区域顺序固定为：tumor_core → tumor_edge → stroma_immune → other
    这个顺序反映了从肿瘤核心向外扩展的空间梯度。

    【生物学意义】
    用于验证空间区域标注的合理性：预期 tumor_edge 和 stroma_immune 区域
    的免疫抑制评分应高于 tumor_core 和 other，体现肿瘤边缘的免疫调控特征。

    【参数】
    df   : 包含 spatial_region、Treg_like_score、immunosuppressive_niche_score 列的 DataFrame
    path : 图片保存路径
    """
    order = ["tumor_core", "tumor_edge", "stroma_immune", "other"] # 定义空间区域类别的显示顺序。
    # 肿瘤核心区域：tumor_core，肿瘤边缘区域：tumor_edge，免疫基质区域：stroma_immune，其他区域：other
    # 只保留有有效区域标签的行
    use = df[df["spatial_region"].isin(order)].copy() # 过滤数据：isin(order) 只保留包含了上述四种有效区域标签的 spot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))  # 1行2列子图
    # 左图：Treg 样评分
    sns.boxplot(data=use, x="spatial_region", y="Treg_like_score", order=order, ax=axes[0])
    # 右图：免疫抑制 niche 评分
    sns.boxplot(data=use, x="spatial_region", y="immunosuppressive_niche_score", order=order, ax=axes[1])
    for ax in axes:
        ax.tick_params(axis="x", rotation=30)  # x 轴标签旋转30度，防止重叠
        ax.set_xlabel("")  # 去掉 x 轴标题，区域名称已足够清晰
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_correlation(df: pd.DataFrame, columns: list[str], path: Path) -> None:
    """
    绘制指定列之间的 Pearson 相关系数热图。

    【图表说明】
    计算选定列（各细胞类型比例 + niche 评分）两两之间的 Pearson 相关系数，
    以热图形式展示。颜色越红表示正相关越强，颜色越蓝表示负相关越强，
    白色表示不相关（相关系数接近0）。每个格子内标注数值（保留2位小数）。

    【生物学意义】
    用于探索细胞类型之间的共定位关系，例如：
    - Treg 与 Myeloid 正相关 → 说明两者倾向于出现在同一空间区域
    - Hepatocyte 与 Treg 负相关 → 说明肝细胞多的区域 Treg 反而少

    【参数】
    df      : 数据 DataFrame
    columns : 需要计算相关性的列名列表
    path    : 图片保存路径
    """
    # 计算 Pearson 相关系数矩阵
    corr = df[columns].corr()
    fig, ax = plt.subplots(figsize=(7, 6))
    # 绘制热图：vlag 是红-白-蓝配色，center=0 使相关系数0对应白色
    sns.heatmap(corr, cmap="vlag", center=0, annot=True, fmt=".2f", square=True, ax=ax)
    ax.set_title("Cell abundance and niche score correlation")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_hep_treg(df: pd.DataFrame, path: Path) -> None:
    """
    绘制肝细胞比例 vs Treg 比例的散点图，颜色编码 niche 评分。

    【图表说明】
    - x 轴：Hepatocyte（肝细胞）比例
    - y 轴：Treg（调节性T细胞）比例
    - 颜色：immunosuppressive_niche_score（免疫抑制 niche 评分），使用 magma 颜色映射
    - 虚线：x 轴 0.75 分位数线（划分高/低肝细胞区域）
             y 轴 0.75 分位数线（划分高/低 Treg 区域）
    - 右侧颜色条：niche 评分的颜色对照

    【生物学意义】
    该图用于直观展示：
    1. 肝细胞和 Treg 在空间上的共存/拮抗关系
    2. 免疫抑制 niche 评分高的 spot 主要分布在散点图的哪个象限
    预期：高 Treg + 低肝细胞的象限（左上角）应有较高的 niche 评分

    【参数】
    df   : 包含 Hepatocyte、Treg、immunosuppressive_niche_score 列的 DataFrame
    path : 图片保存路径
    """
    fig, ax = plt.subplots(figsize=(5.5, 5))
    # 绘制散点，颜色编码 niche 评分，使用 magma（黑-紫-橙-黄）配色
    sc = ax.scatter(
        df["Hepatocyte"],
        df["Treg"],
        c=df["immunosuppressive_niche_score"],
        s=20,
        cmap="magma",
        linewidths=0,
    )
    # 添加 0.75 分位数参考线（虚线），划分四个象限
    ax.axvline(df["Hepatocyte"].quantile(0.75), color="black", linestyle="--", linewidth=1)
    ax.axhline(df["Treg"].quantile(0.75), color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Hepatocyte proportion")
    ax.set_ylabel("Treg proportion")
    fig.colorbar(sc, ax=ax, label="niche score")  # 添加颜色条并标注含义
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _rank_niche_genes(adata: ad.AnnData, labels: pd.Series, top_n: int) -> pd.DataFrame:
    """
    对 niche 高分 spot 进行差异基因分析，筛选出特征性高表达基因。

    【什么是差异基因分析】
    比较"免疫抑制 niche 高分组（niche_high=True）"与"低分组（niche_high=False）"
    之间每个基因的表达量差异，找出在高分组中特异性高表达的基因，
    这些基因构成免疫抑制微生态位的"特征基因签名"（Gene Signature）。

    【排序方法 - log2 Fold Change (log2FC)】
    log2FC = log2(高分组平均表达量 + 1) - log2(低分组平均表达量 + 1)
    加1是为了防止分母为0，+1 是 pseudocount（伪计数）。
    log2FC 越大，说明该基因在高分组中越特异性高表达。

    【参数】
    adata  : AnnData 数据对象，提供基因表达矩阵
    labels : 布尔类型的 pd.Series，True 表示该 spot 为 niche_high（免疫抑制 niche 高分）
    top_n  : 返回排名前 top_n 的差异基因数量

    【返回值】
    pd.DataFrame，包含以下列：
      - gene      : 基因名称
      - mean_high : 高分组中该基因的平均表达量（log1p CP10K）
      - mean_low  : 低分组中该基因的平均表达量（log1p CP10K）
      - log2_fc   : log2 Fold Change，按此列降序排序
    按 log2_fc 降序，返回前 top_n 行（最显著上调基因）
    """
    expr = _expression_frame(adata) # 获取了一个经过标准化的基因表达矩阵，行为 spot，列为基因
    high = labels.astype(bool)  # 确保是布尔类型
    # 安全检查：两组各需至少3个样本才能进行有意义的比较
    if int(high.sum()) < 3 or int((~high).sum()) < 3:
        logging.warning("Too few spots in niche_high/niche_low groups; skipping signature ranking.")
        return pd.DataFrame(columns=["gene", "mean_high", "mean_low", "log2_fc"])
    # 分别计算高分组和低分组中每个基因的平均表达量
    mean_high = expr.loc[high].mean(axis=0)   # 对高分 spot 的行取列均值
    mean_low = expr.loc[~high].mean(axis=0)   # 对低分 spot 的行取列均值
    ranked = pd.DataFrame(
        {
            "gene": expr.columns,
            "mean_high": mean_high.to_numpy(), # .to_numpy() 转成numpy数组
            "mean_low": mean_low.to_numpy(),
        }
    )
    # 计算 log2 Fold Change：高分组 vs 低分组，+1 为伪计数防止除零
    ranked["log2_fc"] = np.log2((ranked["mean_high"] + 1.0) / (ranked["mean_low"] + 1.0))
    # 按 log2FC 降序排列，取前 top_n 个最显著上调基因
    ranked = ranked.sort_values("log2_fc", ascending=False)
    return ranked.head(top_n)


def main() -> int:
    """
    主函数：协调所有分析步骤，完成空间免疫抑制微生态位（Niche）的全流程分析。

    【分析流程概述】
    1. 初始化：设置日志、解析参数、创建输出目录
    2. 数据加载与校验：读取空间转录组数据，验证必要字段存在
    3. 数据预处理：归一化细胞丰度、计算空间邻居、计算免疫抑制基因评分
    4. 特征工程：计算多种评分指标和邻居特征
    5. 区域标注：将每个 spot 归类为 tumor_core / tumor_edge / stroma_immune / other
    6. 结果输出：保存 CSV 结果表、元数据文件
    7. 可视化：生成多种分析图表
    8. 特征基因：差异分析筛选 niche 特征基因，保存供 TCGA 生存分析使用

    【返回值】
    int，进程退出码：0 表示正常完成
    """
    # ---- Step 1: 初始化 ----
    _setup_logging()          # 配置日志格式和输出
    args = _parse_args()      # 解析命令行参数
    np.random.seed(args.seed) # 固定随机种子，确保结果可复现
    # 创建输出目录（如果不存在则自动创建，存在也不报错）
    args.out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.out_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 2: 数据加载与校验 ----
    logging.info("Loading spatial AnnData: %s", args.adata)
    adata = ad.read_h5ad(args.adata)  # 读取 .h5ad 格式的空间转录组数据
    # 必须包含空间坐标，否则无法进行空间分析
    if "spatial" not in adata.obsm:
        raise KeyError("adata.obsm['spatial'] is required for spatial niche analysis.")

    # 提取细胞类型丰度矩阵，并验证所需细胞类型列都存在
    abundance = _get_abundance(adata, args.abundance_key) # 提取细胞丰度矩阵，清洗列名
    required = [args.hepatocyte_col, args.treg_col, args.myeloid_col, args.fibroblast_col, args.tnk_col]
    missing = [col for col in required if col not in abundance.columns]
    if missing:
        raise KeyError(f"Missing cell abundance columns: {missing}; available: {list(abundance.columns)}")

    # ---- Step 3: 数据预处理 ----
    # 将细胞丰度归一化为比例（每个 spot 各细胞类型占比之和为1）
    # replace(0, np.nan) 避免总和为0时除以零，fillna(0.0) 将结果 NaN 填为0
    # 把 abundance 这个表按“行”做归一化，得到每个 spot 内各种细胞类型的相对比例
    proportions = abundance.div(abundance.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0) # 得到的proportions是细胞类型相对丰度矩阵(df)
    # 提取空间坐标数组，形状为 (n_spots, 2)
    coords = np.asarray(adata.obsm["spatial"])
    # 计算空间邻居列表和搜索半径
    neighbors, radius = _spatial_neighbors(coords, args.neighbor_radius_multiplier)
    logging.info("Using spatial neighbor radius %.3f.", radius) # 打印邻域半径
    # 邻域半径：定义哪些算邻居节点

    # 计算免疫抑制基因模块评分
    gene_score, marker_genes = _module_score(adata, IMMUNOSUPPRESSIVE_GENES) # 计算免疫抑制相关基因模块分数，得到每个spot的基因模块得分和参与计算的基因
    logging.info("Immunosuppressive marker genes used: %s", ", ".join(marker_genes) if marker_genes else "none")

    # ---- Step 4: 构建分析 DataFrame ----
    df = proportions.copy()
    df.insert(0, "spot_id", adata.obs_names)  # 在第0列插入 spot ID
    df["spatial_x"] = coords[:, 0]            # x 坐标
    df["spatial_y"] = coords[:, 1]            # y 坐标
    df["neighbor_radius"] = radius            # 记录使用的搜索半径，这是个标量，不同的spot值相同
    # 将基因评分按 spot 索引对齐后写入
    df["immunosuppressive_gene_score"] = gene_score.reindex(df.index).to_numpy() # 把前面算出的免疫抑制基因模块分数写入 df

    # 计算高肝细胞区域标志
    hep = df[args.hepatocyte_col] # 提取肝细胞对应的列
    hep_thr = float(hep.quantile(args.hep_high_quantile))  # 分位数阈值，利用0.75分位数值作为高肝细胞区域
    df["hep_high"] = hep >= hep_thr  # True 表示该 spot 属于高肝细胞（肿瘤核心）区域
    # 计算每个 spot 到高肝细胞区域的最短距离
    df["distance_to_hep_high"] = _distance_to_mask(coords, df["hep_high"].to_numpy()) # 计算每个 spot 到高肝细胞区域的最短空间距离

    # 计算每种细胞类型在空间邻居范围内的平均比例
    for col in [args.hepatocyte_col, "Treg", args.tnk_col, args.myeloid_col, args.fibroblast_col]: # 做邻域统计的几类细胞
        out_col = col.replace("/", "_")  # 列名中的"/"替换为"_"，避免文件命名问题
        df[f"neighbor_{out_col}"] = _neighbor_mean(df[col], neighbors)

    # 判断每个 spot 是否有至少一个高肝细胞邻居（用于标注肿瘤边缘）
    has_hep_high_neighbor = []
    for i, idx in enumerate(neighbors): # neighbors 中每个元素 idx 是某个 spot 的邻居下标集合
        has_hep_high_neighbor.append(bool(len(idx) and df["hep_high"].to_numpy()[idx].any())) # hep_high 表示是否为目标区域
    # 如果这个 spot 有邻居，并且这些邻居里至少有一个属于高肝细胞区域，就记 True，否则记 False
    df["has_hep_high_neighbor"] = has_hep_high_neighbor

    # ---- Step 5: 评分计算 ----
    # 免疫基质评分：Treg + T/NK + 髓系 + 成纤维 的 Z-score 之和
    immune_stroma = (
        _zscore(df["Treg"])
        + _zscore(df[args.tnk_col])
        + _zscore(df[args.myeloid_col])
        + _zscore(df[args.fibroblast_col])
    )
    # Z-score 的作用是把每一列都变成“以均值为 0、以标准差为尺度”的标准化值，这样不同列就能在同一个量纲上比较和求和。
    # Treg 样评分：整合 Treg 比例和免疫抑制基因表达
    df["Treg_like_score"] = _zscore(df["Treg"]) + _zscore(df["immunosuppressive_gene_score"])
    # 如果一个 spot 同时有较高的 Treg 成分，并且免疫抑制基因表达也高，那么它就更像一个“具有 Treg/免疫抑制特征”的区域
    df["immune_stroma_score"] = immune_stroma # 一个 spot 同时具有多少免疫/基质相关成分的综合强度
    # 核心评分：综合 Treg、髓系、成纤维、免疫抑制基因、邻居肝细胞比例的 Z-score 加总
    # neighbor_Hepatocyte 缺失时用自身肝细胞比例补充（孤立 spot 无邻居时的回退策略）
    # 免疫抑制 niche 核心评分
    df["immunosuppressive_niche_score"] = (
        _zscore(df["Treg"])
        + _zscore(df[args.myeloid_col])
        + _zscore(df[args.fibroblast_col])
        + _zscore(df["immunosuppressive_gene_score"])
        + _zscore(df["neighbor_Hepatocyte"].fillna(df[args.hepatocyte_col]))
    )

    # 按分位数阈值划定高 niche 评分区域
    niche_thr = float(df["immunosuppressive_niche_score"].quantile(args.niche_high_quantile))
    df["niche_high"] = df["immunosuppressive_niche_score"] >= niche_thr

    # ---- Step 6: 空间区域标注 ----
    stroma_thr = float(df["immune_stroma_score"].quantile(0.60))  # 免疫基质评分的60%分位数
    region = pd.Series("other", index=df.index)  # 所有 spot 初始化为 "other"
    # 规则1：高肝细胞 spot → 肿瘤核心
    region[df["hep_high"]] = "tumor_core"
    # 规则2：非肿瘤核心 + 高免疫基质评分 → 免疫基质区
    region[(~df["hep_high"]) & (df["immune_stroma_score"] >= stroma_thr)] = "stroma_immune"
    # 规则3：非肿瘤核心 + 有高肝细胞邻居（紧邻肿瘤）+ 高免疫基质评分 → 肿瘤边缘
    # 注意：规则3覆盖规则2，肿瘤边缘比纯免疫基质更精确
    region[(~df["hep_high"]) & df["has_hep_high_neighbor"] & (df["immune_stroma_score"] >= stroma_thr)] = "tumor_edge"
    df["spatial_region"] = region

    # ---- Step 7: 保存结果 ----
    out_csv = args.out_dir / "spatial_niche_scores.csv"
    df.to_csv(out_csv, index=False)  # 保存完整评分表（不保存行索引）
    logging.info("Niche score table saved: %s", out_csv)

    # 保存分析参数元数据，便于后续复现和追踪
    metadata = pd.DataFrame(
        [
            ("abundance_key", args.abundance_key),
            ("neighbor_radius", radius),
            ("hep_high_quantile", args.hep_high_quantile),
            ("hep_high_threshold", hep_thr),
            ("niche_high_quantile", args.niche_high_quantile),
            ("niche_high_threshold", niche_thr),
            ("n_niche_high", int(df["niche_high"].sum())),
            ("marker_genes_used", ",".join(marker_genes)),
        ],
        columns=["parameter", "value"],
    )
    metadata.to_csv(args.out_dir / "spatial_niche_parameters.csv", index=False)

    # ---- Step 8: 可视化 ----
    # 各细胞类型的空间分布图
    _spatial_scatter(df, args.hepatocyte_col, plot_dir / "spatial_hepatocyte.png", "Hepatocyte proportion")
    _spatial_scatter(df, "Treg", plot_dir / "spatial_treg.png", "Treg-like abundance")
    _spatial_scatter(df, args.myeloid_col, plot_dir / "spatial_myeloid.png", "Myeloid proportion")
    _spatial_scatter(df, args.fibroblast_col, plot_dir / "spatial_fibroblast.png", "Fibroblast proportion")
    # 免疫抑制 niche 评分的空间分布图（使用 magma 配色，黑-紫-橙-黄）
    _spatial_scatter(
        df,
        "immunosuppressive_niche_score",
        plot_dir / "spatial_immunosuppressive_niche_score.png",
        "Immunosuppressive niche score",
        cmap="magma",
    )
    # 空间区域标注图（分类变量，使用固定颜色）
    _spatial_scatter(
        df,
        "spatial_region",
        plot_dir / "spatial_region_labels.png",
        "Spatial region labels",
        categorical=True,
    )
    # 距离 vs 评分折线图
    _plot_distance(df, plot_dir / "distance_to_hep_high_vs_niche_score.png")
    # 各区域评分箱线图
    _plot_region_box(df, plot_dir / "region_score_boxplots.png")
    # 肝细胞 vs Treg 散点图
    _plot_hep_treg(df, plot_dir / "hepatocyte_vs_treg_niche_score.png")
    # 细胞类型相关性热图
    corr_cols = [
        args.hepatocyte_col,
        "Treg",
        args.tnk_col,
        args.myeloid_col,
        args.fibroblast_col,
        "Treg_like_score",
        "immunosuppressive_niche_score",
    ]
    _plot_correlation(df, corr_cols, plot_dir / "celltype_niche_correlation.png")

    # ---- Step 9: 特征基因分析 ----
    # 对 niche_high 区域做差异基因分析，筛选前 top_n 个特征基因
    ranked = _rank_niche_genes(adata, df["niche_high"], args.top_niche_genes)
    # 保存带 log2FC 的完整基因排名表
    ranked.to_csv(args.out_dir / "immunosuppressive_niche_signature_genes_ranked.csv", index=False)
    # 保存纯基因名列表（每行一个基因），供下游分析使用
    signature_path = args.out_dir / "immunosuppressive_niche_signature_genes.txt"
    signature_text = chr(10).join(ranked["gene"].tolist()) + chr(10)  # chr(10) = "\n"（换行符）
    signature_path.write_text(signature_text, encoding="utf-8")
    # 同时保存到用户指定的路径（供 TCGA 生存分析直接使用）
    args.signature_out.parent.mkdir(parents=True, exist_ok=True)
    args.signature_out.write_text(signature_text, encoding="utf-8")
    logging.info("Niche signature genes saved: %s", signature_path)
    logging.info("TCGA-ready signature also saved: %s", args.signature_out)

    # 提示用户下一步可以运行的命令
    logging.info("Done. Run DE validation with: python code/run_de_analysis.py --coloc %s --coloc-column niche_high", out_csv)
    return 0  # 返回 0 表示程序正常结束


if __name__ == "__main__":
    # 脚本直接运行时的入口
    # raise SystemExit(main()) 将 main() 的返回值（整数退出码）传递给操作系统
    # 等价于 sys.exit(main())，但无需 import sys
    # 这样 Shell 脚本可以通过 $? 检查脚本是否成功执行
    raise SystemExit(main())
