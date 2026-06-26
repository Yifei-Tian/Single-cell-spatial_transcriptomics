"""
================================================================================
脚本名称: run_de_analysis.py
功能概述: 空间转录组差异表达（DE）分析 —— Step 3
================================================================================

【整体任务说明】
    本脚本将空间转录组数据（spot 级别的基因表达矩阵）与分组标签（共定位标签或
    niche_high 标签）结合，使用 Wilcoxon 秩和检验找出"阳性 spot"相比于其他 spot
    中特异性高表达的基因，最终导出 signature genes（空间特征基因）列表。

    执行步骤：
      步骤 1  加载空间 AnnData（adata_vis_post.h5ad），记录基因数和 spot 数；
      步骤 2  加载共定位/分组标签 CSV，按 spot ID 对齐注入 adata.obs；
              缺失 spot 的标签自动填充为阴性（0），确保 adata 完整性；
      步骤 3  准备表达矩阵：
                - 若指定 layer（默认 "log1p"）存在，直接使用；
                - 否则对当前 adata.X 执行 normalize_total（目标 10000）+ log1p 兜底；
      步骤 4  运行 sc.tl.rank_genes_groups（Wilcoxon 秩和检验），
              对阳性组 vs 其余 spot 进行组间差异分析；
      步骤 5  自动识别"阳性"类别（通常为 1 或 "1"），兼容 scanpy 多版本的结果格式，
              提取该组 Top-N 差异基因写入输出文件。

【输入文件】
    results/adata_vis_post.h5ad       - 经过预处理的空间转录组 AnnData（含基因表达矩阵）
    results/spot_with_coloc_label.csv - 共定位标签文件（含 spot_id 列 + coloc 标签列）
                                        也可传入 spatial_niche_scores.csv 并用 --coloc-column niche_high

【输出文件】
    results/spatial_signature_genes.txt - 阳性 spot 特异高表达基因列表
                                          （每行一个基因名，默认保留前 50 个）

【整体数据流】
    adata_vis_post.h5ad          spot_with_coloc_label.csv（或 spatial_niche_scores.csv）
      （空间基因表达矩阵）                  （每个 spot 的分组标签）
            │                                        │
            └──────────────┬──────────────────────────┘
                           ▼
              注入标签到 adata.obs（缺失 spot 填充阴性）
                           ▼
              准备 log1p 归一化表达矩阵（layer 或兜底归一化）
                           ▼
              Wilcoxon 检验：阳性 spot vs 其余 spot
                           ▼
              提取阳性组 Top-N 差异基因（默认 50 个）
                           ▼
                 spatial_signature_genes.txt

【主要命令行参数】
    --adata          空间 AnnData 文件路径（默认 results/adata_vis_post.h5ad）
    --coloc          分组标签 CSV 文件路径（默认 results/spot_with_coloc_label.csv）
    --out            输出基因列表路径（默认 results/spatial_signature_genes.txt）
    --top-n          导出前 N 个差异基因（默认 50）
    --coloc-column   CSV 中分组标签列名（默认 "coloc"，也可指定 "niche_high"）
    --spot-id-column CSV 中 spot ID 列名（默认 "spot_id"）
    --layer          使用的 AnnData layer 名称（默认 "log1p"，不存在时自动归一化）

【设计亮点】
    大量防御性编程：缺失标签填充、layer 兜底归一化、group 名称兼容多版本 scanpy，
    使整个流程在数据不完整或上游格式变化时也能鲁棒地运行。

【依赖关系】
    上游：run_preprocessing.py（生成 adata_vis_post.h5ad）
          colocation.py（生成 spot_with_coloc_label.csv）
          或 run_spatial_niche_analysis.py（生成 spatial_niche_scores.csv）
    下游：tcga_survival_analysis.R（读取 spatial_signature_genes.txt 进行 ssGSEA 预后分析）
================================================================================
"""

from __future__ import annotations

# ── 标准库导入 ──────────────────────────────────────────────────────────────────
from pathlib import Path   # 用于跨平台的文件路径处理
import argparse             # 用于解析命令行参数
import logging              # 用于输出运行日志（替代 print，更专业）
import sys                  # 用于访问标准输出流和程序退出码

# ── 第三方库导入 ────────────────────────────────────────────────────────────────
import pandas as pd         # 用于读取和处理 CSV 表格数据
import scanpy as sc         # 单细胞/空间转录组分析核心库，处理 AnnData 对象


# ══════════════════════════════════════════════════════════════════════════════
# 函数一：_parse_args()
# 任务：解析命令行参数，让脚本可以灵活地被外部调用，不需要每次修改源码
# ══════════════════════════════════════════════════════════════════════════════
def _parse_args() -> argparse.Namespace:
    """
    解析命令行参数，返回包含所有参数值的命名空间对象。

    支持的参数：
        --adata         输入的 AnnData 文件路径（.h5ad 格式，包含基因表达矩阵）
        --coloc         输入的共定位标签 CSV 文件路径
        --out           输出的 signature gene 列表文件路径（.txt）
        --top-n         导出前 N 个差异最显著的基因，默认 50
        --coloc-column  CSV 中代表共定位结果的列名，默认 "coloc"
        --spot-id-column CSV 中 spot ID 的列名，默认 "spot_id"
        --layer         使用 AnnData 中哪个数据层（layer）做差异分析，默认 "log1p"

    实现方式：
        使用 Python 内置的 argparse 库，每个 add_argument() 调用注册一个参数。
        default_root 通过 __file__（当前脚本路径）向上推两级找到项目根目录，
        使默认路径不依赖于执行时的当前工作目录。
    """
    parser = argparse.ArgumentParser(description="Run Step 3 DE analysis on spatial data.")

    # 通过当前脚本文件的路径，向上推两级找到项目根目录
    # 例如脚本在 /project/scripts/run_de_analysis.py，则 default_root = /project
    default_root = Path(__file__).resolve().parents[1]
    default_results = default_root / "results"   # 默认结果目录：项目根目录/results

    # 参数1：输入的 AnnData 数据文件（.h5ad 是单细胞数据的标准格式）
    parser.add_argument(
        "--adata",
        type=Path,
        default=default_results / "adata_vis_post.h5ad",
        help="Path to adata_vis_post.h5ad",
    )

    # 参数2：输入的共定位标签 CSV 文件（记录每个 spot 是否处于共定位区域）
    parser.add_argument(
        "--coloc",
        type=Path,
        default=default_results / "spot_with_coloc_label.csv",
        help="Path to spot_with_coloc_label.csv",
    )

    # 参数3：输出文件路径，保存筛选出的 signature gene 列表
    parser.add_argument(
        "--out",
        type=Path,
        default=default_results / "spatial_signature_genes.txt",
        help="Output file for signature genes",
    )

    # 参数4：保留多少个差异基因，默认取前 50 个
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top genes to export",
    )

    # 参数5：CSV 文件中存储共定位标签（0/1）的列名
    parser.add_argument(
        "--coloc-column",
        type=str,
        default="coloc",
        help="Column name in CSV to use as coloc labels",
    )

    # 参数6：CSV 文件中存储 spot ID 的列名（用来和 AnnData 对齐）
    parser.add_argument(
        "--spot-id-column",
        type=str,
        default="spot_id",
        help="Column name in CSV to use as spot IDs",
    )

    # 参数7：使用 AnnData 中哪个 layer（数据层）的表达量做差异分析
    # AnnData 可以储存多个版本的表达矩阵（原始计数、归一化、log 变换等），
    # 默认使用预处理阶段保存的 "log1p" 层（已归一化 + log 变换）
    parser.add_argument(
        "--layer",
        type=str,
        default="log1p",
        help="Layer to use for DE (falls back to X if missing)",
    )

    return parser.parse_args()   # 解析并返回所有参数


# ══════════════════════════════════════════════════════════════════════════════
# 函数二：_setup_logging()
# 任务：配置日志系统，让程序运行过程中的关键信息以统一格式输出到屏幕
# ══════════════════════════════════════════════════════════════════════════════
def _setup_logging() -> None:
    """
    初始化 Python 标准日志系统。

    输出格式：  2026-06-06 12:00:00,000 | INFO | 正在加载数据...
                （时间戳）              （级别） （消息内容）

    实现方式：
        使用 logging.basicConfig() 一次性配置：
        - level=INFO：只输出 INFO 及以上级别的消息（DEBUG 级别的会被过滤掉）
        - format：指定每条日志的格式字符串
        - handlers：输出到 sys.stdout（标准输出），而非默认的 stderr，
          方便在 shell 脚本中用 > 重定向日志到文件
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


# ══════════════════════════════════════════════════════════════════════════════
# 函数三：_select_target_group()
# 任务：从共定位标签的类别列表中，自动识别"阳性"类别（通常是 1）
# ══════════════════════════════════════════════════════════════════════════════
def _select_target_group(categories: list) -> object:
    """
    自动识别共定位标签中代表"阳性"的类别值，作为差异分析的 target group。

    背景：
        共定位标签由上游流程产生，不同版本可能用不同类型表示"真/阳性"：
        有的是整数 1，有的是字符串 "1"，有的是布尔值 True，有的是 "1.0"。
        本函数统一处理这些情况，避免因类型不匹配导致找不到目标组。

    参数：
        categories (list): adata.obs 中共定位标签列的所有类别值列表

    返回：
        找到的"阳性"类别值；如果都没有匹配，则返回 categories 中的最后一个

    实现方式：
        按优先级依次检查 "1"、1、"1.0"、True、"True" 是否在 categories 中，
        第一个匹配到的就返回。如果全部未命中，降级使用最后一个类别（兜底策略）。
    """
    # 依次尝试各种"真值"的常见表示形式，兼容不同上游产生的标签格式
    for cand in ("1", 1, "1.0", True, "True"):
        if cand in categories:
            return cand
    # 兜底：如果没有找到任何已知"阳性"标识，使用最后一个类别
    return categories[-1]


# ══════════════════════════════════════════════════════════════════════════════
# 函数四：main()
# 任务：脚本的主执行函数，按顺序完成数据加载→标签融合→差异分析→结果导出
# 返回值：整数退出码（0 = 成功）
# ══════════════════════════════════════════════════════════════════════════════
def main() -> int:
    """
    主流程函数，依次执行以下 5 个步骤：
        1. 加载 AnnData 和共定位标签数据
        2. 将共定位标签对齐并注入 AnnData
        3. 准备用于差异分析的表达矩阵
        4. 运行 Wilcoxon 差异表达检验
        5. 提取阳性组的 Top-N 基因并写入文件
    """
    # 初始化日志系统（在最早时刻调用，确保后续所有日志都能正常输出）
    _setup_logging()
    # 解析命令行参数（获取所有输入/输出路径和配置选项）
    args = _parse_args()

    # ──────────────────────────────────────────────────────────────────────────
    # 步骤 1：加载数据
    # ──────────────────────────────────────────────────────────────────────────
    logging.info("Loading adata: %s", args.adata)
    # sc.read_h5ad() 读取 .h5ad 文件，这是 AnnData 对象的磁盘存储格式。
    # AnnData 是一个"数据框+矩阵"的组合体：
    #   - adata.X          : 基因表达矩阵（行 = spot，列 = 基因）
    #   - adata.obs        : 每个 spot 的元数据（类似 DataFrame）
    #   - adata.var        : 每个基因的元数据
    #   - adata.layers     : 存储多个版本的表达矩阵（如原始计数、归一化、log变换）
    #   - adata.uns        : 非结构化数据（如差异分析结果）
    adata = sc.read_h5ad(args.adata)

    logging.info("Loading coloc labels: %s", args.coloc)
    # 读取共定位标签 CSV，该文件由上一步（空间 niche 分析）产生，
    # 通常包含两列：spot_id（spot 标识符）和 coloc（0=阴性，1=阳性）
    df = pd.read_csv(args.coloc)

    # ──────────────────────────────────────────────────────────────────────────
    # 步骤 2：将共定位标签对齐并注入 AnnData
    # ──────────────────────────────────────────────────────────────────────────

    # 如果 CSV 中有 spot_id 列，将其设为行索引（行标签），
    # 这样后续可以通过 spot ID 直接索引，方便与 AnnData 对齐
    if args.spot_id_column in df.columns:
        df = df.set_index(args.spot_id_column)

    # 安全校验：如果 coloc 列不存在，立即报错，避免后续产生难以追踪的错误
    if args.coloc_column not in df.columns:
        raise ValueError(f"Missing coloc column: {args.coloc_column}")

    # 提取共定位标签列（一个 Series，index 是 spot_id，值是 0/1）
    labels = df[args.coloc_column]

    # 关键操作：reindex() 按照 AnnData 中 spot 的顺序重新排列标签。
    # 这解决了"CSV 和 AnnData 中 spot 顺序可能不一致"的问题。
    # 如果 CSV 中有某个 spot 在 AnnData 里找不到，对应位置会变成 NaN。
    labels = labels.reindex(adata.obs_names)

    # 防御性处理：如果有 spot 没有匹配到标签（NaN），用 0（阴性）填充，
    # 并给出警告，让用户知道数据有缺失，而不是静默地忽略
    if labels.isna().any():
        missing = int(labels.isna().sum())
        logging.warning("%d spots missing coloc labels; filling with 0.", missing)
        labels = labels.fillna(0)

    # 将处理好的标签写入 adata.obs（spot 的元数据表），列名为 "coloc"。
    # 必须转为 category 类型（分类变量），scanpy 的分组分析函数要求这种格式。
    adata.obs[args.coloc_column] = labels.astype("category")

    # ──────────────────────────────────────────────────────────────────────────
    # 步骤 3：准备用于差异分析的表达矩阵（确保数据经过归一化和 log 变换）
    # ──────────────────────────────────────────────────────────────────────────

    if args.layer in adata.layers:
        # 优先路径：直接使用预处理阶段保存的 "log1p" layer。
        # .copy() 是为了防止修改 X 时意外修改 layers 中的原始数据（浅拷贝陷阱）。
        adata.X = adata.layers[args.layer].copy()
    else:
        # 兜底路径：如果指定的 layer 不存在，对当前 X 重新做归一化和 log 变换。
        logging.warning(
            "Layer '%s' missing; normalizing current X with normalize_total + log1p before DE.",
            args.layer,
        )
        # normalize_total：将每个 spot 的总表达量标准化到 target_sum（默认 10,000）。
        # 即每个 spot 的每个基因表达量 = (原始计数 / 该 spot 总计数) × 10000。
        # 目的：消除不同 spot 测序深度（捕获的 RNA 总量）不同带来的偏差。
        sc.pp.normalize_total(adata, target_sum=1e4)

        # log1p：对每个值做 log(x + 1) 变换（+1 是为了避免 log(0) 报错）。
        # 目的：压缩数据范围，使高表达基因不会过于主导统计结果，
        #       同时让数据分布更接近正态，符合统计检验的假设。
        sc.pp.log1p(adata)

    # ──────────────────────────────────────────────────────────────────────────
    # 步骤 4：运行 Wilcoxon 差异表达检验（核心计算步骤）
    # ──────────────────────────────────────────────────────────────────────────

    # sc.tl.rank_genes_groups() 是 scanpy 的标准差异表达函数。
    # 工作原理：
    #   对 coloc 列的每个类别（如 0 和 1），分别找出该类别中
    #   相较于其他所有类别特异性高表达的基因，并按统计显著性排序。
    # 参数说明：
    #   groupby="coloc"      - 按 coloc 列的值（0/1）对 spot 进行分组
    #   method="wilcoxon"    - 使用 Wilcoxon 秩和检验（非参数检验，
    #                          不假设数据服从正态分布，对单细胞数据更合适）
    # 结果存储在 adata.uns["rank_genes_groups"] 中，包含：
    #   - names        : 每个组按排名排列的基因名
    #   - scores       : 统计得分（越高越显著）
    #   - pvals        : p 值
    #   - logfoldchanges: log2 倍数变化（衡量表达量差异大小）
    sc.tl.rank_genes_groups(adata, groupby=args.coloc_column, method="wilcoxon")

    # ──────────────────────────────────────────────────────────────────────────
    # 步骤 5：提取阳性组的 Top-N 差异基因并写入文件
    # ──────────────────────────────────────────────────────────────────────────

    # 获取共定位标签的所有类别（如 [0, 1] 或 ["0", "1"]）
    categories = list(adata.obs[args.coloc_column].cat.categories)

    # 如果只有 1 个类别（所有 spot 标签相同），无法进行组间比较，跳过
    if len(categories) < 2:
        logging.warning("coloc categories < 2; skipping DE signature export.")
        return 0

    # 自动识别"阳性"类别（通常是 1），作为我们关注的 target group
    target_group = _select_target_group(categories)

    # 从 adata.uns 中取出差异分析结果的基因名矩阵
    names = adata.uns.get("rank_genes_groups", {}).get("names")
    if names is None:
        logging.warning("rank_genes_groups missing; cannot export signature.")
        return 0

    # 兼容性处理：scanpy 不同版本返回的 names 结构不同，
    # 可能是 numpy structured array（通过 dtype.names 获取列名），
    # 也可能是 pandas DataFrame（通过 .columns 获取列名）。
    # 这里统一提取可用的组名列表。
    if hasattr(names, "dtype") and getattr(names.dtype, "names", None):
        # numpy structured array 格式：列名存在 dtype.names 中
        available_groups = list(names.dtype.names)
    else:
        # pandas DataFrame 格式：列名直接用 .columns 获取
        available_groups = list(getattr(names, "columns", []))

    # 如果没有找到任何组，说明差异分析结果有问题，跳过
    if not available_groups:
        logging.warning("No DE groups found; cannot export signature.")
        return 0

    # 将 target_group 转为字符串进行比较（统一类型，避免 1 != "1" 的问题）
    if str(target_group) not in available_groups:
        # 如果目标组不在结果中（可能因类型格式不匹配），降级使用最后一个组
        logging.warning("Target group '%s' not in DE groups %s; falling back to '%s'.", target_group, available_groups, available_groups[-1])
        target_group = available_groups[-1]
    else:
        target_group = str(target_group)

    # 用 scanpy 的便捷函数将差异分析结果转成 DataFrame，
    # 结果已按统计得分降序排列（最显著的基因排在最前面）。
    # 列包括：names（基因名）、scores、pvals、pvals_adj、logfoldchanges
    df_de = sc.get.rank_genes_groups_df(adata, group=target_group)

    # 取前 top_n 个基因名，转为 Python 列表
    genes = df_de["names"].head(args.top_n).tolist()

    # 确保输出目录存在（parents=True：同时创建所有缺失的父目录；
    # exist_ok=True：目录已存在时不报错）
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # 将基因列表写入文本文件，每行一个基因名，末尾加换行符
    args.out.write_text("\n".join(genes) + "\n", encoding="utf-8")
    logging.info("Signature genes saved: %s", args.out)

    return 0   # 返回 0 表示程序正常结束


# ══════════════════════════════════════════════════════════════════════════════
# 脚本入口
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    # 当直接用 python run_de_analysis.py 运行时，执行 main() 函数。
    # raise SystemExit(main()) 将 main() 的返回值（整数退出码）传递给操作系统，
    # 使 shell 脚本或 CI/CD 系统可以通过 $? 检查程序是否成功完成。
    # （0 = 成功，非 0 = 失败，这是 Unix 程序的标准约定）
    raise SystemExit(main())
