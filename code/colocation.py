"""
================================================================================
脚本名称: colocation.py
功能概述: 空间共定位（Co-localization）区域定义 —— 分析流程辅助模块
================================================================================

【整体任务说明】
    本脚本基于 Cell2location 反卷积输出的 spot 级别细胞比例数据，使用逻辑回归
    模型识别肿瘤细胞（Hepatocyte）与 Treg 细胞共定位的空间区域，为后续差异
    表达分析提供阳性/阴性标签。

    核心策略：
      - 以肿瘤细胞比例与 Treg 细胞比例之和作为代理变量，通过分位数阈值或
        固定阈值构建二值化训练目标；
      - 训练逻辑回归模型学习两种细胞比例与共定位概率的关系；
      - 输出每个 spot 的共定位概率评分（coloc_score）及二值化标签（coloc）。

【输入文件】
    spot_cell_proportion.csv  - Cell2location 反卷积后每个 spot 的细胞类型比例表
                                （由 run_preprocessing.py 生成）
                                必须包含以下列：
                                  - 肿瘤细胞列（默认 "Hepatocyte"）
                                  - Treg 细胞列（默认 "Treg"）

【输出文件】
    spot_with_coloc_label.csv - 在输入文件基础上新增以下两列：
                                  - coloc_score：逻辑回归模型预测的共定位概率（0~1）
                                  - coloc     ：二值化共定位标签（0=阴性，1=阳性）

【调用方式】
    作为独立脚本运行：
      python colocation.py [--input PATH] [--output PATH] [选项...]
    或被其他脚本导入调用 add_coloc_logistic() 函数。

【依赖关系】
    上游：run_preprocessing.py（生成 spot_cell_proportion.csv）
    下游：run_de_analysis.py   （读取 spot_with_coloc_label.csv 进行 DE 分析）
================================================================================
"""
import argparse
from pathlib import Path

import pandas as pd


def add_coloc_logistic(df, tumor_col, treg_col, sum_threshold, score_col, label_col, threshold_quantile=None, label_threshold=None, label_quantile=0.8):
    """
    使用逻辑回归模型定义共定位区域，并为每个 spot 计算一个“共定位分数”。

    该方法首先基于两种细胞比例之和创建一个临时目标，然后训练一个逻辑回归模型
    来学习这两种细胞比例与共定位可能性的关系，最后用模型预测每个 spot 的共定位分数。

    Args:
        df (pd.DataFrame): 包含 spot 细胞比例的 DataFrame。
        tumor_col (str): 代表肿瘤细胞比例的列名。
        treg_col (str): 代表 Treg 细胞比例的列名。
        sum_threshold (float | None): 用于创建临时训练目标的两种细胞比例之和的阈值。
        score_col (str): 输出的共定位分数（概率）列的名称。
        label_col (str): 输出的最终共定位标签列的名称。
        threshold_quantile (float | None): 当 sum_threshold 为 None 时，用该分位数自动估计阈值。
        label_threshold (float | None): 直接使用该阈值对 score 打标签；若为 None，则用分位数。
        label_quantile (float): 使用 score 的该分位数作为标签阈值。

    Returns:
        pd.DataFrame: 增加了共定位分数和标签列的 DataFrame。
    """
    # 动态导入，仅在使用此函数时才需要 sklearn
    from sklearn.linear_model import LogisticRegression

    # 准备训练数据：X 是特征，y 是临时目标
    x = df[[tumor_col, treg_col]]
    sum_series = df[tumor_col] + df[treg_col]
    if sum_threshold is None:
        quantile_value = 0.8 if threshold_quantile is None else threshold_quantile
        sum_threshold = float(sum_series.quantile(quantile_value))
        print(f"使用分位数 {quantile_value} 自动估计 sum_threshold={sum_threshold:.6f}")
    # 当两种细胞比例之和大于阈值时，我们假设它是一个“正样本” (1)
    y = (sum_series > sum_threshold).astype(int)

    # 训练逻辑回归模型
    model = LogisticRegression().fit(x, y)

    # 使用模型预测每个 spot 属于正样本（共定位）的概率，并存入新列
    # model.predict_proba(x) 返回一个 (n_samples, n_classes) 的数组，[:, 1] 表示取属于类别 "1" 的概率
    df[score_col] = model.predict_proba(x)[:, 1]

    # 基于预测的概率/分数生成标签
    if label_threshold is None:
        label_threshold = float(df[score_col].quantile(label_quantile))
        print(f"使用分位数 {label_quantile} 自动估计 label_threshold={label_threshold:.6f}")
    else:
        print(f"使用固定 label_threshold={label_threshold:.6f}")
    df[label_col] = (df[score_col] >= label_threshold).astype(int)
    return df


def parse_args():
    """
    解析和管理从命令行传入的参数。

    通过 argparse 库，我们可以方便地定义脚本需要哪些输入、它们的类型、默认值，
    并自动生成帮助信息。

    Returns:
        argparse.Namespace: 一个包含所有已解析参数的对象。
    """
    parser = argparse.ArgumentParser(description="根据 spot 的细胞比例定义空间共定位区域（仅使用逻辑回归方法）。")

    # --- 输入/输出文件 ---
    parser.add_argument("--input", default="/home/gyw/R/Python/Single_cell_train/results/spot_cell_proportion.csv", help="输入的 spot_cell_proportion.csv 文件路径。")
    parser.add_argument("--output", default="/home/gyw/R/Python/Single_cell_train/results/spot_with_coloc_label.csv", help="输出的 spot_with_coloc_label.csv 文件路径。")

    # --- 相关列 ---
    parser.add_argument("--tumor-col", default="Hepatocyte", help="代表肿瘤细胞比例的列名。")
    parser.add_argument("--treg-col", default="Treg", help="代表 Treg 细胞比例的列名。")

    # --- 逻辑回归方法参数 ---
    parser.add_argument("--sum-threshold", type=float, default=None, help="用于创建临时训练目标的细胞比例之和的阈值。若不提供，将使用分位数自动估计。")
    parser.add_argument("--threshold-quantile", type=float, default=0.8, help="当 sum_threshold 未提供时，用于估计阈值的分位数（0~1）。")
    parser.add_argument("--label-threshold", type=float, default=None, help="用于生成 coloc 标签的分数阈值。若不提供，将使用分位数自动估计。")
    parser.add_argument("--label-quantile", type=float, default=0.8, help="当 label_threshold 未提供时，用于估计标签阈值的分位数（0~1）。")

    # --- 输出列名 ---
    parser.add_argument("--label-col", default="coloc", help="最终输出的共定位标签列的名称。")
    parser.add_argument("--score-col", default="coloc_score", help="输出的共定位分数（概率）列的名称。")

    return parser.parse_args()


def main():
    """
    脚本的主执行流程。

    负责调度参数解析、文件读取、方法选择、数据处理和结果保存。
    """
    # 1. 获取所有命令行参数
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    # 2. 读取输入的细胞比例数据
    df = pd.read_csv(input_path)

    # 3. 校验关键列是否存在，避免因列名错误导致程序崩溃
    if args.tumor_col not in df.columns or args.treg_col not in df.columns:
        raise KeyError(f"输入文件中缺少指定的列: {args.tumor_col} 或 {args.treg_col}")

    # 4. 使用逻辑回归方法定义共定位区域
    print("使用 'logistic' 方法定义共定位区域...")
    df = add_coloc_logistic(
        df,
        tumor_col=args.tumor_col,
        treg_col=args.treg_col,
        sum_threshold=args.sum_threshold,
        score_col=args.score_col,
        label_col=args.label_col,
        threshold_quantile=args.threshold_quantile,
        label_threshold=args.label_threshold,
        label_quantile=args.label_quantile,
    )

    # 4.1 输出 coloc=1 的比例
    coloc_ratio = float(df[args.label_col].mean()) if len(df) else 0.0
    print(f"coloc=1 的比例: {coloc_ratio:.4f}")

    # 5. 保存处理后的结果
    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"处理完成，结果已保存至: {output_path}")


if __name__ == "__main__":
    # 当该脚本被直接执行时，运行 main 函数
    main()

