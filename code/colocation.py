import argparse
from pathlib import Path

import pandas as pd


def add_coloc_simple(df, tumor_col, treg_col, tumor_threshold, treg_threshold, label_col):
    """Simple rule-based co-location label."""
    df[label_col] = ((df[tumor_col] > tumor_threshold) & (df[treg_col] > treg_threshold)).astype(int)
    return df


def add_coloc_logistic(df, tumor_col, treg_col, sum_threshold, score_col, label_col):
    """Logistic-regression-based co-location score and label."""
    from sklearn.linear_model import LogisticRegression

    x = df[[tumor_col, treg_col]]
    y = (df[tumor_col] + df[treg_col] > sum_threshold).astype(int)
    model = LogisticRegression().fit(x, y)
    df[score_col] = model.predict_proba(x)[:, 1]
    df[label_col] = (df[score_col] > 0.5).astype(int)
    return df


def parse_args():
    parser = argparse.ArgumentParser(description="Define spatial co-location regions from spot cell proportions.")
    parser.add_argument("--input", required=True, help="Path to spot_cell_proportion.csv")
    parser.add_argument("--output", required=True, help="Path to spot_with_coloc_label.csv")
    parser.add_argument("--method", choices=["simple", "logistic"], default="simple")
    parser.add_argument("--tumor-col", default="tumor")
    parser.add_argument("--treg-col", default="Treg")
    parser.add_argument("--tumor-threshold", type=float, default=0.3)
    parser.add_argument("--treg-threshold", type=float, default=0.1)
    parser.add_argument("--sum-threshold", type=float, default=0.4)
    parser.add_argument("--label-col", default="coloc")
    parser.add_argument("--score-col", default="coloc_score")
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    df = pd.read_csv(input_path)
    if args.tumor_col not in df.columns or args.treg_col not in df.columns:
        raise KeyError(f"Missing columns: {args.tumor_col} or {args.treg_col}")

    if args.method == "simple":
        df = add_coloc_simple(
            df,
            tumor_col=args.tumor_col,
            treg_col=args.treg_col,
            tumor_threshold=args.tumor_threshold,
            treg_threshold=args.treg_threshold,
            label_col=args.label_col,
        )
    else:
        df = add_coloc_logistic(
            df,
            tumor_col=args.tumor_col,
            treg_col=args.treg_col,
            sum_threshold=args.sum_threshold,
            score_col=args.score_col,
            label_col=args.label_col,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()

