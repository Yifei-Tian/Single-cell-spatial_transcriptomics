from pathlib import Path

import pandas as pd

from colocation import add_coloc_simple


def main():
    demo_path = Path("demo_spot_cell_proportion.csv")
    out_path = Path("demo_spot_with_coloc_label.csv")

    df = pd.DataFrame(
        {
            "spot_id": ["s1", "s2", "s3"],
            "tumor": [0.4, 0.2, 0.8],
            "Treg": [0.2, 0.05, 0.12],
        }
    )
    df.to_csv(demo_path, index=False)

    df = add_coloc_simple(df, tumor_col="tumor", treg_col="Treg", tumor_threshold=0.3, treg_threshold=0.1, label_col="coloc")
    df.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()

