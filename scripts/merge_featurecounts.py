import sys
import pandas as pd

# Merge fixcol2 featureCounts tables (e.g. PE + SE runs) on Geneid.
ANNOT_COLS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
EXTRA_COL = "uniq_column"


def read_fct(path):
    return pd.read_csv(path, sep="\t", dtype=str).set_index("Geneid", drop=False)


def merge(paths):
    tables = [read_fct(p) for p in paths]
    base = tables[0]
    parts = [base[ANNOT_COLS]]
    for path, tbl in zip(paths, tables):
        if set(tbl.index) != set(base.index):
            sys.exit(f"Gene set in {path} differs from {paths[0]}")
        sample_cols = [c for c in tbl.columns if c not in ANNOT_COLS and c != EXTRA_COL]
        parts.append(tbl.loc[base.index, sample_cols])
    parts.append(base[[EXTRA_COL]])
    return pd.concat(parts, axis=1)


if __name__ == "__main__":
    out_path = sys.argv[1]
    merge(sys.argv[2:]).to_csv(out_path, sep="\t", index=False)
