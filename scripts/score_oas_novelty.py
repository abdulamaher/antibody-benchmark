#!/usr/bin/env python3
import sys, sqlite3, pandas as pd

def main():
    if len(sys.argv) != 4:
        print("Usage: score_oas_novelty.py <h3_csv> <oas_sqlite> <out_csv>", file=sys.stderr)
        sys.exit(2)

    h3_csv, db_path, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

    df = pd.read_csv(h3_csv)
    # try common column names
    h3_col = None
    for c in ["H3_seq", "h3_seq", "cdrh3", "cdr_h3", "CDRH3"]:
        if c in df.columns:
            h3_col = c
            break
    if h3_col is None:
        raise ValueError(f"Could not find H3 column. Columns={list(df.columns)}")

    con = sqlite3.connect(db_path)
    cur = con.cursor()

    def is_known(seq: str) -> int:
        if not isinstance(seq, str) or not seq.strip():
            return 0
        cur.execute("SELECT 1 FROM cdr3 WHERE cdr3_aa = ? LIMIT 1;", (seq.strip(),))
        return 1 if cur.fetchone() else 0

    known = []
    for s in df[h3_col].astype(str).tolist():
        known.append(is_known(s))

    df["oas_known_cdrh3"] = known
    df["novelty_score"] = [0 if k == 1 else 1 for k in known]  # disqualify if known

    con.close()
    df.to_csv(out_csv, index=False)

    print("Wrote:", out_csv)
    print(df["novelty_score"].value_counts(dropna=False))

if __name__ == "__main__":
    main()