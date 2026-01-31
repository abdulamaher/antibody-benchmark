#!/usr/bin/env python3
import os, sys, tarfile, glob, gzip, sqlite3, csv, subprocess

ZENODO_TINY_URL = "https://zenodo.org/record/7502634/files/OAS-aligned-tiny.tar"  # from KA-Search README table

def run(cmd):
    print("+", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def download_if_missing(url, out_path):
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        print(f"[ok] already downloaded: {out_path}")
        return
    ensure_dir(os.path.dirname(out_path))
    run(["wget", "-O", out_path, url])

def extract_tar(tar_path, out_dir):
    marker = os.path.join(out_dir, ".extracted")
    if os.path.exists(marker):
        print(f"[ok] already extracted: {out_dir}")
        return
    ensure_dir(out_dir)
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(out_dir)
    open(marker, "w").write("ok\n")
    print(f"[ok] extracted to: {out_dir}")

def iter_rows_tsv_or_csv(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as f:
        # OAS-aligned datasets often use TSV, but we detect delimiter.
        sample = f.read(4096)
        f.seek(0)
        delim = "\t" if ("\t" in sample and sample.count("\t") > sample.count(",")) else ","
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            yield row

def find_data_tables(root_dir):
    # Common patterns: *.tsv.gz, *.csv.gz, *.tsv, *.csv
    pats = ["**/*.tsv.gz", "**/*.csv.gz", "**/*.tsv", "**/*.csv"]
    out = []
    for pat in pats:
        out.extend(glob.glob(os.path.join(root_dir, pat), recursive=True))
    # Filter out tiny files (like metadata) and keep likely big tables
    out = [p for p in out if os.path.getsize(p) > 1024 * 1024]  # >1MB
    return sorted(out)

def build_sqlite(db_path, tables):
    ensure_dir(os.path.dirname(db_path))
    if os.path.exists(db_path):
        print(f"[warn] db exists, will reuse: {db_path}")
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("PRAGMA journal_mode=WAL;")
    cur.execute("PRAGMA synchronous=NORMAL;")
    cur.execute("CREATE TABLE IF NOT EXISTS cdr3 (cdr3_aa TEXT PRIMARY KEY);")
    con.commit()

    inserted = 0
    batch = []
    BATCH_N = 10000

    for t in tables:
        print(f"[info] scanning: {t}")
        for row in iter_rows_tsv_or_csv(t):
            # AIRR rearrangement schema uses 'cdr3_aa'
            cdr3 = (row.get("cdr3_aa") or "").strip()
            if not cdr3:
                continue
            batch.append((cdr3,))
            if len(batch) >= BATCH_N:
                cur.executemany("INSERT OR IGNORE INTO cdr3(cdr3_aa) VALUES (?)", batch)
                inserted += cur.rowcount
                con.commit()
                batch = []
        if batch:
            cur.executemany("INSERT OR IGNORE INTO cdr3(cdr3_aa) VALUES (?)", batch)
            inserted += cur.rowcount
            con.commit()
            batch = []

    cur.execute("CREATE INDEX IF NOT EXISTS idx_cdr3 ON cdr3(cdr3_aa);")
    con.commit()
    cur.execute("SELECT COUNT(*) FROM cdr3;")
    n = cur.fetchone()[0]
    con.close()
    print(f"[done] db={db_path} unique_cdr3={n} inserted_recently~{inserted}")

def main():
    if len(sys.argv) != 3:
        print("Usage: build_oas_cdr3_sqlite.py <work_dir> <out_sqlite_path>", file=sys.stderr)
        print("Example: build_oas_cdr3_sqlite.py oas_cache oas_cache/oas_tiny_cdr3.sqlite", file=sys.stderr)
        sys.exit(2)

    work_dir = sys.argv[1]
    db_path = sys.argv[2]
    ensure_dir(work_dir)

    tar_path = os.path.join(work_dir, "OAS-aligned-tiny.tar")
    extract_dir = os.path.join(work_dir, "OAS-aligned-tiny")

    download_if_missing(ZENODO_TINY_URL, tar_path)
    extract_tar(tar_path, extract_dir)

    tables = find_data_tables(extract_dir)
    if not tables:
        print("[error] Could not find data tables under extracted dir. Print tree and adjust find patterns.", file=sys.stderr)
        sys.exit(1)

    print("[info] candidate tables:")
    for p in tables[:10]:
        print("  ", p)
    if len(tables) > 10:
        print(f"  ... ({len(tables)} total)")

    build_sqlite(db_path, tables)

if __name__ == "__main__":
    main()
