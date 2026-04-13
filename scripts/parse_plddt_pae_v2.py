import os, glob, json
import numpy as np
import pandas as pd


def find_outputs(cid, outdir):
    scores = sorted(glob.glob(os.path.join(outdir, f"{cid}*scores*.json")))
    pdbs   = sorted(glob.glob(os.path.join(outdir, f"{cid}*.pdb")))
    return (scores[0] if scores else None), (pdbs[0] if pdbs else None)


def res_plddt_from_af_pdb(pdb_path):
    seen = set()
    res_plddt = {}

    with open(pdb_path, "r") as f:
        for ln in f:
            if not ln.startswith("ATOM"):
                continue

            chain = ln[21].strip()
            try:
                resnum = int(ln[22:26].strip())
            except:
                continue

            icode = ln[26].strip()
            key = (chain, resnum, icode)

            if key in seen:
                continue
            seen.add(key)

            b = float(ln[60:66].strip())
            res_plddt[(chain, resnum)] = b

    return res_plddt


def median_bidirectional_pae(pae, h3_idx, t_idx):
    vals1 = pae[np.ix_(h3_idx, t_idx)].ravel()
    vals2 = pae[np.ix_(t_idx, h3_idx)].ravel()
    vals = np.concatenate([vals1, vals2])
    return float(np.median(vals))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--antigen_len", required=True, type=int)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    results = []

    for _, row in df.iterrows():
        cid = row["candidate_id"]
        scores_path, af_pdb_path = find_outputs(cid, args.outdir)

        rec = {
            "candidate_id": cid,
            "af_status": "missing"
        }

        if scores_path and af_pdb_path:
            with open(scores_path, "r") as f:
                scores = json.load(f)

            pae = np.array(scores["pae"], dtype=float)

            lenH = int(row["vh_len"])
            lenL = int(row["vl_len"])
            lenT = args.antigen_len

            offT = lenH + lenL

            # CSV already gives H3 positions directly
            h3_start = int(row["cdrh3_start"]) - 1
            h3_end = int(row["cdrh3_end"]) - 1

            h3_idx = list(range(h3_start, h3_end + 1))
            t_idx = list(range(offT, offT + lenT))

            median_pae_h3_t = median_bidirectional_pae(pae, h3_idx, t_idx)

            # pLDDT from AF PDB
            plddt_map = res_plddt_from_af_pdb(af_pdb_path)

            mean_plddt_H3 = None
            for chain_id in ["A", "H"]:
                vals = [
                    plddt_map.get((chain_id, i + 1), None)
                    for i in h3_idx
                ]
                vals = [v for v in vals if v is not None]
                if vals:
                    mean_plddt_H3 = float(np.mean(vals))
                    break

            rec.update({
                "af_status": "ok",
                "mean_plddt_H3": mean_plddt_H3,
                "median_pae_H3_T_bidirectional": median_pae_h3_t
            })

        results.append(rec)

    out = df.merge(pd.DataFrame(results), on="candidate_id", how="left")
    out.to_csv(args.out_csv, index=False)

    print("Wrote:", args.out_csv)
    print(out[[
        "candidate_id",
        "mean_plddt_H3",
        "median_pae_H3_T_bidirectional"
    ]].head())


if __name__ == "__main__":
    main()