import os, glob, json, re
import numpy as np
import pandas as pd

# Detect H3 residue numbers from PDB REMARKS (your RFdiffusion/PMPNN PDBs)
H3_REMARK_RE = re.compile(r"^REMARK\s+PDBinfo-LABEL:\s+(\d+)\s+H3\s*$")

def parse_h3_resnums_from_remarks(pdb_path):
    res=[]
    if not pdb_path or not isinstance(pdb_path,str) or not os.path.exists(pdb_path):
        return res
    with open(pdb_path,"r") as f:
        for ln in f:
            m=H3_REMARK_RE.match(ln.rstrip("\n"))
            if m:
                res.append(int(m.group(1)))
    return sorted(set(res))

# Extract pLDDT from AF PDB (stored in B-factor column)
def res_plddt_from_af_pdb(pdb_path):
    seen=set()
    res_plddt={}
    with open(pdb_path,"r") as f:
        for ln in f:
            if not ln.startswith("ATOM"):
                continue
            chain = ln[21].strip()
            try:
                resnum = int(ln[22:26].strip())
            except:
                continue
            icode = ln[26].strip()
            key=(chain,resnum,icode)
            if key in seen:
                continue
            seen.add(key)
            b = float(ln[60:66].strip())
            res_plddt[(chain,resnum)] = b
    return res_plddt

# Find AF output files for a candidate
def find_outputs(cid, outdir):
    scores = sorted(glob.glob(os.path.join(outdir, f"{cid}*scores*.json")))
    pdbs   = sorted(glob.glob(os.path.join(outdir, f"{cid}*.pdb")))
    return (scores[0] if scores else None), (pdbs[0] if pdbs else None)

# Compute median PAE (H3 -> T) assuming chain order H:L:T
def median_pae_H3_to_T(pae, lenH, lenL, lenT, h3_resnums):
    offH = 0
    offL = lenH
    offT = lenH + lenL

    # Convert H3 resnums from "RFdiffusion numbering" to AF index
    h_idx = [r-1 for r in h3_resnums if 1 <= r <= lenH]
    if not h_idx:
        return None

    T_idx = np.arange(lenT) + offT

    vals=[]
    for i in (np.array(h_idx) + offH):
        vals.extend(pae[i, T_idx].tolist())

    return float(np.median(vals)) if vals else None

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--antigen_len", required=True, type=int)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    for col in ["candidate_id","VH_seq","VL_seq"]:
        if col not in df.columns:
            raise ValueError(f"Input CSV must contain {col}")

    results=[]

    for _, row in df.iterrows():
        cid = str(row["candidate_id"])
        scores_path, af_pdb_path = find_outputs(cid, args.outdir)

        rec = {"candidate_id": cid, "af_status": "missing"}

        if scores_path and af_pdb_path:

            # Load AF score JSON
            with open(scores_path,"r") as f:
                scores = json.load(f)

            pae = np.array(scores["pae"], dtype=float) if "pae" in scores else None

            lenH = len(str(row["VH_seq"]))
            lenL = len(str(row["VL_seq"]))
            lenT = args.antigen_len

            # Find H3 residue numbers from original complex PDB (if available)
            h3_resnums=[]
            if "pdb_path" in row and isinstance(row["pdb_path"],str):
                h3_resnums = parse_h3_resnums_from_remarks(row["pdb_path"])

            # Parse pLDDT values from AF predicted complex
            plddt_map = res_plddt_from_af_pdb(af_pdb_path)

            # Try chain letters used by AF (H first, then A)
            mean_plddt_H3=None
            for chain_id in ["H","A"]:
                vals=[plddt_map.get((chain_id, r), None) for r in h3_resnums]
                vals=[v for v in vals if v is not None]
                if vals:
                    mean_plddt_H3 = float(np.mean(vals))
                    break

            # Compute median PAE (H3 -> T)
            median_pae = None
            if pae is not None and h3_resnums:
                median_pae = median_pae_H3_to_T(pae, lenH, lenL, lenT, h3_resnums)

            rec.update({
                "af_status": "ok",
                "af_scores_json": scores_path,
                "af_pdb_path": af_pdb_path,
                "mean_plddt_H3": mean_plddt_H3,
                "median_pae_H3_to_T": median_pae,
            })

        results.append(rec)

    out = df.merge(pd.DataFrame(results), on="candidate_id", how="left")
    out.to_csv(args.out_csv, index=False)
    print("Wrote:", args.out_csv)
    print(out[["candidate_id","af_status","mean_plddt_H3","median_pae_H3_to_T"]].head())

if __name__ == "__main__":
    main()
