import os, sys
import pandas as pd

AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
}

def list_chains(pdb_path):
    s=set()
    with open(pdb_path,"r") as f:
        for ln in f:
            if ln.startswith("ATOM"):
                s.add(ln[21].strip())
    return sorted(s)

def pdb_chain_sequence(pdb_path, chain_id):
    seen=set()
    seq=[]
    with open(pdb_path,"r") as f:
        for ln in f:
            if not ln.startswith("ATOM"): 
                continue
            if ln[21].strip()!=chain_id:
                continue
            resname=ln[17:20].strip()
            resseq=ln[22:26].strip()
            icode=ln[26].strip()
            key=(resseq,icode)
            if key in seen:
                continue
            seen.add(key)
            seq.append(AA3_TO_1.get(resname,"X"))
    return "".join(seq)

def main():
    if len(sys.argv)<5:
        print("Usage: python make_fastas.py <csv> <antigen_pdb> <outdir> <chain_or_AUTO>")
        sys.exit(1)

    csv_path, antigen_pdb, outdir, chain = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)

    chains=list_chains(antigen_pdb)
    if chain=="AUTO":
        chain = chains[0] if len(chains)==1 else chains[-1]

    antigen_seq=pdb_chain_sequence(antigen_pdb, chain)
    print(f"Antigen chain={chain} len={len(antigen_seq)} X={antigen_seq.count('X')}")

    df=pd.read_csv(csv_path)
    for c in ["candidate_id","VH_seq","VL_seq"]:
        if c not in df.columns: 
            raise ValueError(f"Missing {c}")

    count=0
    for _,r in df.iterrows():
        cid=str(r["candidate_id"])
        vh=r["VH_seq"]
        vl=r["VL_seq"]
        fout=os.path.join(outdir,f"{cid}.fasta")
        with open(fout,"w") as f:
            f.write(f">{cid}\n{vh}:{vl}:{antigen_seq}\n")
        count+=1

    print("Wrote FASTAs:", count)

if __name__=="__main__":
    main()