cat > scripts/run_task1.sh << 'BASH'
#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./scripts/run_task1.sh inputs/HER2_h3_extracted.csv inputs/antigen_crop_15A.pdb AUTO 75 outputs/HER2
CSV="${1:?CSV path required}"
ANTIGEN_PDB="${2:?antigen PDB path required}"
ANTIGEN_CHAIN="${3:-AUTO}"
ANTIGEN_LEN="${4:?antigen length required}"
OUT_SUBDIR="${5:-outputs/run}"

FASTA_DIR="fastas"
OUT_DIR="$OUT_SUBDIR"

mkdir -p "$FASTA_DIR" "$OUT_DIR"

echo "[1/3] Make multimer FASTAs -> ${FASTA_DIR}"
python3 scripts/make_fastas.py "$CSV" "$ANTIGEN_PDB" "$FASTA_DIR" "$ANTIGEN_CHAIN"

echo "[2/3] Run ColabFold AF2-multimer -> ${OUT_DIR}"
colabfold_batch \
  --model-type alphafold2_multimer_v3 \
  --num-recycle 3 \
  --num-models 5 \
  "$FASTA_DIR" "$OUT_DIR"

echo "[3/3] Parse pLDDT/PAE -> ${OUT_DIR}/task1_structural_scores.csv"
python3 scripts/parse_plddt_pae.py \
  --csv "$CSV" \
  --outdir "$OUT_DIR" \
  --out_csv "$OUT_DIR/task1_structural_scores.csv" \
  --antigen_len "$ANTIGEN_LEN"

echo "DONE: ${OUT_DIR}/task1_structural_scores.csv"
BASH

chmod +x scripts/run_task1.sh