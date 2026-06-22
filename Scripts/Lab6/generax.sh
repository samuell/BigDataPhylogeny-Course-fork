#!/usr/bin/env bash
# DTL reconciliation with GeneRax
#
# Runs three separate GeneRax jobs against the same fixed species tree:
#   1. baseline_UndatedDTL    - the 5 real single-copy core genes, full DTL model
#   2. engineered_UndatedDTL  - the 2 engineered teaching families, full DTL model
#   3. engineered_UndatedDL   - the 2 engineered teaching families, D+L only (no transfers)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SPECIES_TREE="data/species_tree/species_tree_rn.nwk"

echo ">> [1/3] Baseline core gene set, UndatedDTL model"
generax \
    --families data/core_gene_set/families_baseline.txt \
    --species-tree "$SPECIES_TREE" \
    --rec-model UndatedDTL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/baseline_UndatedDTL

echo
echo ">> [2/3] Engineered families, UndatedDTL model (transfers allowed)"
generax \
    --families data/dtl_engineered/families_engineered.txt \
    --species-tree "$SPECIES_TREE" \
    --rec-model UndatedDTL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/engineered_UndatedDTL

echo
echo ">> [3/3] Engineered families, UndatedDL model (no transfers allowed)"
generax \
    --families data/dtl_engineered/families_engineered.txt \
    --species-tree "$SPECIES_TREE" \
    --rec-model UndatedDL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/engineered_UndatedDL

echo
echo "Done. Look inside:"
echo "  results/baseline_UndatedDTL/reconciliations/"
echo "  results/engineered_UndatedDTL/reconciliations/"
echo "  results/engineered_UndatedDL/reconciliations/"
