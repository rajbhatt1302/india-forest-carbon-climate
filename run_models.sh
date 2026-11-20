#!/usr/bin/env bash
# SAFE-CP gap filling model, its variants and the validation experiments (several hours on a CPU).
# Outputs go to data/ml_*; run_all.sh then computes the statistics and tables from them.
set -euo pipefail
cd "$(dirname "$0")"
G=data/gee_outputs
W=$G/grid05_window_results.csv
P=$G/grid05_predictors.csv
X=$G/grid05_predictors_extra.csv
TY=$G/grid05_predictors_type.csv
run() { python -W ignore analysis/safecp.py "$@"; }

# final model: committee + TabPFN, forest type, phenology and soil predictors, nested kriging decision
SAFECP_KRIGE=nested SAFECP_CHAIN=0 SAFECP_TYPE=$TY SAFECP_TAG=finalnest run $W $P $X data/ml_finalnest

# variants (Extended Data Fig. 8)
SAFECP_NO_TABPFN=1 SAFECP_KRIGE=always SAFECP_CHAIN=0 SAFECP_TAG=ablA run $W $P $X data/ml_ablation
SAFECP_NO_TABPFN=1 SAFECP_KRIGE=nested SAFECP_CHAIN=0 SAFECP_TAG=ablB run $W $P $X data/ml_ablation
SAFECP_NO_TABPFN=1 SAFECP_KRIGE=nested SAFECP_CHAIN=1 SAFECP_TARGETS=eN_LSTd,eP_LSTd,eN_LSTn SAFECP_TAG=ablC \
  run $W $P $X data/ml_ablation
SAFECP_NO_TABPFN=1 SAFECP_KRIGE=always SAFECP_CHAIN=0 SAFECP_TYPE=$TY SAFECP_TAG=typeC run $W $P $X data/ml_type
SAFECP_NO_TABPFN=1 SAFECP_KRIGE=always SAFECP_CHAIN=0 SAFECP_TAG=g025base \
  run $G/grid025_window_results.csv $G/grid025_predictors.csv data/ml_g025

# block size of the distance-matched cross-validation and the validation, data and algorithm experiment
NNDM_POT=data/ml_finalnest/safecp_potential_finalnest.csv python analysis/nndm_blocks.py $W $P eN_LSTd
python analysis/data_vs_algorithm.py data/stats/data_vs_algorithm.csv
