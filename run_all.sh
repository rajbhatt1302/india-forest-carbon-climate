#!/usr/bin/env bash
# Statistics and tables of the article from the archived input layers.
# Run from the repository root after copying input_layers/ of the Zenodo archive to data/.
set -euo pipefail
cd "$(dirname "$0")"
G=data/gee_outputs
S=data/stats
W=$G/grid05_window_results.csv
P=$G/grid05_predictors.csv

python analysis/stats_analysis.py $W $P $S
python analysis/spatial_stats.py $W $S
python analysis/gwr_analysis.py $W $P data/ml
python analysis/net_effect.py $W $P $G/grid05_soc.csv $S/net_effect_national.csv
python analysis/time_horizon.py $W $P $G/grid05_soc.csv $S/time_horizon_sensitivity.csv
python analysis/aridity_classes.py $W $P $S/effects_by_aridity.csv
python analysis/projected_warming.py $W data/cmip6_win $S/stats_national.csv $S/projected_warming.csv \
  $S/projected_warming_windows.csv
python analysis/synthesis.py $W $P $G/grid05_soc.csv $G/grid05_predictors_extra.csv \
  data/ml_finalnest/safecp_potential_finalnest.csv $S/synthesis_classes.csv $S/climate_fingerprint.csv
python analysis/sensitivity_threshold_hotdays.py $G $S
python analysis/model_skill_table.py data $S/model_skill_table.csv
python analysis/data_vs_algorithm_summary.py
python analysis/national_summary_table.py
