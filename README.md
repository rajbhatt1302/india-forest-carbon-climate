# Natural forests, tree plantations and the climate of India

Code for the article *Can tree plantations match the cooling and carbon benefits of India's natural forests?*
(Rajarshi Bhattacharjee et al., submitted).

The study compares natural forest, plantation and other tree cover, and cropland in 1,357 windows of 0.5 degree
across India. Within each window, land surface temperature, albedo, evapotranspiration, aboveground biomass and soil
organic carbon of 1 km pixels are regressed on the fractions of natural forest and plantation, with cropland as the
reference and elevation as a covariate. The albedo change is converted to a CO2 equivalent and set against the
carbon stored. Windows without an observed estimate are filled with SAFE-CP, a spatially validated ensemble of tree
models and TabPFN with conformal prediction intervals.

The data set produced by the study (`research_output/`) and the layers read by this code (`input_layers/`) are
archived on Zenodo.

## Repository layout

| Folder | Content |
|---|---|
| `gee/` | Google Earth Engine scripts (JavaScript, Code Editor), numbered in run order |
| `analysis/` | Statistics, net climate effect, gap filling (SAFE-CP), validation and sensitivity analyses (Python) |
| `run_all.sh` | Statistics and tables from the archived input layers |
| `run_models.sh` | SAFE-CP model runs and validation experiments |

## How to reproduce

1. Download the Zenodo archive and copy the content of its `input_layers/` folder into a folder named `data/` at the
   root of this repository.
2. Install Python 3.13 and the packages in `requirements.txt` (`pip install -r requirements.txt`).
3. Run `bash run_all.sh` from the repository root. Results are written to `data/stats/`.
4. Optional: `bash run_models.sh` refits the gap filling model and its variants (several hours on a CPU).
   TabPFN downloads its model weights on first use.

To rebuild the inputs from the satellite data, run the scripts in `gee/` in the Earth Engine Code Editor in numerical
order. Replace `YOUR_PROJECT` with your Earth Engine project. Tables are exported to the Google Drive folder
`india_forest_climate`; place them in `data/gee_outputs/` (CMIP6 window tables in `data/cmip6_win/`).

## Earth Engine scripts

| Script | Output |
|---|---|
| `01_boundary.js` | Official boundary of India (Natural Earth, India point of view) |
| `02_landcover_154m.js`, `03_forest_change_154m.js` | Land cover and forest change fractions on a 1/720 degree grid |
| `04_surface_variables_1km.js` | LST, albedo, ET and biomass on the MODIS 1 km grid, 2019 to 2023 |
| `05_cmip6_annual.js`, `05b_cmip6_annual_by_scenario.js` | Annual temperature and aridity, NEX-GDDP-CMIP6 |
| `06_window_regression.js`, `07_window_regression_soc.js` | Window regressions (0.5 degree) |
| `08_window_regression_025deg.js` | Window regressions on the 0.25 degree grid |
| `09` to `12` | Map layers at 0.05 degree |
| `13_cmip6_window_tables.js`, `14_window_radiation.js` | CMIP6 and radiation per window |
| `15` to `19` | Predictors of the gap filling model |
| `20` to `22` | Sensitivity analyses (natural forest threshold 0.52, hottest days) |

## Analysis scripts

| Script | Output (in `data/stats/` unless stated) |
|---|---|
| `results.py` | Shared module: loads the window table, applies the support rules, bootstrap summaries |
| `stats_analysis.py` | National means, spatial autocorrelation, FDR significance, sensitivity of national means, drivers, attribution |
| `spatial_stats.py` | Local Moran's I (LISA) clusters and Moran's I correlogram |
| `gwr_analysis.py` | Multiscale geographically weighted regression (`data/ml/`) |
| `net_effect.py` | Net climate effect at 100 years and the albedo offset of the carbon benefit by aridity class |
| `time_horizon.py` | Net climate effect at 20, 100 and 500 years |
| `aridity_classes.py` | Effects by aridity class |
| `projected_warming.py` | Projected warming and the year it first exceeds the forest cooling, nationally and per window |
| `synthesis.py` | Classes of local cooling and net climate benefit, population per class, climate fingerprint |
| `sensitivity_threshold_hotdays.py` | Natural forest threshold 0.52 and hottest-day analyses |
| `national_summary_table.py` | National summary table of the article |
| `safecp.py` | SAFE-CP gap filling model (`data/ml_*/`) |
| `nndm_blocks.py` | Block size of the distance-matched cross-validation |
| `data_vs_algorithm.py`, `data_vs_algorithm_summary.py` | Effect of validation design, data and algorithm on model skill |
| `model_skill_table.py` | Skill of SAFE-CP and its variants |



Please cite the Github repository and the Zenodo archive (DOI: 10.5281/zenodo.24395478) when you use this code.
