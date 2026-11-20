"""Effect of validation design, predictors and learner on SAFE-CP skill.
Usage: python data_vs_algorithm.py <out_csv>"""
import sys
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold, KFold
import lightgbm as lgb
import xgboost as xgb
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

SEED = 7
G = os.path.join(ROOT, 'data/gee_outputs')
out = sys.argv[1]
df = R.load(f'{G}/grid05_window_results.csv')
p = pd.read_csv(f'{G}/grid05_predictors.csv').drop(columns=['system:index', '.geo'], errors='ignore')
p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
p['logAI'] = np.log10(p['AI'].clip(0.01, 20))
df = df.merge(p, on='gid', how='left')
e = pd.read_csv(f'{G}/grid05_predictors_extra.csv').drop(columns=['system:index', '.geo'], errors='ignore')
df = df.merge(e, on='gid', how='left')
BASE = ['logAI', 'P', 'PET', 'Tmax', 'Tmin', 'VPD', 'Pjjas', 'SWdown', 'elev', 'elevSD', 'slope',
        'LSTd_bg', 'LSTn_bg', 'ALB_bg', 'ET_bg', 'NDVI_rabi', 'NDVI_dry', 'SM', 'WIND']
ty = pd.read_csv(f'{G}/grid05_predictors_type.csv').drop(columns=['system:index', '.geo'], errors='ignore')
for v in ('EVG', 'FSEAS', 'CNC', 'CGSL', 'CAMP', 'NDVI_kharif'):
    ty[v + '_reg'] = ty[v + '_R100'].fillna(ty[v + '_R250'])
ty['EVG_reg'] = ty['EVG_reg'].clip(0, 1)
TPRED = ['EVG_reg', 'FSEAS_reg', 'CNC_reg', 'CGSL_reg', 'CAMP_reg', 'NDVI_kharif_reg', 'CLAY', 'SAND']
df = df.merge(ty[['gid'] + TPRED], on='gid', how='left')
sc = pd.read_csv(f'{G}/grid05_soc.csv').replace(-9999, np.nan)
df = df.merge(sc[['gid', 'n_r4', 'bN_SOC', 'bP_SOC']], on='gid', how='left')
RS = np.where(df['AI'] > 1.0, 0.207, np.where(df['AI'] >= 0.65, 0.323, 0.440))
K = 100.0 / (5.1e14 * 8.95e-14) * 1e4 / 1e3
okn = df['n_r4'] >= R.MIN_N
for k, cls in (('N', 'nNat'), ('P', 'nPlant')):
    soc = df[f'b{k}_SOC'].where(okn & (df[cls] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS))
    df[f'net{k}'] = df['SWdown'] * df[f'e{k}_ALB'] * 0.854 * K + df[f'e{k}_AGB'] * (1 + RS) * 0.47 * 44 / 12 + soc * 44 / 12


def learners():
    return {'Ridge': make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 20))),
            'RF': RandomForestRegressor(n_estimators=300, min_samples_leaf=5, max_features=0.5, n_jobs=2, random_state=SEED),
            'ET': ExtraTreesRegressor(n_estimators=300, min_samples_leaf=5, max_features=0.6, n_jobs=2, random_state=SEED),
            'LGBM': lgb.LGBMRegressor(n_estimators=400, learning_rate=0.03, num_leaves=8, min_child_samples=15, subsample=0.8,
                                      subsample_freq=1, colsample_bytree=0.7, reg_lambda=1.0, random_state=SEED, verbose=-1),
            'XGB': xgb.XGBRegressor(n_estimators=400, learning_rate=0.03, max_depth=3, min_child_weight=5, subsample=0.8,
                                    colsample_bytree=0.7, reg_lambda=1.0, random_state=SEED, n_jobs=2)}


def r2(y, p):
    return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def blocks(lon, lat, deg):
    return np.array([f'{a}_{b}' for a, b in zip(np.floor(lon / deg).astype(int), np.floor(lat / deg).astype(int))])


TARGETS = ['eN_LSTd', 'eP_LSTd', 'eN_LSTn', 'eN_ALB', 'eN_ET', 'eP_ET', 'eN_AGB', 'netN', 'netP']
RUNS = [('base', 'random'), ('base', 'block3'), ('base', 'block6'), ('extra', 'block3')]
rows = []
for t in TARGETS:
    for pset, design in RUNS:
        PRED = BASE if pset == 'base' else BASE + TPRED
        # same windows for every run
        ok = (df[t].notna() & df[BASE + TPRED].notna().all(1)).to_numpy()
        X, y = df.loc[ok, PRED].to_numpy(float), df.loc[ok, t].to_numpy(float)
        lon, lat = df.loc[ok, 'lon'].to_numpy(), df.loc[ok, 'lat'].to_numpy()
        if design == 'random':
            folds = list(KFold(n_splits=10, shuffle=True, random_state=SEED).split(X))
        else:
            folds = list(GroupKFold(n_splits=10).split(X, y, blocks(lon, lat, 3.0 if design == 'block3' else 6.0)))
        names = list(learners())
        pred = {m: np.zeros(len(y)) for m in names}
        for tr, te in folds:
            L = learners()
            for m in names:
                pred[m][te] = L[m].fit(X[tr], y[tr]).predict(X[te])
        pred['Committee'] = np.mean([pred[m] for m in ('RF', 'ET', 'LGBM', 'XGB')], axis=0)
        for m, pv in pred.items():
            rows.append(dict(target=t, predictors=pset, design=design, model=m, n=len(y), r2=r2(y, pv)))
        print(t, pset, design, {m: round(r2(y, pv), 3) for m, pv in pred.items()}, flush=True)
        pd.DataFrame(rows).to_csv(out, index=False)
print('done')
