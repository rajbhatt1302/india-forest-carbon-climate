"""MGWR of the daytime LST effects on standardized covariates.
Usage: python gwr_analysis.py <window_csv> <predictor_csv> [extra_csv] <out_dir>"""
import sys
import os
import numpy as np
import pandas as pd
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

args = sys.argv[1:]
win, prd = args[0], args[1]
extra = args[2] if len(args) == 4 else None
OUT = args[-1]
os.makedirs(OUT, exist_ok=True)
df = R.load(win)
p = pd.read_csv(prd).drop(columns=['system:index', '.geo'], errors='ignore')
p['logAI'] = np.log10(p['AI'].clip(0.01, 20))
df = df.merge(p, on='gid', how='left')
COV = ['LSTd_bg', 'VPD', 'logAI', 'elev']
if extra:
    e = pd.read_csv(extra)
    df = df.merge(e, on='gid', how='left')
    COV += ['NDVI_rabi']
KM = 111.32
rows, summ = [], []
for t in ['eN_LSTd', 'eP_LSTd', 'eD_LSTd']:
    d = df[df[t].notna() & df[COV].notna().all(1)].reset_index(drop=True)
    coords = np.c_[d['lon'] * KM * np.cos(np.deg2rad(22.0)), d['lat'] * KM]
    y = ((d[t] - d[t].mean()) / d[t].std()).to_numpy().reshape(-1, 1)
    X = ((d[COV] - d[COV].mean()) / d[COV].std()).to_numpy()
    sel = Sel_BW(coords, y, X, multi=True)
    bws = sel.search(multi_bw_min=[30])
    res = MGWR(coords, y, X, sel).fit()
    tcrit = res.critical_tval()
    filt = res.filter_tvals()
    for i in range(len(d)):
        r = dict(target=t, gid=d.loc[i, 'gid'], lon=d.loc[i, 'lon'], lat=d.loc[i, 'lat'], local_r2=np.nan)
        for j, c in enumerate(['intercept'] + COV):
            r['b_' + c] = res.params[i, j]
            r['sig_' + c] = bool(filt[i, j] != 0)
        rows.append(r)
    summ.append(dict(target=t, n=len(d), R2=res.R2, adjR2=res.adj_R2, AICc=res.aicc,
                     **{'bw_' + c: b for c, b in zip(['intercept'] + COV, bws)},
                     **{'tcrit_' + c: tc for c, tc in zip(['intercept'] + COV, tcrit)}))
    print(t, 'MGWR R2', round(res.R2, 3), 'adj', round(res.adj_R2, 3), 'bandwidths', dict(zip(['intercept'] + COV, bws)))
    # global OLS for comparison
    Xa = np.c_[np.ones(len(X)), X]
    b = np.linalg.lstsq(Xa, y, rcond=None)[0]
    print('   global OLS R2', round(1 - ((y - Xa @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum(), 3))
tag = 'extra' if extra else 'base'
pd.DataFrame(rows).to_csv(f'{OUT}/mgwr_local_{tag}.csv', index=False)
pd.DataFrame(summ).to_csv(f'{OUT}/mgwr_summary_{tag}.csv', index=False)
print('saved')
