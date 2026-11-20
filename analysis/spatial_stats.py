"""LISA clusters and Moran's I correlogram of the window effects.
Usage: python spatial_stats.py <window_csv> <out_dir>"""
import sys
import os
import numpy as np
import pandas as pd
from libpysal.weights import KNN, DistanceBand
from esda.moran import Moran, Moran_Local
from statsmodels.stats.multitest import multipletests
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

win, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
df = R.load(win)
KM = 111.32
df['x_km'] = df['lon'] * KM * np.cos(np.deg2rad(22.0))
df['y_km'] = df['lat'] * KM

TARGETS = ['eN_LSTd', 'eP_LSTd', 'eD_LSTd', 'eN_LSTn', 'eN_ET', 'eN_AGB', 'eN_ALB']
lisa_rows, corr_rows = [], []
for t in TARGETS:
    d = df[df[t].notna()]
    xy = d[['x_km', 'y_km']].to_numpy()
    y = d[t].to_numpy()
    w = KNN.from_array(xy, k=8); w.transform = 'r'
    ml = Moran_Local(y, w, permutations=9999, seed=42)
    rej = multipletests(ml.p_sim, alpha=0.05, method='fdr_bh')[0]
    q = np.where(rej, ml.q, 0)                       # 1 HH, 2 LH, 3 LL, 4 HL, 0 not significant
    for gid, qq, ii, pp in zip(d['gid'], q, ml.Is, ml.p_sim):
        lisa_rows.append(dict(target=t, gid=gid, quad=int(qq), local_I=ii, p_sim=pp))
    gm = Moran(y, w, permutations=999)
    print(t, 'global I', round(gm.I, 3), 'p', gm.p_sim, '| significant clusters HH/LL/LH/HL:',
          [(q == k).sum() for k in (1, 3, 2, 4)], 'of', len(y))
    for d0 in range(0, 1500, 100):
        dist = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))
        nb = {i: [int(j) for j in np.where((dist[i] > d0) & (dist[i] <= d0 + 100))[0]] for i in range(len(y))}
        if sum(len(v) for v in nb.values()) < 50:
            continue
        from libpysal.weights import W
        keep = {i: v for i, v in nb.items()}
        ww = W(keep, silence_warnings=True); ww.transform = 'b'
        mi = Moran(y, ww, permutations=999)
        corr_rows.append(dict(target=t, d_lo=d0, d_mid=d0 + 50, I=mi.I, p=mi.p_sim, EI=mi.EI))
pd.DataFrame(lisa_rows).to_csv(f'{OUT}/lisa.csv', index=False)
corr = pd.DataFrame(corr_rows)
corr.to_csv(f'{OUT}/correlogram.csv', index=False)
for t in TARGETS:
    c = corr[corr['target'] == t]
    first_ns = c[(c['p'] > 0.05) | (c['I'] <= 0)]['d_mid']
    print(t, 'correlation range (first band with I not significantly > 0):', first_ns.min() if len(first_ns) else '>1500', 'km')
