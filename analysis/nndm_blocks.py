"""Block size for distance-matched cross-validation (NNDM, Mila et al. 2022).
Usage: python nndm_blocks.py <window_csv> <predictor_csv> <target> [block sizes]"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp
from sklearn.model_selection import GroupKFold
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

win, prd, t = sys.argv[1:4]
sizes = [float(s) for s in (sys.argv[4] if len(sys.argv) > 4 else '1,1.5,2,3,4,6').split(',')]
df = R.load(win)
p = pd.read_csv(prd).drop(columns=['system:index', '.geo'], errors='ignore')
df = df.merge(p, on='gid', how='left')
KM = 111.32
xy = np.column_stack([df['lon'] * KM * np.cos(np.deg2rad(df['lat'])), df['lat'] * KM])
obs = df[t].notna().to_numpy()
gap = (df[t].isna() & df['AI'].notna() & (df['n_r1'] >= 0)).to_numpy()
if os.environ.get('NNDM_POT'):
    pot = pd.read_csv(os.environ['NNDM_POT'])[['gid', 'm_' + t]]
    gap = df[['gid']].merge(pot, on='gid', how='left')['m_' + t].notna().to_numpy()
dpred = cKDTree(xy[obs]).query(xy[gap])[0]
print(f'{t}: observed {obs.sum()}, gap {gap.sum()}; prediction distance median {np.median(dpred):.0f} km, '
      f'p90 {np.percentile(dpred, 90):.0f} km')
xo = xy[obs]
for b in sizes:
    grp = (np.floor(df.loc[obs, 'lon'] / b).astype(int).astype(str) + '_' +
           np.floor(df.loc[obs, 'lat'] / b).astype(int).astype(str)).to_numpy()
    d = np.zeros(obs.sum())
    for tr, te in GroupKFold(n_splits=10).split(xo, groups=grp):
        d[te] = cKDTree(xo[tr]).query(xo[te])[0]
    ks = ks_2samp(d, dpred).statistic
    print(f'  block {b:>4} deg: CV distance median {np.median(d):.0f} km, p90 {np.percentile(d, 90):.0f} km, KS {ks:.3f}')
