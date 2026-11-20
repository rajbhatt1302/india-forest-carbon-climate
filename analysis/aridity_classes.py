"""Effects by aridity class and their Spearman correlation with the aridity index.
Usage: python aridity_classes.py <window_csv> <predictor_csv> <out_csv>"""
import sys
import numpy as np
import pandas as pd
from scipy import stats
import results as R

csv, pcsv, out = sys.argv[1:4]
df = R.load(csv)
p = pd.read_csv(pcsv).drop(columns=['system:index', '.geo'], errors='ignore')
df = df.merge(p[['gid', 'AI']], on='gid', how='left')
w = R.weights(df)
# no arid window has both covers
CL = [(0.2, 0.5, 'Semiarid'), (0.5, 0.65, 'Dry subhumid'), (0.65, 1.0, 'Humid'), (1.0, 99, 'Very humid')]
df['acl'] = np.nan
for i, (lo, hi, _) in enumerate(CL):
    df.loc[(df['AI'] >= lo) & (df['AI'] < hi), 'acl'] = i
GROUPS = [('eN_', 'Natural forest'), ('eP_', 'Plantation'), ('eD_', 'Replacement')]
rows = []
for y in ['LSTd', 'LSTn', 'ET', 'AGB']:
    for g, short in GROUPS:
        for i in range(len(CL)):
            m = df['acl'] == i
            rows.append(dict(var=y, group=short, cls=CL[i][2], **R.weighted_summary(df.loc[m, g + y], w[m])))
        d = df[[g + y, 'AI']].dropna()
        rho, pv = stats.spearmanr(d['AI'], d[g + y])
        rows.append(dict(var=y, group=short, cls='spearman', mean=rho, lo=pv, hi=np.nan, n=len(d)))
pd.DataFrame(rows).to_csv(out, index=False)
