"""National means, Moran's I, FDR significance, sensitivity, drivers and attribution.
Usage: python stats_analysis.py <window_csv> <predictor_csv> <out_dir>"""
import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

win, prd, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(11)
NBOOT = 2000
BLOCK = 6.0

df = R.load(win)
p = pd.read_csv(prd).drop(columns=['system:index', '.geo'], errors='ignore')
p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
df = df.merge(p, on='gid', how='left')
df['w'] = df['fIndia'] * np.cos(np.deg2rad(df['lat']))
df['block'] = (np.floor(df['lon'] / BLOCK).astype(int).astype(str) + '_' + np.floor(df['lat'] / BLOCK).astype(int).astype(str))
df['logAI'] = np.log10(df['AI'].clip(0.01, 20))
VARS = ['LSTd', 'LSTn', 'ALB', 'ET', 'AGB', 'LSTd_MAM', 'LSTn_MAM', 'LSTd_JJAS', 'LSTd_ONDJF']
GR = [('N', 'Natural forest vs cropland'), ('P', 'Plantation vs cropland'), ('D', 'Natural forest replaced by plantation')]


def wmean(v, w):
    return float((v * w).sum() / w.sum())


def boot_iid(v, w):
    idx = rng.integers(0, len(v), (NBOOT, len(v)))
    b = (v[idx] * w[idx]).sum(1) / w[idx].sum(1)
    return np.percentile(b, [2.5, 97.5])


def boot_block(v, w, blocks):
    ub = np.unique(blocks)
    pos = {b: np.where(blocks == b)[0] for b in ub}
    out = np.empty(NBOOT)
    for i in range(NBOOT):
        pick = rng.choice(ub, len(ub), replace=True)
        ii = np.concatenate([pos[b] for b in pick])
        out[i] = (v[ii] * w[ii]).sum() / w[ii].sum()
    return np.percentile(out, [2.5, 97.5])


def morans_i(v, lon, lat, k=8, nperm=999):
    x = np.c_[lon * np.cos(np.deg2rad(lat)), lat]
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    nb = np.argsort(d, 1)[:, :k]
    z = v - v.mean()
    def I(zz):
        return (len(zz) / (zz ** 2).sum()) * (zz[:, None] * zz[nb]).mean(1).sum() / len(zz)
    i0 = I(z)
    perm = np.array([I(rng.permutation(z)) for _ in range(nperm)])
    pval = (np.sum(perm >= i0) + 1) / (nperm + 1)
    return float(i0), float(pval)


# national statistics
rows = []
for y in VARS:
    for g, gl in GR:
        e, z = df[f'e{g}_{y}'], df[f'z{g}_{y}']
        m = e.notna() & df['w'].notna()
        v, w, bl = e[m].to_numpy(), df.loc[m, 'w'].to_numpy(), df.loc[m, 'block'].to_numpy()
        if len(v) < 20:
            continue
        mi, mp = morans_i(v, df.loc[m, 'lon'].to_numpy(), df.loc[m, 'lat'].to_numpy())
        lo_i, hi_i = boot_iid(v, w)
        lo_b, hi_b = boot_block(v, w, bl)
        pv = 2 * stats.norm.sf(np.abs(z[m].to_numpy()))
        rej = multipletests(pv, alpha=0.05, method='fdr_bh')[0]
        rows.append(dict(var=y, group=g, label=gl, n=len(v), n_blocks=len(np.unique(bl)), mean=wmean(v, w),
                         ci_iid_lo=lo_i, ci_iid_hi=hi_i, ci_block_lo=lo_b, ci_block_hi=hi_b,
                         morans_I=mi, morans_p=mp, n_eff=len(v) * (1 - mi) / (1 + mi),
                         share_sig_neg=float(np.mean(rej & (v < 0))), share_sig_pos=float(np.mean(rej & (v > 0))),
                         share_sig_raw=float(np.mean(pv < 0.05))))
nat = pd.DataFrame(rows)
nat.to_csv(f'{OUT}/stats_national.csv', index=False)
print(nat[['var', 'group', 'n', 'mean', 'ci_iid_lo', 'ci_iid_hi', 'ci_block_lo', 'ci_block_hi', 'morans_I', 'n_eff',
           'share_sig_neg', 'share_sig_pos']].round(3).to_string())

# sensitivity of national means
SENS = []
base = R.load(win)


def recompute(min_n, min_class):
    d = base.copy()
    out = {}
    for y in ['LSTd', 'ET', 'ALB', 'AGB']:
        n = d['n_' + R.REG_OF[y]] >= min_n
        out[('N', y)] = d['bN_' + y].where(n & (d['nNat'] >= min_class) & (d['nCrop'] >= min_class))
        out[('P', y)] = d['bP_' + y].where(n & (d['nPlant'] >= min_class) & (d['nCrop'] >= min_class))
        out[('D', y)] = (d['bP_' + y] - d['bN_' + y]).where(n & (d['nNat'] >= min_class) & (d['nPlant'] >= min_class))
    return out


cases = [('Baseline (n >= 200, class >= 30)', 200, 30, 'area', None),
         ('Class support >= 10 pixels', 200, 10, 'area', None),
         ('Class support >= 60 pixels', 200, 60, 'area', None),
         ('Regression n >= 400', 400, 30, 'area', None),
         ('Equal window weights', 200, 30, 'equal', None),
         ('Without windows above 1,000 m', 200, 30, 'area', 'elev'),
         ('Without the Northeast', 200, 30, 'area', 'ne'),
         ('Only windows fully inside India', 200, 30, 'area', 'full')]
for name, mn, mc, wt, sub in cases:
    eff = recompute(mn, mc)
    keep = pd.Series(True, index=df.index)
    if sub == 'elev':
        keep = df['elev'] <= 1000
    elif sub == 'ne':
        keep = ~((df['lon'] >= 89.5) & (df['lat'] >= 21.5))
    elif sub == 'full':
        keep = df['fIndia'] >= 0.99
    ww = df['w'] if wt == 'area' else pd.Series(1.0, index=df.index)
    for (g, y), e in eff.items():
        m = e.notna() & keep & ww.notna()
        v, w, bl = e[m].to_numpy(), ww[m].to_numpy(), df.loc[m, 'block'].to_numpy()
        lo, hi = boot_block(v, w, bl)
        SENS.append(dict(case=name, group=g, var=y, n=len(v), mean=wmean(v, w), lo=lo, hi=hi))
sens = pd.DataFrame(SENS)
sens.to_csv(f'{OUT}/stats_sensitivity.csv', index=False)
print(sens[sens['var'] == 'LSTd'].round(3).to_string())

# drivers
DRV = ['logAI', 'SWdown', 'VPD', 'LSTd_bg', 'elev', 'Pjjas']
DLAB = {'logAI': 'Aridity index (log)', 'SWdown': 'Solar radiation', 'VPD': 'Vapour pressure deficit',
        'LSTd_bg': 'Background daytime LST', 'elev': 'Elevation', 'Pjjas': 'Monsoon share of rainfall'}
drows, vifs = [], []
for y in ['LSTd', 'LSTn', 'ET', 'AGB']:
    for g, gl in GR:
        col = f'e{g}_{y}'
        d = df[[col, 'block', 'w'] + DRV].dropna()
        X = (d[DRV] - d[DRV].mean()) / d[DRV].std()
        Y = (d[col] - d[col].mean()) / d[col].std()
        Xa = np.c_[np.ones(len(X)), X.to_numpy()]
        beta = np.linalg.lstsq(Xa, Y.to_numpy(), rcond=None)[0][1:]
        r2 = 1 - np.sum((Y.to_numpy() - Xa @ np.linalg.lstsq(Xa, Y.to_numpy(), rcond=None)[0]) ** 2) / np.sum((Y - Y.mean()) ** 2)
        ub = d['block'].unique(); pos = {b: np.where(d['block'].to_numpy() == b)[0] for b in ub}
        bs = []
        for _ in range(1000):
            ii = np.concatenate([pos[b] for b in rng.choice(ub, len(ub), replace=True)])
            bs.append(np.linalg.lstsq(Xa[ii], Y.to_numpy()[ii], rcond=None)[0][1:])
        bs = np.array(bs)
        for j, nm in enumerate(DRV):
            lo, hi = np.percentile(bs[:, j], [2.5, 97.5])
            drows.append(dict(var=y, group=g, predictor=nm, label=DLAB[nm], beta=beta[j], lo=lo, hi=hi, r2=r2, n=len(d),
                              robust=bool(lo > 0 or hi < 0)))
        if y == 'LSTd' and g == 'N':
            for j, nm in enumerate(DRV):
                vifs.append(dict(predictor=nm, vif=variance_inflation_factor(Xa, j + 1)))
drv = pd.DataFrame(drows)
drv.to_csv(f'{OUT}/stats_drivers.csv', index=False)
pd.DataFrame(vifs).to_csv(f'{OUT}/stats_drivers_vif.csv', index=False)
print(drv[drv['var'] == 'LSTd'].round(3).to_string())
print(pd.DataFrame(vifs).round(2))

# attribution of LSTd to ET and albedo
arows = []
for g, gl in GR:
    d = df[[f'e{g}_LSTd', f'e{g}_ET', f'e{g}_ALB', 'block']].dropna()
    y = d[f'e{g}_LSTd'].to_numpy()
    X = np.c_[np.ones(len(d)), d[f'e{g}_ET'].to_numpy() / 100.0, d[f'e{g}_ALB'].to_numpy() / 0.01]
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    ub = d['block'].unique(); pos = {k: np.where(d['block'].to_numpy() == k)[0] for k in ub}
    bs = np.array([np.linalg.lstsq(X[ii], y[ii], rcond=None)[0] for ii in
                   (np.concatenate([pos[k] for k in rng.choice(ub, len(ub), replace=True)]) for _ in range(1000))])
    # contribution = slope x mean change
    mET, mAL = X[:, 1].mean(), X[:, 2].mean()
    cET, cAL = bs[:, 1] * mET, bs[:, 2] * mAL
    arows.append(dict(group=g, label=gl, n=len(d), mean_LSTd=y.mean(),
                      slope_ET_per100mm=b[1], slope_ET_lo=np.percentile(bs[:, 1], 2.5), slope_ET_hi=np.percentile(bs[:, 1], 97.5),
                      slope_ALB_per001=b[2], slope_ALB_lo=np.percentile(bs[:, 2], 2.5), slope_ALB_hi=np.percentile(bs[:, 2], 97.5),
                      contrib_ET=b[1] * mET, contrib_ET_lo=np.percentile(cET, 2.5), contrib_ET_hi=np.percentile(cET, 97.5),
                      contrib_ALB=b[2] * mAL, contrib_ALB_lo=np.percentile(cAL, 2.5), contrib_ALB_hi=np.percentile(cAL, 97.5),
                      intercept=b[0], r=np.corrcoef(d[f'e{g}_ET'], d[f'e{g}_ALB'])[0, 1]))
att = pd.DataFrame(arows)
att.to_csv(f'{OUT}/stats_attribution.csv', index=False)
print(att.round(3).T.to_string())
