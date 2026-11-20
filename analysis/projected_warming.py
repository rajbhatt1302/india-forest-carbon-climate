"""Projected warming (NEX-GDDP-CMIP6, 10 models) against the daytime cooling of natural forest.
Usage: python projected_warming.py <window_csv> <cmip6_dir> <stats_national_csv> <out_csv> <out_windows_csv>"""
import sys
import glob
import numpy as np
import pandas as pd
import results as R

win, cdir, natcsv, out, outw = sys.argv[1:6]
df = R.load(win).set_index('gid')
w = R.weights(df.reset_index()).set_axis(df.index)
MODELS = sorted(set(f.split('/')[-1].rsplit('_', 1)[0] for f in glob.glob(cdir + '/*_hist.csv')))
assert len(MODELS) == 10, MODELS
YH, YF = range(1995, 2015), range(2015, 2100)
SCEN = ['ssp245', 'ssp585']


def tab(m, s, years):
    d = pd.read_csv(f'{cdir}/{m}_{s}.csv').set_index('gid')
    return pd.DataFrame({y: d[f'T_{s}_{y}'] for y in years})


base = {m: tab(m, 'hist', YH) for m in MODELS}
Tcm = pd.concat([base[m].mean(1) for m in MODELS], axis=1).mean(1)
nat_series, local = {}, {}
for s in SCEN:
    A = [tab(m, s, YF).sub(base[m].mean(1), axis=0) for m in MODELS]
    loc = sum(A) / len(A)
    local[s] = loc.T.rolling(11, center=True, min_periods=6).mean().T
    ok = A[0].notna().all(1)
    ww = w.reindex(A[0].index).where(ok, 0).fillna(0)
    nat_series[s] = pd.DataFrame({m: (a.mul(ww, axis=0).sum() / ww.sum()) for m, a in zip(MODELS, A)})
nat = pd.read_csv(natcsv)
nat = nat[nat['var'] == 'LSTd'].drop_duplicates('group').set_index('group')

# space for time
d = pd.DataFrame({'cool': -df['eN_LSTd'], 'T': Tcm.reindex(df.index), 'lon': df['lon'], 'lat': df['lat']}).dropna()
X1 = np.column_stack([np.ones(len(d)), d['T']])
beta = np.linalg.lstsq(X1, d['cool'], rcond=None)[0]
blk = (np.floor(d['lon'] / 6).astype(int).astype(str) + '_' + np.floor(d['lat'] / 6).astype(int).astype(str)).to_numpy()
ub = np.unique(blk); rng = np.random.default_rng(1); bs = []
for _ in range(2000):
    pick = np.concatenate([np.where(blk == b)[0] for b in rng.choice(ub, len(ub))])
    bs.append(np.linalg.lstsq(X1[pick], d['cool'].to_numpy()[pick], rcond=None)[0][1])
slo, shi = np.percentile(bs, [2.5, 97.5])
dT = {s: nat_series[s].loc[2081:2099].mean().mean() for s in SCEN}

# local crossing years
years = pd.DataFrame(index=df.index)
for s in SCEN:
    loc = local[s].reindex(df.index)
    cool = -df['eN_LSTd']
    yr = pd.Series(np.nan, index=df.index)
    for g in df.index[cool.notna() & loc.notna().all(1)]:
        hit = loc.columns[(loc.loc[g] >= cool[g]).to_numpy()]
        yr[g] = hit.min() if len(hit) else 9999
    years[f'year_warming_exceeds_cooling_{s}'] = yr

rows = []
for s in SCEN:
    X = nat_series[s]
    for p0, p1 in ((2041, 2060), (2081, 2099)):
        v = X.loc[p0:p1].mean()
        rows.append(dict(item=f'warming_{s}_{p0}_{p1}', mean=v.mean(), lo=v.min(), hi=v.max()))
    rm = X.mean(1).rolling(11, center=True).mean()
    for g in 'NPD':
        thr = abs(nat.loc[g, 'mean']); yrs = rm[rm >= thr].index
        rows.append(dict(item=f'year_national_warming_exceeds_{g}_{s}', mean=yrs.min() if len(yrs) else np.nan))
for s in SCEN:
    yr = years[f'year_warming_exceeds_cooling_{s}']
    n = int(yr.notna().sum())
    rows.append(dict(item=f'windows_{s}', mean=n))
    rows.append(dict(item=f'share_windows_not_exceeded_{s}', mean=float((yr == 9999).sum() / n)))
    rows.append(dict(item=f'area_share_not_exceeded_{s}', mean=float(w[yr == 9999].sum() / w[yr.notna()].sum())))
    rows.append(dict(item=f'median_year_exceeded_{s}', mean=float(yr[(yr < 9999)].median())))
rows.append(dict(item='slope_cooling_per_degC', mean=beta[1], lo=slo, hi=shi))
for s in SCEN:
    rows.append(dict(item=f'implied_change_in_cooling_2081_2099_{s}', mean=beta[1] * dT[s], lo=slo * dT[s], hi=shi * dT[s]))
pd.DataFrame(rows).to_csv(out, index=False)
years.dropna(how='all').reset_index().to_csv(outw, index=False)
print(pd.DataFrame(rows).round(3).to_string(index=False))
