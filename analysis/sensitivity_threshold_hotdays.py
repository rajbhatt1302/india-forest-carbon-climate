"""Sensitivity analyses: natural forest threshold 0.52, and hot pre-monsoon days.
Usage: python sensitivity_threshold_hotdays.py <gee_dir> <out_dir>"""
import sys
import os
import numpy as np
import results as R
import pandas as pd

gd, od = sys.argv[1:3]
main = R.load(f'{gd}/grid05_window_results.csv')
geo = main[['gid', 'lon', 'lat', 'fIndia', 'nNat', 'nPlant', 'nCrop']]
NAME = {'N': 'Natural forest vs cropland', 'P': 'Plantation vs cropland', 'D': 'Natural forest replaced by plantation'}


def summ(df, col):
    return R.weighted_summary(df[col], R.weights(df))


# threshold 0.52
f52 = f'{gd}/grid05_window_results_nf052.csv'
import os
if os.path.exists(f52):
    t = pd.read_csv(f52).drop(columns=['system:index', '.geo'], errors='ignore')
    t = t.merge(main[['gid', 'fIndia']], on='gid', how='left')
    for y in R.Y_SEASON:
        for p in ('bN_', 'bP_', 'sN_', 'sP_', 'sD_'):
            t[p + y] = -9999
    t['n_r2'] = 0
    tmp = f'{od}/_nf052_tmp.csv'
    out = []
    R.load(f'{gd}/grid05_window_results.csv')
    s_main = {(y, k): R.weighted_summary(main[f'e{k}_{y}'], R.weights(main)) for y in ['LSTd', 'LSTn', 'ALB', 'ET', 'AGB'] for k in 'NPD'}
    t.to_csv(tmp, index=False); d52 = R.load(tmp); os.remove(tmp)
    for (y, k), s0 in s_main.items():
        s1 = R.weighted_summary(d52[f'e{k}_{y}'], R.weights(d52))
        out.append(dict(var=y, transition=NAME[k], main_mean=s0['mean'], main_lo=s0['lo'], main_hi=s0['hi'], main_n=s0['n'],
                        nf052_mean=s1['mean'], nf052_lo=s1['lo'], nf052_hi=s1['hi'], nf052_n=s1['n']))
    # net climate effect, 100 yr
    p = pd.read_csv(f'{gd}/grid05_predictors.csv').drop(columns=['system:index', '.geo'], errors='ignore')
    p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
    sc = pd.read_csv(f'{gd}/grid05_soc_nf052.csv').replace(-9999, np.nan)
    d = d52.merge(p[['gid', 'SWdown', 'AI']], on='gid', how='left').merge(sc[['gid', 'n_r4', 'bN_SOC', 'bP_SOC']], on='gid', how='left')
    d.index = d52.index
    okn = d['n_r4'] >= R.MIN_N
    d['eN_SOC'] = d['bN_SOC'].where(okn & (d['nNat'] >= R.MIN_CLASS) & (d['nCrop'] >= R.MIN_CLASS))
    d['eP_SOC'] = d['bP_SOC'].where(okn & (d['nPlant'] >= R.MIN_CLASS) & (d['nCrop'] >= R.MIN_CLASS))
    d['eD_SOC'] = (d['bP_SOC'] - d['bN_SOC']).where(okn & (d['nNat'] >= R.MIN_CLASS) & (d['nPlant'] >= R.MIN_CLASS))
    RS = np.where(d['AI'] > 1.0, 0.207, np.where(d['AI'] >= 0.65, 0.323, 0.440))
    K = 100.0 / (5.1e14 * 8.95e-14) * 1e4 / 1e3
    f9 = pd.read_csv(f'{od}/net_effect_national.csv')
    F9 = {'N': 'Natural forest', 'P': 'Plantation', 'D': 'Replacement'}
    for k in 'NPD':
        alb = d['SWdown'] * d[f'e{k}_ALB'] * 0.854 * K
        carb = d[f'e{k}_AGB'] * (1 + RS) * 0.47 * 44 / 12 + d[f'e{k}_SOC'] * 44 / 12
        net = (alb + carb).where(alb.notna() & carb.notna())
        s1 = R.weighted_summary(net, R.weights(d))
        s0 = f9[(f9['transition'] == F9[k]) & (f9['component'] == 'Net')].iloc[0]
        out.append(dict(var='Net CO2-eq 100 yr', transition=NAME[k], main_mean=s0['mean'], main_lo=s0['lo'], main_hi=s0['hi'],
                        main_n=s0['n'], nf052_mean=s1['mean'], nf052_lo=s1['lo'], nf052_hi=s1['hi'], nf052_n=s1['n']))
    out = pd.DataFrame(out); out.to_csv(f'{od}/sensitivity_nf052.csv', index=False)
    # same cases in the layout of stats_sensitivity.csv
    ext = out[out['var'].isin(['LSTd', 'ALB', 'ET', 'AGB'])]
    code = {v: k for k, v in NAME.items()}
    pd.DataFrame({'case': 'Natural forest threshold 0.52', 'group': ext['transition'].map(code), 'var': ext['var'],
                  'n': ext['nf052_n'], 'mean': ext['nf052_mean'], 'lo': ext['nf052_lo'], 'hi': ext['nf052_hi']}
                 ).to_csv(f'{od}/stats_sensitivity_extra.csv', index=False)
    print(out.round(3).to_string(index=False))

# hot days
fh = f'{gd}/grid05_hotdays.csv'
if os.path.exists(fh):
    R.load(f'{gd}/grid05_window_results.csv')
    h = pd.read_csv(fh).drop(columns=['system:index', '.geo', 'lon', 'lat'], errors='ignore').replace(-9999, np.nan)
    d = main[['gid', 'lon', 'lat', 'fIndia', 'nNat', 'nPlant', 'nCrop']].merge(h, on='gid', how='left')
    d.index = main.index
    res = []
    for reg, ys in (('h1', ['LSTd_hot90', 'LSTd_allMAM']), ('h2', ['LSTd_hot95'])):
        okn = d['n_' + reg] >= R.MIN_N
        for y in ys:
            d['eN_' + y] = d['bN_' + y].where(okn & (d['nNat'] >= R.MIN_CLASS) & (d['nCrop'] >= R.MIN_CLASS))
            d['eP_' + y] = d['bP_' + y].where(okn & (d['nPlant'] >= R.MIN_CLASS) & (d['nCrop'] >= R.MIN_CLASS))
            d['eD_' + y] = (d['bP_' + y] - d['bN_' + y]).where(okn & (d['nNat'] >= R.MIN_CLASS) & (d['nPlant'] >= R.MIN_CLASS))
    w = R.weights(d)
    for k in 'NPD':
        for y in ['LSTd_allMAM', 'LSTd_hot90', 'LSTd_hot95']:
            res.append(dict(transition=NAME[k], quantity=y, **R.weighted_summary(d[f'e{k}_{y}'], w)))
        for y in ['LSTd_hot90', 'LSTd_hot95']:
            # paired difference, same windows
            diff = d[f'e{k}_{y}'] - d[f'e{k}_LSTd_allMAM']
            s = R.weighted_summary(diff, w)
            v = diff.dropna()
            res.append(dict(transition=NAME[k], quantity=f'{y} minus allMAM', **s,
                            share_windows_more_negative=float((v < 0).mean()) if len(v) else np.nan))
    res = pd.DataFrame(res); res.to_csv(f'{od}/hotdays.csv', index=False)
    d.to_csv(f'{od}/hotdays_windows.csv', index=False)
    print(res.round(3).to_string(index=False))
