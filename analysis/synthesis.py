"""Cooling and net benefit classes, population per class and climate fingerprint.
Usage: python synthesis.py <window_csv> <predictor_csv> <soc_csv> <extra_csv> <potential_csv> <out_classes_csv> <out_fingerprint_csv>"""
import sys
import numpy as np
import pandas as pd
import results as R

win, prd, socf, extf, potf, outc, outf = sys.argv[1:8]
TA, K = 0.854, 100.0 / (5.1e14 * 8.95e-14) * 1e4 / 1e3
CF, C2CO2 = 0.47, 44.0 / 12.0

df = R.load(win)
p = pd.read_csv(prd)[['gid', 'AI', 'SWdown']]
p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
sc = pd.read_csv(socf).replace(-9999, np.nan)
ex = pd.read_csv(extf)
pot = pd.read_csv(potf)
df = df.merge(p, on='gid', how='left').merge(sc[['gid', 'n_r4', 'bN_SOC', 'bP_SOC']], on='gid', how='left') \
       .merge(ex[['gid', 'POP']], on='gid', how='left').merge(pot.drop(columns=['lon', 'lat']), on='gid', how='left')
RS = np.where(df['AI'] > 1.0, 0.207, np.where(df['AI'] >= 0.65, 0.323, 0.440))
okn = df['n_r4'] >= R.MIN_N
socN = df['bN_SOC'].where(okn & (df['nNat'] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS))
df['netN_obs'] = df['SWdown'] * df['eN_ALB'] * TA * K + df['eN_AGB'] * (1 + RS) * CF * C2CO2 + socN * C2CO2
# observed, else SAFE-CP
df['cool'] = (-df['eN_LSTd']).fillna(-df['m_eN_LSTd'])
df['net'] = df['netN_obs'].fillna(df['m_netN'])
df['modelled'] = df['eN_LSTd'].isna() | df['netN_obs'].isna()
unc_cool = df['eN_LSTd'].isna() & (df['lo_eN_LSTd'] < 0) & (df['hi_eN_LSTd'] > 0)
unc_net = df['netN_obs'].isna() & (df['lo_netN'] < 0) & (df['hi_netN'] > 0)
df['uncertain'] = unc_cool | unc_net
df['cx'] = np.digitize(df['cool'], [2.0, 4.0])            # 0 low, 1 mid, 2 high cooling
df['ny'] = np.digitize(df['net'], [-25.0, 25.0])          # 0 net warming, 1 neutral, 2 net cooling
ok = df['cool'].notna() & df['net'].notna()
conf = ok & ~df['uncertain']
tab = df[conf].groupby(['ny', 'cx']).agg(windows=('gid', 'size'), pop_millions=('POP', lambda s: s.sum() / 1e6),
                                       modelled=('modelled', 'mean'))
tab.to_csv(outc)
print(tab)
print('excluded (uncertain sign):', int(df.loc[ok, 'uncertain'].sum()), 'windows,',
      round(df.loc[ok & df['uncertain'], 'POP'].sum() / 1e6), 'million people')

# climate fingerprint (observed windows only)
w = R.weights(df)
soc = {'N': socN,
       'P': df['bP_SOC'].where(okn & (df['nPlant'] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS)),
       'D': (df['bP_SOC'] - df['bN_SOC']).where(okn & (df['nNat'] >= R.MIN_CLASS) & (df['nPlant'] >= R.MIN_CLASS))}
AXES = ['Daytime cooling', 'Night-time cooling', 'Evapotranspiration', 'Albedo (brightening, x 0.01)',
        'Biomass carbon', 'Soil carbon', 'Net climate benefit']


def es(v):
    m = v.notna() & w.notna()
    mu = float((v[m] * w[m]).sum() / w[m].sum())
    sd = float(np.sqrt(((v[m] - mu) ** 2 * w[m]).sum() / w[m].sum()))
    return mu / sd, mu


rows = []
for k, name in (('N', 'Natural forest vs cropland'), ('P', 'Plantation vs cropland'),
                ('D', 'Natural forest replaced by plantation')):
    bio = df[f'e{k}_AGB'] * (1 + RS) * CF * C2CO2
    socc = soc[k] * C2CO2
    alb = df['SWdown'] * df[f'e{k}_ALB'] * TA * K
    series = [-df[f'e{k}_LSTd'], -df[f'e{k}_LSTn'], df[f'e{k}_ET'], 100 * df[f'e{k}_ALB'], bio, socc, alb + bio + socc]
    for ax, sv in zip(AXES, series):
        e, mu = es(sv)
        rows.append(dict(transition=name, effect=ax, standardized_effect=e, mean=mu))
pd.DataFrame(rows).to_csv(outf, index=False)
