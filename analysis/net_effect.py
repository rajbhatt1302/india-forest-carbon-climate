"""Net climate effect (t CO2-eq/ha, 100 yr): albedo forcing as CO2 equivalent plus biomass and soil carbon.
RF = -SWdown * dAlbedo * 0.854; AGWP_CO2,100 = 8.95e-14 (IPCC AR6); root:shoot by aridity (IPCC 2019).
Positive = climate benefit.
Usage: python net_effect.py <window_csv> <predictor_csv> <soc_csv> <out_csv>"""
import sys
import numpy as np
import pandas as pd
import results as R

TA, TH, AE, AGWP = 0.854, 100.0, 5.1e14, 8.95e-14
K_EESF = TH / (AE * AGWP) * 1e4 / 1e3
CF, C2CO2 = 0.47, 44.0 / 12.0
NAME = {'N': 'Natural forest', 'P': 'Plantation', 'D': 'Replacement'}

csv, pcsv, scsv, out = sys.argv[1:5]
df = R.load(csv)
p = pd.read_csv(pcsv).drop(columns=['system:index', '.geo'], errors='ignore')
p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
df = df.merge(p[['gid', 'SWdown', 'AI']], on='gid', how='left')
sc = pd.read_csv(scsv).replace(-9999, np.nan)
df = df.merge(sc[['gid', 'n_r4', 'bN_SOC', 'bP_SOC']], on='gid', how='left')
okn = df['n_r4'] >= R.MIN_N
df['eN_SOC'] = df['bN_SOC'].where(okn & (df['nNat'] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS))
df['eP_SOC'] = df['bP_SOC'].where(okn & (df['nPlant'] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS))
df['eD_SOC'] = (df['bP_SOC'] - df['bN_SOC']).where(okn & (df['nNat'] >= R.MIN_CLASS) & (df['nPlant'] >= R.MIN_CLASS))
RS = np.where(df['AI'] > 1.0, 0.207, np.where(df['AI'] >= 0.65, 0.323, 0.440))
w = R.weights(df)
for k in 'NPD':
    rf = -df['SWdown'] * df[f'e{k}_ALB'] * TA
    alb = -rf * K_EESF
    carb = df[f'e{k}_AGB'] * (1 + RS) * CF * C2CO2 + df[f'e{k}_SOC'] * C2CO2
    ok = alb.notna() & carb.notna()
    df[f'rf{k}'] = rf.where(ok); df[f'alb{k}'] = alb.where(ok); df[f'carb{k}'] = carb.where(ok)
    df[f'net{k}'] = (alb + carb).where(ok)

rows = []
for k in 'NPD':
    for v, lab in (('alb', 'Albedo'), ('carb', 'Carbon'), ('net', 'Net')):
        rows.append(dict(transition=NAME[k], component=lab, **R.weighted_summary(df[f'{v}{k}'], w)))
    rows.append(dict(transition=NAME[k], component='RF_Wm2', **R.weighted_summary(df[f'rf{k}'], w)))

# albedo offset of the carbon benefit by aridity class
CL = [(0.2, 0.5, 'Semiarid'), (0.5, 0.65, 'Dry subhumid'), (0.65, 1.0, 'Humid'), (1.0, 99, 'Very humid')]
for k in 'NP':
    for a, b, lab in CL:
        m = (df['AI'] >= a) & (df['AI'] < b)
        sa = R.weighted_summary(df.loc[m, f'alb{k}'], w[m]); sc = R.weighted_summary(df.loc[m, f'carb{k}'], w[m])
        if sa['n'] >= 8 and sc['mean'] > 0:
            mm = m & df[f'alb{k}'].notna()
            v_a, v_c, ww = df.loc[mm, f'alb{k}'].to_numpy(), df.loc[mm, f'carb{k}'].to_numpy(), w[mm].to_numpy()
            bl = R._BLOCK.reindex(df.index[mm]).to_numpy()
            rng = np.random.default_rng(1); ub = np.unique(bl); pos = {b: np.where(bl == b)[0] for b in ub}
            boot = []
            for _ in range(2000):
                idx = np.concatenate([pos[b] for b in rng.choice(ub, len(ub), replace=True)])
                boot.append(-(v_a[idx] * ww[idx]).sum() / (v_c[idx] * ww[idx]).sum() * 100)
            r = -(v_a * ww).sum() / (v_c * ww).sum() * 100
            rows.append(dict(transition=NAME[k], component='offset_pct_' + lab, n=len(v_a), mean=r,
                             lo=np.percentile(boot, 2.5), hi=np.percentile(boot, 97.5)))
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False))
