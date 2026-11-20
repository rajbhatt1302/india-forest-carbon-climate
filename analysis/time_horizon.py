"""Net climate effect at 20, 100 and 500 yr (AGWP_CO2 from IPCC AR6).
Usage: python time_horizon.py <window_csv> <predictor_csv> <soc_csv> <out_csv>"""
import sys
import numpy as np
import results as R
import pandas as pd

TA, AE, CF, C2CO2 = 0.854, 5.1e14, 0.47, 44.0 / 12.0
AGWP = {20: 2.43e-14, 100: 8.95e-14, 500: 3.14e-13}
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
CL = [(0.2, 0.5, 'Semi-arid'), (0.5, 0.65, 'Dry sub-humid'), (0.65, 1.0, 'Humid'), (1.0, 99, 'Very humid')]
NAME = {'N': 'Natural forest vs cropland', 'P': 'Plantation vs cropland', 'D': 'Natural forest replaced by plantation'}
rows = []
for TH, ag in AGWP.items():
    K = TH / (AE * ag) * 1e4 / 1e3
    for k in 'NPD':
        alb = df['SWdown'] * df[f'e{k}_ALB'] * TA * K
        carb = df[f'e{k}_AGB'] * (1 + RS) * CF * C2CO2 + df[f'e{k}_SOC'] * C2CO2
        ok = alb.notna() & carb.notna()
        for comp, v in (('Albedo', alb), ('Carbon', carb), ('Net', alb + carb)):
            s = R.weighted_summary(v.where(ok), w)
            rows.append(dict(TH=TH, K=K, transition=NAME[k], component=comp, **s))
        if k in 'NP':
            for a, b, lab in CL:
                m = ok & (df['AI'] >= a) & (df['AI'] < b)
                va, vc, ww = alb[m].to_numpy(), carb[m].to_numpy(), w[m].to_numpy()
                if len(va) < 8:
                    continue
                bl = R._BLOCK.reindex(df.index[m]).to_numpy(); rng = np.random.default_rng(1); ub = np.unique(bl)
                pos = {q: np.where(bl == q)[0] for q in ub}; boot = []
                for _ in range(2000):
                    idx = np.concatenate([pos[q] for q in rng.choice(ub, len(ub), replace=True)])
                    boot.append(-(va[idx] * ww[idx]).sum() / (vc[idx] * ww[idx]).sum() * 100)
                rows.append(dict(TH=TH, K=K, transition=NAME[k], component='Albedo offset of carbon (%), ' + lab, n=len(va),
                                 mean=-(va * ww).sum() / (vc * ww).sum() * 100, lo=np.percentile(boot, 2.5),
                                 hi=np.percentile(boot, 97.5)))
res = pd.DataFrame(rows)
res.to_csv(out, index=False)
pd.set_option('display.width', 200)
print(res.round(2).to_string())
