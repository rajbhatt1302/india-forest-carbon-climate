"""National summary table (observed windows, block bootstrap 95% intervals)."""
import pandas as pd

st = pd.read_csv('data/stats/stats_national.csv').drop_duplicates(['var', 'group']).set_index(['var', 'group'])
f9 = pd.read_csv('data/stats/net_effect_national.csv')
th = pd.read_csv('data/stats/time_horizon_sensitivity.csv')
hot = pd.read_csv('data/stats/hotdays.csv')
G = {'N': 'Natural forest', 'P': 'Plantation', 'D': 'Replacement'}
F9 = {'N': 'Natural forest', 'P': 'Plantation', 'D': 'Replacement'}
TH = {'N': 'Natural forest vs cropland', 'P': 'Plantation vs cropland', 'D': 'Natural forest replaced by plantation'}


def fmt(m, lo, hi, n, d):
    if d == 0 and abs(m) < 10:
        d = 1
    s = f'{m:.{d}f} ({lo:.{d}f} to {hi:.{d}f})'
    return s.replace('-', '−'), int(n)


ROWS = [('Surface temperature', None), ('Daytime LST, annual (°C)', ('LSTd', 2)), ('Night-time LST, annual (°C)', ('LSTn', 2)),
        ('Daytime LST, pre-monsoon MAM (°C)', ('LSTd_MAM', 2)), ('Daytime LST, monsoon JJAS (°C)', ('LSTd_JJAS', 2)),
        ('Daytime LST, post-monsoon and winter ONDJF (°C)', ('LSTd_ONDJF', 2)),
        ('Night-time LST, pre-monsoon MAM (°C)', ('LSTn_MAM', 2)),
        ('Daytime LST, hottest pre-monsoon days (°C)', ('hot', 'LSTd_hot90')),
        ('Surface properties', None), ('White-sky shortwave albedo', ('ALB', 3)), ('Evapotranspiration (mm yr⁻¹)', ('ET', 0)),
        ('Aboveground biomass (Mg ha⁻¹)', ('AGB', 1)),
        ('Net climate effect (t CO₂-eq ha⁻¹; positive = climate benefit)', None),
        ('Albedo, 100 yr', ('f9', 'Albedo')), ('Ecosystem carbon (biomass, roots, soil 0 to 30 cm)', ('f9', 'Carbon')),
        ('Net, 20 yr', ('th', 20)), ('Net, 100 yr', ('f9', 'Net')), ('Net, 500 yr', ('th', 500))]
out = []
for lab, spec in ROWS:
    row = {'Variable': lab}
    if spec is None:
        out.append(row); continue
    for k in 'NPD':
        if spec[0] == 'f9':
            r = f9[(f9['transition'] == F9[k]) & (f9['component'] == spec[1])].iloc[0]; d = 0
        elif spec[0] == 'th':
            r = th[(th['TH'] == spec[1]) & (th['transition'] == TH[k]) & (th['component'] == 'Net')].iloc[0]; d = 0
        elif spec[0] == 'hot':
            r = hot[(hot['transition'] == TH[k]) & (hot['quantity'] == spec[1])].iloc[0]; d = 2
        else:
            r = st.loc[(spec[0], k)].rename({'ci_block_lo': 'lo', 'ci_block_hi': 'hi'}); d = spec[1]
        s, n = fmt(r['mean'], r['lo'], r['hi'], r['n'], d)
        row[G[k]] = s; row[G[k] + ' n'] = n
    out.append(row)
t = pd.DataFrame(out)
for c in [c for c in t if c.endswith(' n')]:
    t[c] = t[c].astype('Int64')
t.to_csv('data/stats/national_summary.csv', index=False)
print(t.astype(object).where(t.notna(), '').to_string(index=False))
