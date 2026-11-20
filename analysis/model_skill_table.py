"""Cross-validated skill of SAFE-CP and its variants (3 deg blocks).
Usage: python model_skill_table.py <data_dir> <out_csv>"""
import sys
import numpy as np
import pandas as pd

dd, out = sys.argv[1:3]
rd = lambda f: pd.read_csv(f'{dd}/{f}').set_index('target')
A, C, T = rd('ml_ablation/safecp_skill_ablA.csv'), rd('ml_ablation/safecp_skill_ablC.csv'), rd('ml_type/safecp_skill_typeC.csv')
F, G = rd('ml_finalnest/safecp_skill_finalnest.csv'), rd('ml_g025/safecp_skill_g025base.csv')
TG = [('eN_LSTd', 'Daytime LST', 'N'), ('eP_LSTd', 'Daytime LST', 'P'), ('eN_LSTn', 'Night-time LST', 'N'),
      ('eN_ALB', 'Albedo', 'N'), ('eN_ET', 'ET', 'N'), ('eP_ET', 'ET', 'P'), ('eN_AGB', 'Biomass', 'N'),
      ('netN', 'Net climate effect', 'N'), ('netP', 'Net climate effect', 'P')]
GN = {'N': 'natural forest', 'P': 'plantation'}
COLS = [('Residual kriging', lambda t: A.loc[t, 'r2_final_dm']),
        ('Chained ET', lambda t: C.loc[t, 'r2_final_dm'] if t in C.index else np.nan),
        ('Extra predictors', lambda t: T.loc[t, 'r2_committee_dm']),
        ('Final model', lambda t: F.loc[t, 'r2_final_dm']),
        ('0.25° grid', lambda t: G.loc[t, 'r2_committee_dm'])]
rows = []
for t, nm, g in TG:
    r = {'Target': f'{nm}, {GN[g]}', 'Base model': A.loc[t, 'r2_committee_dm']}
    for lab, f in COLS:
        r[lab] = f(t)
    r['Final: RMSE'] = F.loc[t, 'rmse_final_dm']; r['Final: 90% interval coverage'] = F.loc[t, 'coverage90_dm']
    r['Windows (0.5°)'] = int(F.loc[t, 'n'])
    rows.append(r)
pd.DataFrame(rows).round(3).to_csv(out, index=False)
