"""Summary of data_vs_algorithm.csv: change in R2 from each choice."""
import pandas as pd
import numpy as np
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

D = os.path.join(ROOT, 'data')
d = pd.read_csv(f'{D}/stats/data_vs_algorithm.csv')
g = lambda p, des, m: d[(d.predictors == p) & (d.design == des) & (d.model == m)].set_index('target')['r2']
F = pd.read_csv(f'{D}/ml_finalnest/safecp_skill_finalnest.csv').set_index('target')
G = pd.read_csv(f'{D}/ml_g025/safecp_skill_g025base.csv').set_index('target')
A = pd.read_csv(f'{D}/ml_ablation/safecp_skill_ablA.csv').set_index('target')
trees = pd.concat([g('base', 'block3', m) for m in ('RF', 'ET', 'LGBM', 'XGB')], axis=1)
out = pd.DataFrame({
    'Validation design (random vs blocks)': g('base', 'random', 'Committee') - g('base', 'block3', 'Committee'),
    'Stricter blocks (6 vs 3 degree)': g('base', 'block6', 'Committee') - g('base', 'block3', 'Committee'),
    'Tree algorithm (best vs worst)': trees.max(1) - trees.min(1),
    'Machine learning vs linear': g('base', 'block3', 'Committee') - g('base', 'block3', 'Ridge'),
    'Ensembling (committee vs best single)': g('base', 'block3', 'Committee') - trees.max(1),
    'Extra predictors (committee)': g('extra', 'block3', 'Committee') - g('base', 'block3', 'Committee'),
    'Extra predictors (linear)': g('extra', 'block3', 'Ridge') - g('base', 'block3', 'Ridge'),
    'TabPFN added (final model)': F['r2_final_dm'] - F['r2_committee_dm'],
    'Window size (0.25 vs 0.5 degree)': G['r2_committee_dm'] - A['r2_committee_dm'],
})
out = out.loc[['eN_LSTd', 'eP_LSTd', 'eN_LSTn', 'eN_ALB', 'eN_ET', 'eP_ET', 'eN_AGB', 'netN', 'netP']]
s = out.T
s['median'] = out.median(); s['min'] = out.min(); s['max'] = out.max()
s.round(3).to_csv(f'{D}/stats/data_vs_algorithm_summary.csv')
pd.set_option('display.width', 250)
print(s.round(3).to_string())
print('\ncommittee R2 random / block3 / block6 (base):')
print(pd.DataFrame({'random': g('base', 'random', 'Committee'), 'block3': g('base', 'block3', 'Committee'),
                    'block6': g('base', 'block6', 'Committee'), 'ridge_block3': g('base', 'block3', 'Ridge')}).round(3))
