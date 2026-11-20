"""Load the window regression table and derive supported effect estimates."""
import numpy as np
import pandas as pd

Y_ANNUAL = ['LSTd', 'LSTn', 'ALB', 'ET']
Y_SEASON = ['LSTd_MAM', 'LSTn_MAM', 'LSTd_JJAS', 'LSTd_ONDJF']
Y_AGB = ['AGB']
REG_OF = {**{y: 'r1' for y in Y_ANNUAL}, **{y: 'r2' for y in Y_SEASON}, 'AGB': 'r3'}

# minimum support per window
MIN_N = 200          # pixels in the regression
MIN_CLASS = 30       # pixels dominated (>50 %) by a class

# OLS standard errors are inflated for spatial dependence of pixels
SE_INFLATION = 3.0


BLOCK_DEG = 6.0      # bootstrap block size
_BLOCK = None        # block label of every window, indexed like the loaded table


def load(path):
    global _BLOCK
    df = pd.read_csv(path)
    df = df.drop(columns=[c for c in df.columns if c in ('system:index', '.geo')], errors='ignore')
    df = df.replace(-9999, np.nan)
    new = {}
    for y in Y_ANNUAL + Y_SEASON + Y_AGB:
        r = REG_OF[y]
        n = df.get('n_' + r)
        ok_n = n >= MIN_N if n is not None else False
        nat_ok = ok_n & (df['nNat'] >= MIN_CLASS) & (df['nCrop'] >= MIN_CLASS)
        pl_ok = ok_n & (df['nPlant'] >= MIN_CLASS) & (df['nCrop'] >= MIN_CLASS)
        d_ok = ok_n & (df['nNat'] >= MIN_CLASS) & (df['nPlant'] >= MIN_CLASS)
        # effect = tree cover minus cropland
        df['eN_' + y] = df['bN_' + y].where(nat_ok)
        df['eP_' + y] = df['bP_' + y].where(pl_ok)
        df['eD_' + y] = (df['bP_' + y] - df['bN_' + y]).where(d_ok)       # natural forest replaced by plantation
        # 0 estimate, 1 cover absent, 2 too few pixels
        for k, ok in (('N', nat_ok), ('P', pl_ok), ('D', d_ok)):
            df['why' + k + '_' + y] = np.where(ok, 0, np.where(ok_n, 1, 2))
        for k, s in (('N', 'sN_'), ('P', 'sP_'), ('D', 'sD_')):
            se = df[s + y] * SE_INFLATION
            df['z' + k + '_' + y] = (df['e' + k + '_' + y] / se)
            df['sig' + k + '_' + y] = df['z' + k + '_' + y].abs() > 1.96
    _BLOCK = (np.floor(df['lon'] / BLOCK_DEG).astype(int).astype(str) + '_' +
              np.floor(df['lat'] / BLOCK_DEG).astype(int).astype(str))
    return df


def weighted_summary(values, weights, n_boot=2000, seed=1):
    """Area-weighted mean with a 6 deg block bootstrap 95% interval."""
    m = values.notna() & weights.notna() & (weights > 0)
    v, w = values[m].to_numpy(), weights[m].to_numpy()
    if len(v) < 5:
        return dict(n=len(v), mean=np.nan, lo=np.nan, hi=np.nan)
    blocks = _BLOCK.reindex(values[m].index).to_numpy()
    assert not pd.isna(blocks).any(), 'values lost the window index; cannot assign spatial blocks'
    rng = np.random.default_rng(seed)
    ub = np.unique(blocks)
    pos = {b: np.where(blocks == b)[0] for b in ub}
    boots = np.empty(n_boot)
    for i in range(n_boot):
        ii = np.concatenate([pos[b] for b in rng.choice(ub, len(ub), replace=True)])
        boots[i] = (v[ii] * w[ii]).sum() / w[ii].sum()
    return dict(n=len(v), mean=float((v * w).sum() / w.sum()), lo=float(np.percentile(boots, 2.5)),
                hi=float(np.percentile(boots, 97.5)))


def weights(df):
    """Area weight of a window: share inside India x cos(latitude)."""
    return df['fIndia'] * np.cos(np.deg2rad(df['lat']))


MIN_R2_FILL = 0.3


def attach_safecp(df, pot_csv, skill_csv, min_r2=MIN_R2_FILL):
    """Add SAFE-CP estimates for targets with cross-validated R2 >= min_r2."""
    pot = pd.read_csv(pot_csv)
    sk = pd.read_csv(skill_csv).set_index('target')['r2_final_dm']
    ren = {}
    for c in pot.columns:
        if c.startswith('m_e') and c[3] in 'NP' and sk.get(c[2:], 0) >= min_r2:
            ren[c] = 'm' + c[3:]
    return df.merge(pot[['gid'] + list(ren)].rename(columns=ren), on='gid', how='left')


