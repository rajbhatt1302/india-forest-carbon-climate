"""SAFE-CP gap filling: tree committee (+ TabPFN), residual kriging, conformal 90% intervals.
Options: SAFECP_KRIGE (nested/always/never), SAFECP_TYPE=<csv>, SAFECP_CHAIN, SAFECP_NO_TABPFN, SAFECP_TAG.
Usage: python safecp.py <window_csv> <predictor_csv> [extra_csv] <out_dir>"""
import sys
import os
import json
import time
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
from sklearn.model_selection import GroupKFold, KFold
import lightgbm as lgb
import xgboost as xgb
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
import results as R

SEED = 7
ALPHA = 0.10
KM = 111.32
args = sys.argv[1:]
win, prd = args[0], args[1]
extra = args[2] if len(args) == 4 else None
OUT = args[-1]
os.makedirs(OUT, exist_ok=True)
def make_tabpfn():
    # TabPFN v2 regressor
    from tabpfn import TabPFNRegressor
    from tabpfn.constants import ModelVersion
    m = TabPFNRegressor.create_default_for_version(ModelVersion.V2, device='cpu')
    m.set_params(random_state=SEED, n_estimators=4)
    return m


try:
    _t = make_tabpfn()
    _t.fit(np.random.rand(30, 3), np.random.rand(30)); _t.predict(np.random.rand(5, 3))
    HAS_TABPFN = os.environ.get('SAFECP_NO_TABPFN') is None
except Exception as e:
    print('TabPFN unavailable:', str(e)[:200])
    HAS_TABPFN = False
print('TabPFN available:', HAS_TABPFN)

# data
df = R.load(win)
p = pd.read_csv(prd).drop(columns=['system:index', '.geo'], errors='ignore')
p.loc[p['SWdown'] < 150, 'SWdown'] = np.nan
p['logAI'] = np.log10(p['AI'].clip(0.01, 20))
df = df.merge(p, on='gid', how='left')
PRED = ['logAI', 'P', 'PET', 'Tmax', 'Tmin', 'VPD', 'Pjjas', 'SWdown', 'elev', 'elevSD', 'slope',
        'LSTd_bg', 'LSTn_bg', 'ALB_bg', 'ET_bg']
if extra:
    e = pd.read_csv(extra).drop(columns=['system:index', '.geo'], errors='ignore')
    df = df.merge(e, on='gid', how='left')
    PRED += ['NDVI_rabi', 'NDVI_dry', 'SM', 'WIND']                  # available in every window
TYPE = os.environ.get('SAFECP_TYPE')
if TYPE:
    ty = pd.read_csv(TYPE).drop(columns=['system:index', '.geo'], errors='ignore')
    for v in ('EVG', 'FSEAS', 'CNC', 'CGSL', 'CAMP', 'NDVI_kharif'):
        ty[v + '_reg'] = ty[v + '_R100'].fillna(ty[v + '_R250'])
    ty['EVG_reg'] = ty['EVG_reg'].clip(0, 1)
    TPRED = ['EVG_reg', 'FSEAS_reg', 'CNC_reg', 'CGSL_reg', 'CAMP_reg', 'NDVI_kharif_reg', 'CLAY', 'SAND']
    df = df.merge(ty[['gid'] + TPRED], on='gid', how='left')
    PRED += TPRED
if 'bN_SOC' not in df.columns:                      # SOC regression is a separate table
    SOC_CSV = os.path.join(ROOT, 'data/gee_outputs/grid05_soc.csv')
    sc = pd.read_csv(SOC_CSV).replace(-9999, np.nan)
    df = df.merge(sc[['gid', 'n_r4', 'bN_SOC', 'bP_SOC']], on='gid', how='left')
RS = np.where(df['AI'] > 1.0, 0.207, np.where(df['AI'] >= 0.65, 0.323, 0.440))
K = 100.0 / (5.1e14 * 8.95e-14) * 1e4 / 1e3
okn = df['n_r4'] >= R.MIN_N
for k, cls in (('N', 'nNat'), ('P', 'nPlant')):
    soc = df[f'b{k}_SOC'].where(okn & (df[cls] >= R.MIN_CLASS) & (df['nCrop'] >= R.MIN_CLASS))
    df[f'net{k}'] = df['SWdown'] * df[f'e{k}_ALB'] * 0.854 * K + df[f'e{k}_AGB'] * (1 + RS) * 0.47 * 44 / 12 + soc * 44 / 12
df['x_km'] = df['lon'] * KM * np.cos(np.deg2rad(22.0))
df['y_km'] = df['lat'] * KM


def blocks(deg):
    return (np.floor(df['lon'] / deg).astype(int).astype(str) + '_' + np.floor(df['lat'] / deg).astype(int).astype(str)).to_numpy()


# learners
def committee():
    return [RandomForestRegressor(n_estimators=300, min_samples_leaf=5, max_features=0.5, n_jobs=2, random_state=SEED),
            ExtraTreesRegressor(n_estimators=300, min_samples_leaf=5, max_features=0.6, n_jobs=2, random_state=SEED),
            lgb.LGBMRegressor(n_estimators=400, learning_rate=0.03, num_leaves=8, min_child_samples=15, subsample=0.8,
                              subsample_freq=1, colsample_bytree=0.7, reg_lambda=1.0, random_state=SEED, verbose=-1),
            xgb.XGBRegressor(n_estimators=400, learning_rate=0.03, max_depth=3, min_child_weight=5, subsample=0.8,
                             colsample_bytree=0.7, reg_lambda=1.0, random_state=SEED, n_jobs=2)]


class Base:
    """Committee of tree ensembles, optionally averaged with TabPFN."""
    def __init__(self, use_tabpfn):
        self.use_tabpfn = use_tabpfn

    def fit(self, X, y):
        self.ms = [m.fit(X, y) for m in committee()]
        if self.use_tabpfn:
            self.tp = make_tabpfn().fit(X, y)
        return self

    def predict(self, X, parts=False):
        c = np.mean([m.predict(X) for m in self.ms], axis=0)
        t = self.tp.predict(X) if self.use_tabpfn else None
        out = c if t is None else 0.5 * (c + t)
        return (out, c, t) if parts else out


def kriger():
    k = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=300.0, length_scale_bounds=(20.0, 5000.0), nu=1.5) \
        + WhiteKernel(1.0, (1e-4, 1e3))
    return GaussianProcessRegressor(kernel=k, normalize_y=True, n_restarts_optimizer=0, random_state=SEED)


def fit_srk(X, y, xy, grp, use_tabpfn, inner_k=4):
    """Base model plus kriging of its out-of-fold residuals, if a nested test says it helps."""
    folds = list(GroupKFold(n_splits=min(inner_k, len(np.unique(grp)))).split(X, y, grp))
    oof = np.zeros(len(y))
    for tr, te in folds:
        oof[te] = Base(use_tabpfn).fit(X[tr], y[tr]).predict(X[te])
    res = y - oof
    kr = np.zeros(len(y))
    for tr, te in folds:
        kr[te] = kriger().fit(xy[tr], res[tr]).predict(xy[te])
    gp = kriger().fit(xy, res)
    mode = os.environ.get('SAFECP_KRIGE', 'nested')
    gp.use_kriging = {'always': True, 'never': False}.get(mode, bool(np.sum((res - kr) ** 2) < np.sum(res ** 2)))
    base = Base(use_tabpfn).fit(X, y)
    return base, gp


def predict_srk(base, gp, X, xy):
    b, c, t = base.predict(X, parts=True)
    r, s = gp.predict(xy, return_std=True)
    return (b + r if gp.use_kriging else b), s, b, c, t


CHAIN = {'eN_LSTd': 'eN_ET', 'eN_LSTn': 'eN_ET', 'eP_LSTd': 'eP_ET'} if os.environ.get('SAFECP_CHAIN', '1') == '1' else {}


def chain_feature(Xtr, yc, grp_tr, Xte, k=4):
    """Out-of-fold ET predictions for training rows and fitted predictions for new rows."""
    m = ~np.isnan(yc)
    oof = np.full(len(Xtr), np.nan)
    for a, b in GroupKFold(n_splits=min(k, len(np.unique(grp_tr)))).split(Xtr, yc, grp_tr):
        a = a[m[a]]
        oof[b] = Base(False).fit(Xtr[a], yc[a]).predict(Xtr[b])
    full = Base(False).fit(Xtr[m], yc[m])
    return oof, full.predict(Xte)


def r2(y, p):
    return float(1 - np.sum((y - p) ** 2) / np.sum((y - np.mean(y)) ** 2))


# run
TARGETS = os.environ.get('SAFECP_TARGETS', 'eN_LSTd,eP_LSTd,eN_LSTn,eN_ALB,eN_ET,eP_ET,eN_AGB,netN,netP').split(',')
skill, cvrows, potcols = [], [], {}
TAG = os.environ.get('SAFECP_TAG', 'extra' if extra else 'base')
BDM, BST = float(os.environ.get('SAFECP_BLOCK_DM', 3.0)), float(os.environ.get('SAFECP_BLOCK_ST', 6.0))
print('krige', os.environ.get('SAFECP_KRIGE', 'nested'), 'tag', TAG, 'predictors', len(PRED), 'chain', CHAIN, 'blocks', BDM, BST, flush=True)
PARTS = f'{OUT}/safecp_parts_{TAG}'
os.makedirs(PARTS, exist_ok=True)
for t in TARGETS:
    if os.path.exists(f'{PARTS}/{t}_skill.json'):
        print('done already', t, flush=True)
        continue
    t0 = time.time()
    ok = (df[t].notna() & df[PRED].notna().all(1)).to_numpy()
    X, y = df.loc[ok, PRED].to_numpy(float), df.loc[ok, t].to_numpy(float)
    xy = df.loc[ok, ['x_km', 'y_km']].to_numpy()
    rec = dict(target=t, n=int(ok.sum()), tabpfn=HAS_TABPFN, chain=CHAIN.get(t, ''), n_pred=len(PRED) + (t in CHAIN))
    yc = df.loc[ok, CHAIN[t]].to_numpy(float) if t in CHAIN else None
    for deg, tag in ((BDM, 'dm'), (BST, 'stress')):
        grp = blocks(deg)[ok]
        final = np.zeros(len(y)); sd = np.zeros(len(y)); basep = np.zeros(len(y)); comm = np.zeros(len(y)); kuse = []
        tabp = np.full(len(y), np.nan); fold = np.zeros(len(y), int)
        use_tp = HAS_TABPFN and tag == 'dm'          # committee only for the stress test
        for f, (tr, te) in enumerate(GroupKFold(n_splits=10).split(X, y, grp)):
            Xtr, Xte = X[tr], X[te]
            if yc is not None:
                ctr, cte = chain_feature(X[tr], yc[tr], grp[tr], X[te])
                Xtr, Xte = np.column_stack([Xtr, ctr]), np.column_stack([Xte, cte])
            base, gp = fit_srk(Xtr, y[tr], xy[tr], grp[tr], use_tp, inner_k=3 if use_tp else 4)
            final[te], sd[te], basep[te], comm[te], tt = predict_srk(base, gp, Xte, xy[te])
            kuse.append(gp.use_kriging)
            if tt is not None:
                tabp[te] = tt
            fold[te] = f
        rec[f'r2_committee_{tag}'] = r2(y, comm)
        if use_tp:
            rec[f'r2_tabpfn_{tag}'] = r2(y, tabp)
            rec[f'r2_avg_{tag}'] = r2(y, basep)
        rec[f'r2_final_{tag}'] = r2(y, final)
        rec[f'kriging_used_share_{tag}'] = float(np.mean(kuse))
        rec[f'rmse_final_{tag}'] = float(np.sqrt(np.mean((y - final) ** 2)))
        # cross-conformal coverage: calibrate on the other folds, evaluate on the held-out fold
        score = np.abs(y - final) / sd
        cover, width = [], []
        for f in range(10):
            nf = (fold != f).sum()
            q = np.quantile(score[fold != f], min(1.0, np.ceil((1 - ALPHA) * (nf + 1)) / nf))
            lo, hi = final[fold == f] - q * sd[fold == f], final[fold == f] + q * sd[fold == f]
            cover.append(np.mean((y[fold == f] >= lo) & (y[fold == f] <= hi))); width.append(np.mean(hi - lo))
        rec[f'coverage90_{tag}'] = float(np.mean(cover)); rec[f'width90_{tag}'] = float(np.mean(width))
        if tag == 'dm':
            qdm = float(np.quantile(score, min(1.0, np.ceil((1 - ALPHA) * (len(score) + 1)) / len(score))))
            cvrows.append(pd.DataFrame({'target': t, 'gid': df.loc[ok, 'gid'].to_numpy(), 'obs': y, 'final': final,
                                        'committee': comm, 'tabpfn': tabp, 'sd': sd}))
    # final fit on all observed windows -> gap windows (inside the predictor envelope)
    grp = blocks(BDM)[ok]
    todo = (df[t].isna() & df[PRED].notna().all(1) & (df['n_r1'] >= 0)).to_numpy()
    Xn, xyn = df.loc[todo, PRED].to_numpy(float), df.loc[todo, ['x_km', 'y_km']].to_numpy()
    XF, XnF, FEAT = X, Xn, list(PRED)
    if yc is not None:
        ctr, cte = chain_feature(X, yc, grp, Xn)
        XF, XnF, FEAT = np.column_stack([X, ctr]), np.column_stack([Xn, cte]), FEAT + ['ET_effect_pred']
    base, gp = fit_srk(XF, y, xy, grp, HAS_TABPFN, inner_k=3 if HAS_TABPFN else 4)
    # area of applicability in predictor space (Meyer and Pebesma 2021)
    mu, s_ = X.mean(0), X.std(0) + 1e-9
    Xs, Xns = (X - mu) / s_, (Xn - mu) / s_
    from scipy.spatial.distance import cdist, pdist
    mean_d = pdist(Xs).mean()
    cvdi = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=10).split(X, y, grp):
        cvdi[te] = cdist(Xs[te], Xs[tr]).min(1) / mean_d
    q1, q3 = np.percentile(cvdi, [25, 75]); thr = min(q3 + 1.5 * (q3 - q1), cvdi.max())
    din = cdist(Xns, Xs).min(1) / mean_d
    pr, sdn, _, _, _ = predict_srk(base, gp, XnF, xyn)
    inside = din <= thr
    for nm, val in (('m', pr), ('lo', pr - qdm * sdn), ('hi', pr + qdm * sdn)):
        col = pd.Series(np.nan, index=df.index)
        col.loc[np.where(todo)[0][inside]] = val[inside]
        potcols[f'{nm}_{t}'] = col
    # TreeSHAP averaged over the committee (TabPFN not included)
    if t in ('eN_LSTd', 'eP_LSTd', 'eN_ET', 'eN_AGB', 'netN'):
        import shap
        sv = np.mean([shap.TreeExplainer(m).shap_values(XF) for m in base.ms], axis=0)
        np.save(f'{OUT}/shap_{t}_{TAG}.npy', sv)
        pd.DataFrame(XF, columns=FEAT).assign(gid=df.loc[ok, 'gid'].to_numpy()).to_csv(f'{OUT}/shapX_{t}_{TAG}.csv', index=False)
    rec['kriging_used_final'] = gp.use_kriging
    rec.update(q90_dm=qdm, aoa_thr=float(thr), n_potential=int(inside.sum()), n_outside_aoa=int((~inside).sum()),
               minutes=round((time.time() - t0) / 60, 1))
    pd.concat([df[['gid']], pd.DataFrame({k: v for k, v in potcols.items() if k.endswith('_' + t)})], axis=1) \
        .to_csv(f'{PARTS}/{t}_potential.csv', index=False)
    cvrows[-1].to_csv(f'{PARTS}/{t}_cv.csv', index=False)
    json.dump({k: (v if not isinstance(v, (np.floating, np.integer, np.bool_)) else v.item()) for k, v in rec.items()},
              open(f'{PARTS}/{t}_skill.json', 'w'))
    print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items()}), flush=True)

# combine all finished targets
tag = TAG
done = [t for t in TARGETS if os.path.exists(f'{PARTS}/{t}_skill.json')]
pd.DataFrame([json.load(open(f'{PARTS}/{t}_skill.json')) for t in done]).to_csv(f'{OUT}/safecp_skill_{tag}.csv', index=False)
pd.concat([pd.read_csv(f'{PARTS}/{t}_cv.csv') for t in done]).to_csv(f'{OUT}/safecp_cv_{tag}.csv', index=False)
pot = df[['gid', 'lon', 'lat']].copy()
for t in done:
    pot = pot.merge(pd.read_csv(f'{PARTS}/{t}_potential.csv'), on='gid', how='left')
pot.to_csv(f'{OUT}/safecp_potential_{tag}.csv', index=False)
json.dump(PRED, open(f'{OUT}/safecp_predictors_{tag}.json', 'w'))
print('targets done:', done, 'of', TARGETS)
print('saved')
