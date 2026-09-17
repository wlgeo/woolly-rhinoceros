"""
d18O_water.py -- reconstruct d18O_water (VSMOW) from measured d18O_CO2
(90 degC acid) and clumped-isotope T, with error propagation via ufloats.

Carbonate d18O is reported on both VSMOW and VPDB.  Note that only the
*carbonate* value is put on VPDB: d18O of the CO2 itself is kept on VSMOW,
since the "CO2 vs VPDB" column written by D4Xgui applies the carbonate
conversion to a gas and therefore sits on no standard scale.

The carbonate standards are carried through steps 1-2 as a scale check
against Bernasconi et al. 2018 (Table 4: ETH-1 -2.19, ETH-2 -18.69,
ETH-3 -1.78 permil VPDB), but no d18O_water is derived for them -- they are
heated, re-equilibrated or kinetically altered materials, so the
equilibrium calcite-water fractionation of step 2 does not apply.

Step 1: alpha_acid(90 degC) = exp(3.59/T_K - 1.79e-3) = 1.008129
        (Kim, Mucci & Taylor 2007, calcite)
Step 2: 1000 ln alpha_cc-w = A*10^3*(1/T-1/T0) + B
        Kim & O'Neil 1997 (OLS refit of Table 1 5 mM subset via refit_KO97,
            in the Daëron-2019 uncorrelated/centred form; n=9, dof=7):
            A = 18.03 +/- 0.36 , B = 27.81 +/- 0.05 , T0 = 299.4 K
            (equivalent to the paper's 18.03*10^3/T_K - 32.42; r(A,B)_paper
            form = +0.9993, centred here so parameter errors are independent.)

References
    Kim & O'Neil 1997        doi:10.1016/S0016-7037(97)00169-5
    Kim, Mucci & Taylor 2007 doi:10.1016/j.chemgeo.2007.08.005
    Daëron et al. 2019       doi:10.1038/s41467-019-08336-5
                             (centred/uncorrelated form used to refit KO97)
    Fiebig et al. 2024       doi:10.1016/j.chemgeo.2024.122382
"""

import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.path as mpath
from matplotlib.markers import MarkerStyle
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from uncertainties import ufloat, umath, covariance_matrix

warnings.filterwarnings('ignore',
    message='Using UFloat objects with std_dev==0 may give unexpected results.')


PATH = 'Supplementary Dataset 5.xlsx'
SUMMARY_SHEET = 'Summary of clumped isotope '
CONF_LEVEL = 0.95

STANDARDS = ['ETH-1', 'ETH-2', 'ETH3oxi', 'GU1']

SAMPLE_MASK = lambda s: (s.str.contains('17NTS') | s.str.contains('XSC')
                         | s.isin(STANDARDS))
GROUPS = [('17NTS', lambda s: s.str.contains('NTS')),
          ('XSC',   lambda s: s.str.contains('XSC'))]


# calibration constants -- see module docstring for sources
T_ACID_K   = 90.0 + 273.15
ACID_ALPHA = ufloat(round(np.exp(3.59 / T_ACID_K - 1.79e-3), 6), 0.0)  # 1.008129

# d18O of VPDB on the VSMOW scale: d18O_VSMOW = 1.03092*d18O_VPDB + 30.92
# (Kim, Coplen & Horita 2015, as used by Daeron et al. 2019)
D18O_VPDB_VSMOW = 30.92

# Kim & O'Neil 1997 Table 1, 5 mM Ca2+/HCO3- subset (n=9).  Columns are
# the paper's published 10^3 ln alpha_cc-w values (VSMOW); they match
# 1000*ln((1000+d18O_carb)/(1000+d18O_water)) to <0.005 permil.
KO97_TABLE1_5mM = dict(
    sample=np.array(['ST-52-10', 'ST-54-10',
                     'ST-45-25-N', 'ST-50-25', 'ST-52-25', 'ST-54-25',
                     'ST-50-40', 'ST-52-40', 'ST-54-40']),
    T_C=np.array([10., 10., 25., 25., 25., 25., 40., 40., 40.]),
    ln1000=np.array([31.35, 31.21, 27.87, 28.31, 27.96, 28.10,
                     25.15, 25.25, 25.13]),
)


def refit_KO97(table=KO97_TABLE1_5mM):
    """OLS refit of Kim & O'Neil 1997 Table 1 (5 mM subset) in the
    Daëron-2019 centred form

        1000 ln alpha = A * 10^3 * (1/T - 1/T0) + B

    with T0 = 1/mean(1/T) so that A and B are uncorrelated.  Equivalent
    to the paper's 18.03*10^3/T - 32.42 (r(A,B) ~ +0.9993 in that form).
    Returns (A, B, T0) with A, B as ufloats (1SE)."""
    T_K = table['T_C'] + 273.15
    y = np.asarray(table['ln1000'], dtype=float)
    T0 = 1.0 / np.mean(1.0 / T_K)
    x = 1000.0 * (1.0 / T_K - 1.0 / T0)
    X = np.column_stack([x, np.ones(len(y))])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - 2
    cov = np.sum(resid**2) / dof * np.linalg.inv(X.T @ X)
    A = ufloat(round(float(beta[0]), 2), round(float(np.sqrt(cov[0, 0])), 2))
    B = ufloat(round(float(beta[1]), 2), round(float(np.sqrt(cov[1, 1])), 2))
    return A, B, round(float(T0), 1)


KO97_A, KO97_B, KO97_T0 = refit_KO97()   # 18.03+/-0.36, 27.81+/-0.05, 299.4 K


def d18O_carb_from_CO2(d18O_CO2):
    return (1000.0 + d18O_CO2) / ACID_ALPHA - 1000.0


def to_VPDB(d18O_VSMOW):
    return ((d18O_VSMOW - D18O_VPDB_VSMOW)
            / (1.0 + D18O_VPDB_VSMOW / 1000.0))


def d18O_water_KO97(d18O_carb, T_K):
    alpha = umath.exp(
        (KO97_A * 1000.0 * (1.0/T_K - 1.0/KO97_T0) + KO97_B) / 1000.0)
    return (1000.0 + d18O_carb) / alpha - 1000.0


def process():
    df = pd.read_excel(PATH, sheet_name=SUMMARY_SHEET)
    df = df.loc[SAMPLE_MASK(df['Sample name']) & (df['N'] > 1)].copy()

    rows, ufloats = [], {}
    for _, r in df.iterrows():
        is_std = r['Sample name'] in STANDARDS
        T_lo, T_hi, T_mean = (r['T(min, 1SE), Fiebig24 (original)'],
                              r['T(max, 1SE), Fiebig24 (original)'],
                              r['T(mean), Fiebig24 (original)'])
        has_T = not any(pd.isna(x) for x in (T_lo, T_hi, T_mean))
        if not is_std and not has_T:
            continue

        d18O_CO2       = ufloat(r['d18O_CO2_VSMOW'], r['SD_d18O'] / np.sqrt(r['N']))
        d18O_carb      = d18O_carb_from_CO2(d18O_CO2)
        d18O_carb_pdb  = to_VPDB(d18O_carb)

        row = dict(
            Sample=r['Sample name'], N=int(r['N']),
            T_C=np.nan, SE_T_C=np.nan,
            d18O_CO2_VSMOW=round(d18O_CO2.nominal_value, 2),
            SE_d18O_CO2=round(d18O_CO2.std_dev, 2),
            d18O_carb_VSMOW=round(d18O_carb.nominal_value, 2),
            SE_d18O_carb=round(d18O_carb.std_dev, 2),
            d18O_carb_VPDB=round(d18O_carb_pdb.nominal_value, 2),
            SE_d18O_carb_VPDB=round(d18O_carb_pdb.std_dev, 2),
            d18O_water_VSMOW_KO97=np.nan, SE_d18O_water_KO97=np.nan,
        )

        # standards are heated / re-equilibrated / kinetically altered, so
        # the equilibrium calcite-water fractionation does not apply to them
        if not is_std:
            T_K    = ufloat(T_mean + 273.15, 0.5 * (T_hi - T_lo))
            T_C    = T_K - 273.15
            d18O_w = d18O_water_KO97(d18O_carb, T_K)
            row.update(
                T_C=T_C.nominal_value, SE_T_C=T_C.std_dev,
                d18O_water_VSMOW_KO97=round(d18O_w.nominal_value, 2),
                SE_d18O_water_KO97=round(d18O_w.std_dev, 2),
            )
            ufloats[r['Sample name']] = dict(T_C=T_C, d18O_w_KO97=d18O_w)
        rows.append(row)

    result = pd.DataFrame(rows)
    result['_std'] = result['Sample'].isin(STANDARDS)
    result = (result.sort_values('_std', kind='stable')
                    .drop(columns='_std')
                    .reset_index(drop=True))
    return result, ufloats


# --- plot ----------------------------------------------------------------

plt.rc('font', family='Times New Roman', size=8)
plt.rc('axes', titlesize=10, labelsize=12, labelpad=7)
plt.rc('xtick', labelsize=10); plt.rc('ytick', labelsize=10)
plt.rc('legend', fontsize=8); plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,wasysym}'
plt.rcParams['axes.unicode_minus']  = False

_c, _s4, _s5, _s6 = (mpath.Path.unit_circle(),
                     mpath.Path.unit_regular_star(4),
                     mpath.Path.unit_regular_star(5),
                     mpath.Path.unit_regular_star(6))
MARKER = [_s5, 's', '^', MarkerStyle('o', fillstyle='right'),
          '<', _s4, '>', '*', MarkerStyle('o', fillstyle='top'),
          mpath.Path(vertices=np.concatenate([_c.vertices, _s4.vertices[::-1]]),
                     codes=np.concatenate([_c.codes, _s4.codes])),
          'h', _s6, 'd', 'H', '8', MarkerStyle('o', fillstyle='left'),
          'D', 'v', 'P', 'o', '^'] * 10
COLORS = [[v/255 for v in c] for c in [
    (230,159,0,0), (86,180,233,0), (0,158,115,0), (240,228,66,0),
    (0,114,178,0), (213,94,0,0),   (204,121,167,0), (255,0,0,0)]] * 20

_natural_key = lambda s: [int(t) if t.isdigit() else t.lower()
                          for t in re.split(r'(\d+)', str(s)) if t]


def conf_ellipse(x, y, p, **kwargs):
    """Joint p-level confidence ellipse for two (possibly correlated) ufloats,
    from the 2x2 covariance matrix: axes = 2*sqrt(lambda_i * chi2_p(2)),
    rotation from the eigenvectors. Returns a matplotlib Ellipse patch."""
    val, vec = np.linalg.eigh(np.array(covariance_matrix((x, y))))
    width, height = 2.0 * np.sqrt(val * chi2.ppf(p, 2))
    angle = np.degrees(np.arctan2(vec[1, 0], vec[0, 0]))
    return Ellipse(xy=(x.nominal_value, y.nominal_value),
                   width=width, height=height, angle=angle, **kwargs)


def plot(result, ufloats):
    saved = []
    k = float(np.sqrt(chi2.ppf(CONF_LEVEL, df=2)))
    for name, match in GROUPS:
        samples = sorted(result.loc[match(result['Sample']), 'Sample'].unique(),
                         key=_natural_key)
        if not samples:
            continue
        sub = result.loc[match(result['Sample'])]

        fig, ax = plt.subplots(figsize=(6.5, 5))
        for i, s in enumerate(samples):
            u = ufloats[s]; col = mcolors.to_hex(COLORS[i])
            ax.add_patch(conf_ellipse(u['T_C'], u['d18O_w_KO97'], p=CONF_LEVEL,
                                      fc='none', ec=col, lw=0.8, zorder=10))
            ax.plot(u['T_C'].nominal_value, u['d18O_w_KO97'].nominal_value,
                    marker=MARKER[i], color=col, linestyle='',
                    markeredgewidth=0.25, markeredgecolor='k',
                    label=s, zorder=15)
        ax.set_xlabel(r'$T\ (\text{\textdegree C})$')
        ax.set_ylabel(r'$\delta^{18}\text{O}_{\text{water, VSMOW}}\ (\permil)$')
        ax.grid(True, ls=':', lw=0.4, alpha=0.6, zorder=0)

        # axis limits = actual 95% ellipse bounding box + small margin
        xh = k * sub['SE_T_C']
        yh = k * sub['SE_d18O_water_KO97']
        x0, x1 = (sub['T_C'] - xh).min(), (sub['T_C'] + xh).max()
        y0 = (sub['d18O_water_VSMOW_KO97'] - yh).min()
        y1 = (sub['d18O_water_VSMOW_KO97'] + yh).max()
        ax.set_xlim(x0 - 0.05*(x1-x0), x1 + 0.05*(x1-x0))
        ax.set_ylim(y0 - 0.08*(y1-y0), y1 + 0.08*(y1-y0))

        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0),
                  edgecolor='k', framealpha=1, borderaxespad=0,
                  fancybox=False, title=name)

        fig.tight_layout()
        out = f'd18Owater_{name}.pdf'
        fig.savefig(out, bbox_inches='tight')
        saved.append(out)
    return saved


if __name__ == '__main__':
    n = len(KO97_TABLE1_5mM['ln1000'])
    print(f"KO97 refit (Table 1, 5 mM, n={n}, dof={n-2}; Daëron-2019 form): "
          f"A={KO97_A}, B={KO97_B}, T0={KO97_T0} K\n")

    result, ufloats = process()
    pd.set_option('display.width', 220); pd.set_option('display.max_columns', None)
    print(result.round(2).to_string(index=False))

    out = PATH.replace('.xlsx', '_d18Owater.xlsx')
    result.to_excel(out, index=False)
    print(f"\nwrote {out}")
    plot(result, ufloats)
