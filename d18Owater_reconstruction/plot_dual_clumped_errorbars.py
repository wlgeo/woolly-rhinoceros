"""
plot_dual_clumped_errorbars.py
------------------------------

D47 vs D48 figure with per-sample uncertainties shown as ±2SE x/y error bars,
relative to the Fiebig-24 equilibrium calibration curve and its 95% CI band.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.markers import MarkerStyle
import matplotlib.path as mpath


def _natural_key(name):
	"""Sort key that interprets embedded integer runs as numbers, so
	`17NTS05oxi` < `17NTS20oxi` < `17NTS111oxi` instead of the default
	lexical order `17NTS05oxi` < `17NTS111oxi` < `17NTS20oxi`."""
	return [int(tok) if tok.isdigit() else tok.lower()
	        for tok in re.split(r'(\d+)', str(name)) if tok]

circle = mpath.Path.unit_circle()
star6 = mpath.Path.unit_regular_star(6)
star4 = mpath.Path.unit_regular_star(4)
star10 = mpath.Path.unit_regular_star(5)
MARKER = [
	star10, 's', '^', MarkerStyle("o", fillstyle="right"),
	'<', star4, '>', '*', MarkerStyle("o", fillstyle="top"),
	mpath.Path(vertices=np.concatenate([circle.vertices, star4.vertices[::-1, ...]]),
	           codes=np.concatenate([circle.codes, star4.codes])),
	'h', star6, 'd', 'H', '8', MarkerStyle("o", fillstyle="left"),
	'D', 'v', 'P', 'o', '^'
] * 10

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc("font", size=SMALL_SIZE)
plt.rc("axes", titlesize=MEDIUM_SIZE)
plt.rc("axes", labelsize=BIGGER_SIZE, labelpad=7)
plt.rc("xtick", labelsize=MEDIUM_SIZE)
plt.rc("ytick", labelsize=MEDIUM_SIZE)
plt.rc("legend", fontsize=SMALL_SIZE)
plt.rc("figure", titlesize=BIGGER_SIZE)
plt.rcParams["font.family"] = "Times New Roman"
plt.rc("text", usetex=True)
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath,wasysym}"


PATH = 'Supplementary Dataset 5.xlsx'
SUMMARY_SHEET = 'Summary of clumped isotope '
REPLICATES_SHEET = 'Replicates of all analysis'


SHOW_CALIB_CI = True        # 95 % CI band of the equilibrium curve
CALIB_CI_LEVEL = 0.95


COLORS = [
	[_/255 for _ in (230, 159, 0, 0)],
	[_/255 for _ in (86, 180, 233, 0)],
	[_/255 for _ in (0, 158, 115, 0)],
	[_/255 for _ in (240, 228, 66, 0)],
	[_/255 for _ in (0, 114, 178, 0)],
	[_/255 for _ in (213, 94, 0, 0)],
	[_/255 for _ in (204, 121, 167, 0)],
	[_/255 for _ in (255, 0, 0, 0)],  # red
] * 20


SAMPLE_MASK = lambda s: (s.str.contains('XSC') | s.str.contains('NTS'))
DROP_MASK = lambda s: (s.str.contains('760') | s.str.contains('3150') | s.str.contains('3450'))


def _load_data():
	summary = pd.read_excel(PATH, sheet_name=SUMMARY_SHEET)
	reps = pd.read_excel(PATH, sheet_name=REPLICATES_SHEET)

	m = SAMPLE_MASK(summary['Sample name']) & ~DROP_MASK(summary['Sample name'])
	summary = summary.loc[m].copy()

	m = SAMPLE_MASK(reps['Sample']) & ~DROP_MASK(reps['Sample'])
	reps = reps.loc[m].copy()

	return summary, reps


def scatterSample(summary, ax):
	"""Draw samples in `summary` on `ax` with 2SE x/y error bars."""
	for idx, sample in enumerate(sorted(summary['Sample name'].unique(), key=_natural_key)):
		df_sample = summary.loc[summary['Sample name'] == sample]
		x_mean = df_sample['D48'].mean()
		y_mean = df_sample['D47'].mean()
		color = mcolors.to_hex(COLORS[idx])

		ax.errorbar(
			x_mean, y_mean,
			xerr=df_sample['2SE_D48'].mean(),
			yerr=df_sample['2SE_D47'].mean(),
			markeredgewidth=0.25,
			markeredgecolor='k',
			color=color,
			lw=0,
			capsize=4,
			label=sample,
			marker=MARKER[idx],
			zorder=15,
		)

		LINES = ax.errorbar(
			x_mean, y_mean,
			xerr=df_sample['2SE_D48'].mean(),
			yerr=df_sample['2SE_D47'].mean(),
			linewidth=.15,
			marker=None,
			color='k',
			fmt='--',
			capsize=0.01,
			zorder=13,
		)
		for L in LINES[-1]:
			L.set_linestyle('--')


## Fiebig 2024 dual-clumped calibration in the compact "D63/D64 * scaling +
## offset" form. For each panel:
##
##     D4x(T) = b0 + b1 * H(1/T_K)                                  (1)
##
## with H(x) the theoretical Hill (2014) polynomial (fixed, no uncertainty)
## and (b0, b1) the Fiebig-24 empirical affine parameters (rounded to the same
## 3–4 sig figs used in D4Xgui). Below we build the 5-coefficient
## representation expected by `D95eq.D4x_calib_function` --
##
##     coefs = [b0, b1*h1, b1*h2, b1*h3, b1*h4]                    (2)
##
## with (b0, b1) drawn as *correlated* ufloats via `uncertainties`, then let
## D95eq's `confidence_band` produce the joint 95% envelope by
## taking the union of joint (Δ47, Δ48) confidence ellipses along the
## parametric curve. That gives ONE closed region -- not two overlapping
## rectilinear bands -- which is the proper joint 95% CI in dual-clumped
## space.
##
## The (b0, b1) covariance matrices below were recovered from the Fiebig
## (2024) 95% prediction bands shipped with D95thermo
## (`calib_coefs/D{47,48}calib_Fiebig2024.xlsx`) using the same lstsq recipe
## as `d95_integration._fit_affine_cov_from_band`.

import uncertainties as _uc
import correldata as _cd
import D95eq as _D95eq

## Hill (2014) theoretical polynomials, 4 sig figs (matching D4Xgui).
_HILL_D63 = (-5.897, -3521.0, 2.391e7, -3.541e9)   # -> D47
_HILL_D64 = ( 6.002, -12990.0, 8.996e6, -7.423e8)  # -> D48

## Fiebig 2024 (b0, b1) nominal values and their 2x2 covariance matrices
## (D47: corr = -0.712 ; D48: corr = -0.664).
F24_D47_B0, F24_D47_B1 = 0.1848, 1.038
F24_D47_COV = [
	[ 4.102e-06, -1.212e-05],
	[-1.212e-05,  7.067e-05],
]
F24_D48_B0, F24_D48_B1 = 0.1214, 1.038
F24_D48_COV = [
	[ 7.397e-06, -6.905e-05],
	[-6.905e-05,  1.460e-03],
]


def _build_F24_calib_coefs(b0_nv, b1_nv, cov, hill):
	"""Build the correldata uarray of coefficients expected by
	`D95eq.D4x_calib_function` -- see equation (2) above."""
	b0, b1 = _uc.correlated_values([b0_nv, b1_nv], cov)
	return _cd.uarray([b0, b1 * hill[0], b1 * hill[1], b1 * hill[2], b1 * hill[3]])


F24_D47_COEFS = _build_F24_calib_coefs(F24_D47_B0, F24_D47_B1, F24_D47_COV, _HILL_D63)
F24_D48_COEFS = _build_F24_calib_coefs(F24_D48_B0, F24_D48_B1, F24_D48_COV, _HILL_D64)


def _D47_calib(T):
	return _D95eq.D4x_calib_function(T, F24_D47_COEFS)


def _D48_calib(T):
	return _D95eq.D4x_calib_function(T, F24_D48_COEFS)


def _Fiebig24D47D48(temps_C):
	"""Nominal-only evaluator (for °C markers)."""
	T = np.asarray(temps_C, dtype=float)
	D47_eq = _D95eq.D4x_calib_function(T, F24_D47_COEFS, return_without_uncertainties=True)
	D48_eq = _D95eq.D4x_calib_function(T, F24_D48_COEFS, return_without_uncertainties=True)
	return D47_eq, D48_eq


def plot_calib(ax, show_ci=SHOW_CALIB_CI, ci_level=CALIB_CI_LEVEL):
	annotation = [5, 10, 15, 20, 25]
	ann_eq47, ann_eq48 = _Fiebig24D47D48(np.array(annotation))
	OFFSETS = [0.0035, 0.0000022]

	eq_range = np.arange(3, 28, .1)
	eq47_nom, eq48_nom = _Fiebig24D47D48(eq_range)

	if show_ci:
		## Joint 95% confidence region in (Δ48, Δ47) space, built by D95eq
		## as the *union* of joint 2D confidence ellipses along the curve.
		band_xy = _D95eq.confidence_band(
			eq_range,
			fx=_D48_calib,  # x axis is Δ48
			fy=_D47_calib,  # y axis is Δ47
			p=ci_level,
		)
		ax.add_patch(
			MplPolygon(
				band_xy,
				closed=True,
				facecolor='k', alpha=0.15, lw=0,
				zorder=0,
			)
		)

	ax.plot(eq48_nom, eq47_nom, c='k', label='Equilibrium (Fiebig24)', zorder=1)

	for i, txt in enumerate(annotation):
		ax.scatter(
			ann_eq48[i], ann_eq47[i],
			c='w', marker='o', edgecolors='black', zorder=99,
		)
		if txt == 15:
			continue
		ax.annotate(f"{txt} °C", (ann_eq48[i] + OFFSETS[0], ann_eq47[i] + OFFSETS[1]))


## Sub-group definitions for the split panels.
GROUPS = [
	('17NTS', lambda s: s.str.contains('NTS')),
	('XSC',   lambda s: s.str.contains('XSC')),
]


def main(show_calib_ci=SHOW_CALIB_CI, calib_ci_level=CALIB_CI_LEVEL):
	summary, _reps = _load_data()

	for group_name, group_match in GROUPS:
		sub_summary = summary.loc[group_match(summary['Sample name'])].copy()

		fig, ax = plt.subplots(figsize=(6, 5))

		scatterSample(sub_summary, ax)
		plot_calib(ax, show_ci=show_calib_ci, ci_level=calib_ci_level)

		ax.set_xlim(0.21, 0.31)
		ax.set_ylim(0.6, 0.66)
		ax.set_aspect('equal', adjustable='box')
		#ax.set_title(group_name, fontsize=BIGGER_SIZE)
		ax.set_xlabel(r"$\mathit{\Delta_{48},\,_{CDES90} \,  (\permil)}$")
		ax.set_ylabel(r"$\mathit{\Delta_{47},\,_{CDES90} \,  (\permil)}$")

		ax.legend(
			loc='upper left',
			bbox_to_anchor=(1.02, 1.0),
			ncol=1,
			edgecolor='k',
			framealpha=1,
			borderaxespad=0,
			fancybox=False,
		)

		fig.tight_layout()

		out = f'D47D48_errorbar_{group_name}.pdf'
		fig.savefig(out)
		print(f"Saved {out}")

	plt.show()


if __name__ == "__main__":
	main()
