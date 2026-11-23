"""Online-to-offline pT scaling derived from the efficiency turn-ons.

For every online pT cut the efficiency turn-on vs gen pT is fitted and the gen
pT at a target efficiency (the so-called "95% WP") is extracted; a linear fit of
those points gives the scaling offline = a * online + b.
"""

import os

import hist
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
from hist.intervals import ratio_uncertainty
from scipy.optimize import curve_fit, root_scalar
from scipy.stats import norm


def eff_func(x, par0, par1, par2, par3, par4):
    """Efficiency turn-on curve (error function + exponential tail)."""
    term1_x = par0 * (x - par1)
    term1 = norm.cdf(term1_x, loc=0, scale=par0 * par2)
    exponent = -par0 * (x - par1) + par0**2 * par2**2 / 2
    term2 = np.exp(exponent)
    term3 = norm.cdf(term1_x, loc=par0**2 * par2**2, scale=par0 * par2)
    return (term1 - term2 * term3) * (par3 - par4) + par4


def pt_response_kernels(gen_pt, reco_pt, matched, bins=None, min_entries=200):
    """Measured pT response of the online object, in bins of the gen pT.

    Returns the gen pT centres and, for each of them, the reco/gen pT ratios of
    the matched objects. The full measured distribution is kept as the kernel:
    no functional form (e.g. a Gaussian) is assumed, so the bremsstrahlung tail
    and any other non-Gaussian feature enters the smearing as measured. Bins
    with fewer than min_entries matched objects are merged forward until the
    kernel is statistically usable, and the centre is the median gen pT of the
    merged entries (not the bin midpoint).
    """
    gen_pt = np.asarray(gen_pt, dtype=float)
    reco_pt = np.asarray(reco_pt, dtype=float)
    matched = np.asarray(matched, dtype=bool)
    edges = np.arange(0, 102, 2) if bins is None else np.asarray(bins, dtype=float)

    centres, kernels = [], []
    lo = edges[0]
    for hi in edges[1:]:
        sel = matched & (gen_pt >= lo) & (gen_pt < hi)
        if sel.sum() < min_entries and hi < edges[-1]:
            continue  # widen this bin with the next one
        ratio = reco_pt[sel] / gen_pt[sel]
        if ratio.size:
            centres.append(float(np.median(gen_pt[sel])))
            kernels.append(ratio)
        lo = hi
    return np.array(centres), kernels


def sample_pt_ratio(gen_pt, centres, kernels, rng):
    """Sample a reco/gen pT ratio per gen electron from the measured kernels.

    The kernels are linearly interpolated in gen pT. Interpolating the CDFs is
    the two-component mixture (1 - w) K_left + w K_right, drawn here by picking
    the left kernel with probability 1 - w: no shape is assumed for the
    response, only that it varies smoothly with gen pT. Below the lowest
    measured centre the lowest kernel is used, above the highest the last one.
    """
    pt = np.asarray(gen_pt, dtype=float)
    centres = np.asarray(centres, dtype=float)
    if not kernels:
        return np.ones_like(pt)
    if len(kernels) == 1:
        return rng.choice(kernels[0], size=pt.size)

    x = np.clip(pt, centres[0], centres[-1])
    idx = np.clip(np.searchsorted(centres, x, side="right") - 1, 0, len(centres) - 2)
    w = (x - centres[idx]) / (centres[idx + 1] - centres[idx])
    use_right = rng.random(pt.size) < w

    out = np.empty(pt.size)
    for k in range(len(centres) - 1):
        for take_right, kernel in ((False, kernels[k]), (True, kernels[k + 1])):
            sel = (idx == k) & (use_right == take_right)
            if sel.any():
                out[sel] = rng.choice(kernel, size=int(sel.sum()))
    return out


def find_inverse(f, y_target, xmin, xmax, *params):
    sol = root_scalar(lambda x: f(x, *params) - y_target, bracket=[xmin, xmax])
    if not sol.converged:
        raise ValueError("Could not find the inverse for the given value.")
    return sol.root


def derive_scaling_params(
    reco_pt,
    gen_pt,
    matched,
    pt_cuts=None,
    target_eff=0.95,
    bins=(200, 3, 80),
    savefolder=None,
    label=None,
):
    """Derive the (a, b) of offline_pT = a * online_pT + b from the turn-ons.

    reco_pt/gen_pt/matched: per gen-electron arrays; unmatched entries (or
    entries outside the acceptance) must be masked out by the caller or left
    with matched=False. The turn-on denominator is the matched gen electrons,
    the numerator those with reco pT above the cut.

    If savefolder is given, the fitted turn-on of every cut and the linear fit
    of the 95% WP points are saved there (as in the TkEleRegression studies).
    """
    reco_pt = np.asarray(reco_pt)
    gen_pt = np.asarray(gen_pt)
    matched = np.asarray(matched, dtype=bool)
    if pt_cuts is None:
        pt_cuts = np.arange(15, 56, 10)
    if savefolder is not None:
        os.makedirs(savefolder, exist_ok=True)

    den_h = hist.Hist(hist.axis.Regular(*bins))
    den_h.fill(gen_pt[matched])
    y95 = []
    for cut in pt_cuts:
        pt_h = hist.Hist(hist.axis.Regular(*bins))
        pt_h.fill(gen_pt[matched & (reco_pt > cut)])
        eff_h = pt_h / den_h
        x = eff_h.axes[0].centers
        y = eff_h.values()
        y_err = ratio_uncertainty(
            pt_h.values(), den_h.values(), uncertainty_type="efficiency"
        )
        y_err = (y_err[0] + y_err[1]) / 2

        # empty bins have no efficiency and no error: drop them from the fit
        good = np.isfinite(y) & np.isfinite(y_err) & (y_err > 0)
        x, y, y_err = x[good], y[good], y_err[good]
        # data-driven start and physical bounds: a fixed start far from the
        # turn-on makes curve_fit wander into degenerate minima (overflowing
        # exp, plateau ~ 0.1) and the 95% WP bracket then has no sign change
        p0 = [0.4, cut, 0.8, y.max(), 0.0]
        bounds = ([1e-3, cut - 30, 1e-3, 0.5, -0.5], [10.0, cut + 30, 50.0, 1.05, 0.5])
        popt, _ = curve_fit(eff_func, x, y, p0=p0, sigma=y_err, bounds=bounds)
        y95.append(find_inverse(eff_func, target_eff, cut - 10, cut + 20, *popt))

        if savefolder is not None:
            fig, ax = plt.subplots()
            ax.errorbar(
                x,
                y,
                yerr=y_err,
                fmt="o",
                color="C0",
                markersize=5,
                label=label,
                markeredgecolor="black",
            )
            xx = np.linspace(x[0], x[-1], 100)
            ax.plot(
                xx,
                eff_func(xx, *popt),
                color="C0",
                alpha=0.35,
                linewidth=3,
                linestyle="--",
            )
            ax.plot(
                y95[-1],
                target_eff,
                "*",
                color="white",
                markeredgecolor="C0",
                markersize=20,
                zorder=999,
            )
            ax.axvline(
                cut,
                color="black",
                linestyle="--",
                linewidth=1,
                alpha=0.5,
                label=f"Online $p_T$ cut: {cut:.2f} GeV",
            )
            ax.axhline(target_eff, color="black", linestyle=":", linewidth=1, alpha=0.3)
            ax.set_xlim(cut - 10, cut + 20)
            ax.set_ylim(-0.05, 1.09)
            ax.set_xlabel("Gen $p_T$ [GeV]")
            ax.set_ylabel("Efficiency")
            ax.legend(fontsize=14)
            hep.cms.label(
                llabel="Phase-2 Simulation Preliminary",
                rlabel="PU 200 (14 TeV)",
                ax=ax,
                fontsize=13,
            )
            fig.savefig(f"{savefolder}/turn_on_{cut:.2f}.png", bbox_inches="tight")
            fig.savefig(f"{savefolder}/turn_on_{cut:.2f}.pdf", bbox_inches="tight")
            plt.close(fig)
    y95 = np.array(y95)

    def linear(x, a, b):
        return a * x + b

    (a, b), _ = curve_fit(linear, pt_cuts, y95)

    if savefolder is not None:
        fig, ax = plt.subplots()
        ax.plot(pt_cuts, y95, ".", color="C0", markersize=14, label=label)
        ax.plot(
            pt_cuts,
            linear(pt_cuts, a, b),
            color="C0",
            alpha=0.35,
            linewidth=3,
            label=f"{label}: {a:.2f} $p_T^{{L1}}$ + {b:.2f} GeV",
        )
        ax.legend(fontsize=16)
        ax.set_xlabel("Online $p_T$ cut [GeV]")
        ax.set_ylabel(f"$p_T^{{{target_eff * 100:.0f}\\%}}$ [GeV]")
        hep.cms.label(
            llabel="Phase-2 Simulation Preliminary",
            rlabel="PU 200 (14 TeV)",
            ax=ax,
            fontsize=13,
        )
        fig.savefig(f"{savefolder}/online_to_offline_scaling.png", bbox_inches="tight")
        fig.savefig(f"{savefolder}/online_to_offline_scaling.pdf", bbox_inches="tight")
        plt.close(fig)

    return float(a), float(b)
