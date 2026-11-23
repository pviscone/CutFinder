from cycler import cycler
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib as mpl
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import numpy as np

from .scaling import pt_response_kernels, sample_pt_ratio

cms10 = [
    "#3f90da",
    "#ffa90e",
    "#bd1f01",
    "#94a4a2",
    "#832db6",
    "#a96b59",
    "#e76300",
    "#b9ac70",
    "#92dadd",
    "#717581",
]

acab_palette = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)

ggplot_palette = (
    "#348ABD",
    "#E24A33",
    "#988ED5",
    "#777777",
    "#FBC15E",
    "#8EBA42",
    "#FFB5B8",
)

palette = ggplot_palette

hep.styles.cms.CMS["patch.linewidth"] = 2
hep.styles.cms.CMS["lines.linewidth"] = 2
hep.styles.cms.CMS["axes.prop_cycle"] = cycler("color", palette)
hep.styles.cms.CMS["legend.frameon"] = False
hep.styles.cms.CMS["figure.autolayout"] = True
hep.style.use(hep.style.CMS)

# mpl.use('agg')
mplstyle.use("fast")
mpl.rcParams["path.simplify"] = True
# mpl.rcParams['path.simplify_threshold'] = 0.5
# mpl.rcParams['agg.path.chunksize'] = 10000

markers = ["v", "^", "X", "P", "d", "*", "p", "o"]


def display_names(*configs):
    """Map config names to their display names for plot labels (fallback: the name)."""
    return {cfg.name: getattr(cfg, "display_name", None) or cfg.name for cfg in configs}


class Plotter:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def plot_rates(self, glob, obj, refs, output=None, online_axis=False):
        isScaled = obj.isScaled  # check if both are scaled already performed in main

        names = display_names(obj, *refs)
        obj_label = names[obj.name]

        records = obj.records
        # the rate comparison is on the offline pT scale: the x axis is the
        # offline threshold (obj scaling); online_axis=True shows the same
        # curves vs the online threshold (the configured pt_bins)
        if online_axis:
            pt_bins = glob.pt_bins
            xlabel = "Online pT [GeV]"
            name = "rates_vs_online"
        else:
            pt_bins = (
                obj.forward_scaling(np.array(glob.pt_bins)) if isScaled else glob.pt_bins
            )
            xlabel = "Offline pT [GeV]" if isScaled else "Online pT [GeV]"
            name = "rates"

        fig, ax = plt.subplots(
            2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0}
        )
        main_ax = ax[0]
        ratio_ax = ax[1]

        handles = []
        labels = []

        for idx, (ref, ref_dict) in enumerate(records.items()):
            # same colors as the efficiency plots: ref, obj (full), obj (fitted)
            c_ref = palette[(3 * idx) % len(palette)]
            c_full = palette[(3 * idx + 1) % len(palette)]
            c_fitted = palette[(3 * idx + 2) % len(palette)]
            # plot obj rate
            (full_line,) = main_ax.plot(
                pt_bins,
                ref_dict["full"]["rate"],
                marker=markers[idx % len(markers)],
                linestyle="-",
                markeredgecolor="black",
                color=c_full,
                alpha=1,
                linewidth=3,
            )

            # plot ref rate
            (ref_line,) = main_ax.plot(
                pt_bins,
                ref_dict["ref_rate"],
                marker=markers[idx % len(markers)],
                linestyle="--",
                color=c_ref,
                alpha=0.666,
            )
            ref_line.set_path_effects(
                [pe.Stroke(linewidth=3, foreground="black"), pe.Normal()]
            )

            (fitted_line,) = main_ax.plot(
                pt_bins,
                ref_dict["fitted"]["rate"],
                marker=markers[idx % len(markers)],
                linestyle="-.",
                markeredgecolor="black",
                color=c_fitted,
                linewidth=3,
            )

            fitted_line.set_path_effects(
                [pe.Stroke(linewidth=1.5, foreground="magenta", alpha=0.7), pe.Normal()]
            )

            ratio = np.array(ref_dict["full"]["rate"]) / np.array(ref_dict["ref_rate"])
            ratio_fitted = np.array(ref_dict["fitted"]["rate"]) / np.array(
                ref_dict["ref_rate"]
            )

            ratio_ax.plot(
                pt_bins,
                ratio,
                marker=markers[idx % len(markers)],
                linestyle="-",
                markeredgecolor="black",
                color=c_full,
            )

            ratio_ax.plot(
                pt_bins,
                ratio_fitted,
                marker=markers[idx % len(markers)],
                linestyle="-.",
                markeredgecolor="black",
                color=c_fitted,
            )

            handles += [full_line, fitted_line, ref_line]
            labels += [
                f"{obj_label} (full)",
                f"{obj_label} (fitted)",
                names.get(ref, ref),
            ]

        main_ax.legend(handles, labels, loc="upper right", fontsize="small")

        # scale
        main_ax.set_yscale("log")
        # labels
        hep.cms.label(
            llabel="Phase-2 Simulation Preliminary",
            rlabel="PU 200 (14 TeV)",
            ax=main_ax,
            fontsize=13,
        )
        main_ax.set_ylabel("Rate [kHz]")
        ratio_ax.set_ylabel("Obj/Ref")
        ratio_ax.set_xlabel(xlabel)

        # grids
        ratio_ax.axhline(
            1.0, color="black", linestyle="--", linewidth=1, zorder=-99, alpha=0.7
        )
        main_ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)

        # lims
        main_ax.set_ylim(0.3, glob.maxRate * 10)
        ratio_ax.set_ylim(0, 2)
        if isScaled and not online_axis:
            # the offline threshold axis starts at the scaling offset: show it from 0
            main_ax.set_xlim(left=0)
        fig.savefig(f"{output}/{obj.name}/{name}.png", bbox_inches="tight")
        fig.savefig(f"{output}/{obj.name}/{name}.pdf", bbox_inches="tight")
        plt.close(fig)

    def plot_cuts(self, obj, refs, output=None, offline_axis=False):
        records = obj.records

        names = display_names(obj, *refs)
        obj_label = names[obj.name]

        fig, ax = plt.subplots()
        handles = []
        labels = []
        for idx, (ref, ref_dict) in enumerate(records.items()):
            fitted_bins = np.array(ref_dict["fitted"]["bins"])
            fitted_cuts = np.array(ref_dict["fitted"]["cuts"])
            fitted_bins = fitted_bins[fitted_cuts > -9999].tolist()
            fitted_cuts = fitted_cuts[fitted_cuts > -9999].tolist()

            if fitted_bins[-1] != ref_dict["full"]["bins"][-1]:
                # extend fitted to the end of the full cuts for better visualization
                fitted_bins.append(ref_dict["full"]["bins"][-1])
                fitted_cuts.append(fitted_cuts[-1])

            if offline_axis:
                # record bins are in the online scale (obj scaling)
                fitted_bins = obj.forward_scaling(np.asarray(fitted_bins, dtype=float))

            (fitted_line,) = ax.step(
                fitted_bins,
                fitted_cuts,
                linestyle="-",
                color=palette[idx % len(palette)],
                marker=markers[idx % len(markers)],
                markeredgecolor="black",
                alpha=1,
                linewidth=3,
                where="post",
            )

            full_bins = np.array(ref_dict["full"]["bins"])
            full_cuts = np.array(ref_dict["full"]["cuts"])
            full_err = np.array(ref_dict["full"]["cuts_err"])
            full_bins = full_bins[full_cuts > -9999]
            full_cuts = full_cuts[full_cuts > -9999]
            full_err = full_err[full_cuts > -9999]
            if offline_axis:
                full_bins = obj.forward_scaling(full_bins.astype(float))

            # plot full cuts
            line, _, _ = ax.errorbar(
                full_bins,
                full_cuts,
                yerr=full_err,
                marker=markers[idx % len(markers)],
                linestyle="--",
                color=palette[idx % len(palette)],
                alpha=0.666,
            )
            line.set_path_effects(
                [pe.Stroke(linewidth=3, foreground="black"), pe.Normal()]
            )

            handles += [line, fitted_line]
            suffix = "" if len(records) == 1 else f" ({names.get(ref, ref)})"
            labels += [f"{obj_label} (full){suffix}", f"{obj_label} (fitted){suffix}"]

        ax.legend(handles, labels, loc="upper right", fontsize="small")
        hep.cms.label(
            llabel="Phase-2 Simulation Preliminary",
            rlabel="PU 200 (14 TeV)",
            ax=ax,
            fontsize=13,
        )
        ax.set_xlabel("Offline pT [GeV]" if offline_axis else "Online pT [GeV]")
        ax.set_ylabel("Cut Value (score)")

        y_min, y_max = ax.get_ylim()
        if y_min > -1 and y_max < 1:
            ax.set_ylim(-1, 1)
        if offline_axis:
            # the offline threshold axis starts at the scaling offset: show it from 0
            ax.set_xlim(left=0)
        ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
        name = "cuts_vs_offline" if offline_axis else "cuts"
        fig.savefig(f"{output}/{obj.name}/{name}.png", bbox_inches="tight")
        fig.savefig(f"{output}/{obj.name}/{name}.pdf", bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def _eff_curve(x, matched, bins):
        """Efficiency histogram: count(matched)/count(all) with completeness errors."""
        from hist.intervals import ratio_uncertainty

        num, _ = np.histogram(x[matched], bins)
        den, _ = np.histogram(x, bins)
        eff = np.divide(num, den, out=np.full_like(den, np.nan, dtype=float), where=den > 0)
        err_lo, err_hi = ratio_uncertainty(num, den, uncertainty_type="efficiency")
        with np.errstate(invalid="ignore"):
            err = np.nanmean(np.stack([err_lo, err_hi]), axis=0)
        centers = 0.5 * (bins[:-1] + bins[1:])
        return centers, eff, err

    def plot_efficiencies(self, obj, refs, output=None):
        """Gen-level efficiency of the WPs: reference vs obj(full) vs obj(fitted).

        One figure per (ref, variable): variable in {gen pt, gen eta, online pt,
        offline pt}. Every variable bins the same set of accepted gen electrons,
        so the denominator is the full gen set in all of them. The online and
        offline pT axes replace the gen pT of every electron by pt * ratio, the
        ratio sampled from the measured reco/gen distribution of the matched
        probes (no functional form assumed), interpolated in gen pT. The
        efficiency turn-on is therefore diluted by the pT resolution and moved
        to the scale on which the object is actually selected, while gen
        electrons that never reconstruct stay in the denominator.
        The reference curve is the ref's own ID selection; obj curves apply the
        per-bin ("full") and regression-fitted ("fitted") score WPs found against
        that ref.
        """
        ref_by_name = {ref.name: ref for ref in refs}
        names = display_names(obj, *refs)
        obj_label = names[obj.name]
        rng = np.random.default_rng(42)  # reproducible kernel sampling
        for ref_name, rec in obj.records.items():
            ref = ref_by_name.get(ref_name)
            curves = {}
            if ref is not None and ref.eff_config is not None:
                curves[names.get(ref_name, ref_name)] = (
                    ref,
                    ref.eff_config,
                    None,
                )  # ref: no score WP
            curves[f"{obj_label} (full)"] = (
                obj,
                obj.eff_config,
                (rec["full"]["bins"], rec["full"]["cuts"]),
            )
            curves[f"{obj_label} (fitted)"] = (
                obj,
                obj.eff_config,
                (rec["fitted"]["bins"], rec["fitted"]["cuts"]),
            )
            # one event loop per curve, shared by all the x variables
            datasets = {
                label: eff_config.get_gen_data(wp=wp)
                for label, (_, eff_config, wp) in curves.items()
            }

            # measured pT response of every curve: the reco/gen ratio kernel per
            # gen pT bin, from the matched probes (before any ID selection)
            responses = {}
            for label, (_, eff_config, _) in curves.items():
                data = datasets[label]
                responses[label] = pt_response_kernels(
                    data["pt"],
                    data["probe_pt"],
                    data["probe_matched"],
                    np.asarray(eff_config.pt_bins, dtype=float),
                )

            # the response itself, as measured: median and 68% interval of the
            # reco/gen ratio per gen pT bin
            fig, ax = plt.subplots()
            for idx, (label, (centres, kernels)) in enumerate(responses.items()):
                if not kernels:
                    continue
                med = np.array([np.median(k) for k in kernels])
                lo = np.array([np.percentile(k, 16) for k in kernels])
                hi = np.array([np.percentile(k, 84) for k in kernels])
                ax.errorbar(
                    centres,
                    med,
                    yerr=[med - lo, hi - med],
                    marker=markers[idx % len(markers)],
                    linestyle="-",
                    linewidth=3,
                    color=palette[idx % len(palette)],
                    markeredgecolor="black",
                    label=label,
                )
            ax.axhline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.5)
            ax.set_xlabel("Gen $p_T$ [GeV]")
            ax.set_ylabel("$p_T^{\\mathrm{online}} / p_T^{\\mathrm{gen}}$")
            ax.set_xlim(0, 60)
            ax.legend(fontsize="small")
            hep.cms.label(
                llabel="Phase-2 Simulation Preliminary",
                rlabel="PU 200 (14 TeV)",
                ax=ax,
                fontsize=13,
            )
            response_name = "onlinePt_response"
            if len(obj.records) > 1:
                response_name = f"{response_name}_vs_{ref_name}"
            fig.savefig(f"{output}/{obj.name}/{response_name}.png", bbox_inches="tight")
            fig.savefig(f"{output}/{obj.name}/{response_name}.pdf", bbox_inches="tight")
            plt.close(fig)

            # the offline pT axis needs the online-to-offline scaling of every
            # config: skip it in --no-scaling runs
            variables = ["pt", "eta", "online_pt"]
            if all(cfg.forward_scaling is not None for cfg, _, _ in curves.values()):
                variables.append("offline_pt")
            else:
                print(f"Skipping the offline pT efficiency for {ref_name}: no scaling")

            for var in variables:
                fig, (main_ax, ratio_ax) = plt.subplots(
                    2,
                    1,
                    sharex=True,
                    gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0},
                )
                base_eff = base_centers = None
                for idx, (label, (cfg, eff_config, _)) in enumerate(curves.items()):
                    data = datasets[label]
                    if var == "pt":
                        x = data["pt"]
                        bins = eff_config.pt_bins
                    elif var == "eta":
                        x = data["eta"]
                        bins = eff_config.eta_bins
                    else:
                        # gen pT moved to the online (or offline) scale with the
                        # measured resolution kernel of the curve
                        x = data["pt"] * sample_pt_ratio(
                            data["pt"], *responses[label], rng
                        )
                        bins = np.asarray(eff_config.pt_bins, dtype=float)
                        if var == "offline_pt":
                            x = cfg.forward_scaling(x)
                            bins = cfg.forward_scaling(bins)
                    matched = data["matched"]
                    centers, eff, err = self._eff_curve(x, matched, bins)
                    linestyle = {"full": "-", "fitted": "-."}.get(
                        label.split("(")[-1].rstrip(")"), "--"
                    )
                    main_ax.errorbar(
                        centers,
                        eff,
                        yerr=err,
                        marker=markers[idx % len(markers)],
                        linestyle=linestyle,
                        linewidth=3,
                        color=palette[idx % len(palette)],
                        markeredgecolor="black",
                        label=label,
                    )
                    if base_eff is None:
                        base_eff = eff  # reference curve for the ratio panel
                        base_centers = centers
                    elif base_eff is not None:
                        with np.errstate(divide="ignore", invalid="ignore"):
                            ratio = eff / base_eff
                        ratio_ax.plot(
                            base_centers,
                            ratio,
                            marker=markers[idx % len(markers)],
                            linestyle=linestyle,
                            color=palette[idx % len(palette)],
                            markeredgecolor="black",
                        )

                main_ax.legend(loc="lower right", fontsize="small")
                hep.cms.label(
                    llabel="Phase-2 Simulation Preliminary",
                    rlabel="PU 200 (14 TeV)",
                    ax=main_ax,
                    fontsize=13,
                )
                main_ax.set_ylabel("Efficiency")
                main_ax.set_ylim(0, 1.09)
                main_ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
                ratio_ax.set_ylabel("Obj/Ref")
                ratio_ax.set_ylim(0.6, 1.4)
                ratio_ax.axhline(
                    1.0, color="black", linestyle="--", linewidth=1, zorder=-99, alpha=0.7
                )
                if var == "pt":
                    ratio_ax.set_xlabel("Gen $p_T$ [GeV]")
                    ratio_ax.set_xlim(0, 100)
                    suffix = "efficiency_vs_genPt"
                elif var == "eta":
                    ratio_ax.set_xlabel("Gen $\\eta$")
                    suffix = "efficiency_vs_genEta"
                elif var == "online_pt":
                    ratio_ax.set_xlabel("Online $p_T$ [GeV]")
                    ratio_ax.set_xlim(0, 100)
                    suffix = "efficiency_vs_onlinePt"
                else:
                    ratio_ax.set_xlabel("Offline $p_T$ [GeV]")
                    ratio_ax.set_xlim(left=0)
                    suffix = "efficiency_vs_offlinePt"
                if len(obj.records) > 1:
                    suffix = f"{suffix}_vs_{ref_name}"
                fig.savefig(f"{output}/{obj.name}/{suffix}.png", bbox_inches="tight")
                fig.savefig(f"{output}/{obj.name}/{suffix}.pdf", bbox_inches="tight")
                plt.close(fig)
