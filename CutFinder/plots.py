from cycler import cycler
import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import matplotlib as mpl
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
import numpy as np

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

ggplot_palette=('#348ABD','#E24A33', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8')

palette = ggplot_palette

hep.styles.cms.CMS["patch.linewidth"] = 2
hep.styles.cms.CMS["lines.linewidth"] = 2
hep.styles.cms.CMS["axes.prop_cycle"] = cycler("color", palette)
hep.styles.cms.CMS["legend.frameon"] = False
hep.styles.cms.CMS["figure.autolayout"] = True
hep.style.use(hep.style.CMS)

#mpl.use('agg')
mplstyle.use('fast')
mpl.rcParams['path.simplify'] = True
#mpl.rcParams['path.simplify_threshold'] = 0.5
#mpl.rcParams['agg.path.chunksize'] = 10000

markers = ["v", "^", "X", "P", "d", "*", "p", "o"]

class Plotter:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def plot_rates(self, glob, obj, output=None):
        isScaled = obj.isScaled #check if both are scaled already performed in main

        records = obj.records
        pt_bins = glob.pt_bins

        fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.})
        main_ax = ax[0]
        ratio_ax = ax[1]

        fitted_line_legend = Line2D([0], [0], color='black', lw=3, linestyle='-.')
        fitted_line_legend.set_path_effects([
            pe.Stroke(linewidth=1.5, foreground='magenta', alpha=0.7),
            pe.Normal()
        ])
        custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='--'),
                        Line2D([0], [0], color='black', lw=3, linestyle='-'),
                        fitted_line_legend]
        labels = ["Reference", obj.name + "(full)", obj.name + "(fitted)"]

        for idx, (ref, ref_dict) in enumerate(records.items()):
            #plot obj rate
            main_ax.plot(pt_bins, ref_dict["full"]["rate"],
                marker=markers[idx % len(markers)],
                linestyle="-",
                markeredgecolor='black',
                color=palette[idx % len(palette)],
                alpha=1,
                linewidth=3
            )


            #plot ref rate
            ref_line, = main_ax.plot(
                pt_bins,
                ref_dict["ref_rate"],
                marker=markers[idx % len(markers)],
                linestyle="--",
                color=palette[idx % len(palette)],
                alpha=0.666,
            )
            ref_line.set_path_effects([
                pe.Stroke(linewidth=3, foreground='black'),
                pe.Normal()
            ])

            fitted_line, = main_ax.plot(pt_bins, ref_dict["fitted"]["rate"],
                marker=markers[idx % len(markers)],
                linestyle="-.",
                markeredgecolor='black',
                color=palette[idx % len(palette)],
                linewidth=3
            )

            fitted_line.set_path_effects([
                pe.Stroke(linewidth=1.5, foreground='magenta', alpha=0.7),
                pe.Normal()
            ])

            ratio = np.array(ref_dict["full"]["rate"])/np.array(ref_dict["ref_rate"])
            ratio_fitted = np.array(ref_dict["fitted"]["rate"])/np.array(ref_dict["ref_rate"])

            ratio_ax.plot(pt_bins, ratio,
                marker=markers[idx % len(markers)],
                linestyle="-",
                markeredgecolor='black',
                color=palette[idx % len(palette)],
            )

            ratio_ax.plot(pt_bins, ratio_fitted,
                marker=markers[idx % len(markers)],
                linestyle="-.",
                markeredgecolor='black',
                color=palette[idx % len(palette)],
            )

            labels.append(f"{ref}")
            custom_lines.append(Line2D([0], [0], color=palette[idx % len(palette)], lw=3, linestyle='-', marker=markers[idx % len(markers)], markeredgecolor='black'))

        main_ax.legend(custom_lines, labels, loc='upper right', fontsize='small')

        #scale
        main_ax.set_yscale("log")
        #labels
        main_ax.set_title(f"{obj.name}")
        main_ax.set_ylabel("Rate [kHz]")
        ratio_ax.set_ylabel("Obj/Ref")
        if isScaled:
            ratio_ax.set_xlabel("Offline pT [GeV]")
        else:
            ratio_ax.set_xlabel("Online pT [GeV]")

        #grids
        ratio_ax.axhline(1.0, color="black", linestyle="--", linewidth=1, zorder = -99, alpha=0.7)
        main_ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)

        #lims
        main_ax.set_ylim(0.3, glob.maxRate*10)
        ratio_ax.set_ylim(0,2)
        fig.savefig(f"{output}/{obj.name}/rates.png")
        fig.savefig(f"{output}/{obj.name}/rates.pdf")
        plt.close(fig)

    def plot_cuts(self, obj, output=None):
        records = obj.records

        fig, ax = plt.subplots()
        custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='--'),
                        Line2D([0], [0], color='black', lw=2, linestyle='-')]
        labels = [obj.name, f"{obj.name} (fitted)"]
        for idx, (ref, ref_dict) in enumerate(records.items()):

            fitted_bins = np.array(ref_dict["fitted"]["bins"])
            fitted_cuts = np.array(ref_dict["fitted"]["cuts"])
            fitted_bins = fitted_bins[fitted_cuts>-9999].tolist()
            fitted_cuts = fitted_cuts[fitted_cuts>-9999].tolist()

            if fitted_bins[-1] != ref_dict["full"]["bins"][-1]:
                #extend fitted to the end of the full cuts for better visualization
                fitted_bins.append(ref_dict["full"]["bins"][-1])
                fitted_cuts.append(fitted_cuts[-1])

            ax.step(fitted_bins,
                    fitted_cuts,
                    linestyle="-",
                    color=palette[idx % len(palette)],
                    marker=markers[idx % len(markers)],
                    markeredgecolor='black',
                    alpha=1,
                    linewidth=3,
                    where="post",
            )

            full_bins = np.array(ref_dict["full"]["bins"])
            full_cuts = np.array(ref_dict["full"]["cuts"])
            full_bins = full_bins[full_cuts>-9999]
            full_cuts = full_cuts[full_cuts>-9999]

            #plot full cuts
            line, = ax.plot(full_bins,
                        full_cuts,
                        marker=markers[idx % len(markers)],
                        linestyle="--",
                        color=palette[idx % len(palette)],
                        alpha=0.666,)
            line.set_path_effects([
                pe.Stroke(linewidth=3, foreground='black'),
                pe.Normal()
            ])

            custom_lines.append(Line2D([0], [0], color=palette[idx % len(palette)], lw=3, linestyle='-', marker=markers[idx % len(markers)]))
            #labels.append(ref + f'\n$\chi^2$/ndf: {ref_dict["fitted"]["chi2"]:.4f}')
            labels.append(ref)

        ax.legend(custom_lines, labels, loc='upper right', fontsize='small')
        ax.set_title(f"{obj.name} Cuts")
        ax.set_xlabel("Online pT [GeV]")
        ax.set_ylabel("Cut Value (score)")

        y_min, y_max = ax.get_ylim()
        if y_min > -1 and y_max < 1:
            ax.set_ylim(-1, 1)
        ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.7)
        fig.savefig(f"{output}/{obj.name}/cuts.png")
        fig.savefig(f"{output}/{obj.name}/cuts.pdf")
        plt.close(fig)