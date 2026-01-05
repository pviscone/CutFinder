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

    def plot_rates(self, glob, obj, record_name, output=None):
        isScaled = obj.isScaled #check if both are scaled already performed in main

        records = obj.records
        pt_bins = glob.pt_bins

        fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.})
        main_ax = ax[0]
        ratio_ax = ax[1]
        custom_lines = [Line2D([0], [0], color='black', lw=2, linestyle='-')]
        labels = [obj.name]

        for idx, (ref, ref_dict) in enumerate(records.items()):
            #plot obj rate
            main_ax.plot(pt_bins, ref_dict[record_name]["rate"],
                marker=markers[idx % len(markers)],
                linestyle="-",
                color=palette[idx % len(palette)],
                alpha=1,
                linewidth=3
            )


            #plot ref rate
            line, = main_ax.plot(
                pt_bins,
                ref_dict["ref_rate"],
                marker=markers[idx % len(markers)],
                linestyle="--",
                color=palette[idx % len(palette)],
                alpha=0.666,
            )
            line.set_path_effects([
                pe.Stroke(linewidth=4, foreground='black'),
                pe.Normal()
            ])



            ratio = np.array(ref_dict[record_name]["rate"])/np.array(ref_dict["ref_rate"])
            ratio_ax.plot(pt_bins, ratio,
                marker=markers[idx % len(markers)],
                linestyle="--",
                color=palette[idx % len(palette)],
            )

            labels.append(f"{ref}")
            custom_lines.append(Line2D([0], [0], color=palette[idx % len(palette)], lw=3, linestyle='--', marker=markers[idx % len(markers)]))

        main_ax.legend(custom_lines, labels, loc='upper right', fontsize='small')

        #scale
        main_ax.set_yscale("log")
        #labels
        main_ax.set_title(f"{obj.name}: {record_name}")
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
        fig.savefig(f"{output}/{obj.name}/{record_name}.png")
        fig.savefig(f"{output}/{obj.name}/{record_name}.pdf")
        plt.close(fig)

    def plot_cuts(self):
        pass