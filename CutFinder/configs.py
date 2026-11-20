from CutFinder.algorithms import iterative_bin_cutter
from CutFinder.functions import applyWP, wp_mask_strings
from CutFinder.extractor import extract_columns
from CutFinder.regressors import bayesian_blocks_gaussian
from CutFinder.scaling import derive_scaling_params

from typing import Optional
from collections.abc import Callable, Iterable

from array import array
import hist
import numpy as np
from rich import print as pprint
import ROOT
import sympy as sp


class GlobalConf:
    def __init__(
        self,
        pt_bins,
        maxRate=31038.96,
        algo=iterative_bin_cutter,
        algo_kwargs={},
        regressor=bayesian_blocks_gaussian,
        regressor_kwargs={},
        fitrange: Optional[tuple[float, float]] = None,
        scaling_kwargs={},
    ):  # (PU200)
        self.pt_bins = pt_bins
        self.maxRate = maxRate
        self.algo = algo
        self.regressor = regressor
        self.fitrange = fitrange
        self.algo_kwargs = algo_kwargs
        self.regressor_kwargs = regressor_kwargs
        self.scaling_kwargs = scaling_kwargs


# Keepalive registry: RDataFrame does NOT take ownership of the TChain; if the
# chain is garbage-collected, every lazy action on the rdf dereferences freed
# memory (segfault, nondeterministic). Keep every chain alive for the process.
_chains: list = []


def open_rdf(samples_path: str, tree: str = "Events"):
    """Open an RDataFrame on a TChain; supports {start..end} numeric ranges in the path."""
    chain = ROOT.TChain(tree)
    if "{" in samples_path or "}" in samples_path:
        range_part = samples_path.split("{")[1].split("}")[0]
        start, end = map(int, range_part.split(".."))
        for i in range(start, end + 1):
            chain.Add(samples_path.replace("{" + range_part + "}", str(i)))
    else:
        chain.Add(samples_path)
    _chains.append(chain)
    return ROOT.RDataFrame(chain)


class Config:
    def __init__(
        self,
        *,
        samples_path: Optional[str] = None,
        pt_branch: str,
        score_branch: str,
        preprocess_function: Optional[Callable] = None,
        scaling_function: Optional[Callable] = None,
        rate: Optional[np.array] = None,
        tree: str = "Events",
        eff_config: Optional["ConfigEff"] = None,
    ):
        if samples_path is None:
            assert rate is not None

        self.samples_path = samples_path
        self.pt_branch = pt_branch
        self.score_branch = score_branch
        self.scaling_function = scaling_function

        self.func = preprocess_function
        self.rate = rate
        self.tree = tree
        self.eff_config = eff_config

        self.scaling = None
        self.forward_scaling = None
        self.inverse_scaling = None
        self.scaling_params = None  # (a, b) when derived from the eff_config turn-on
        if scaling_function is not None:
            self.set_scaling(scaling_function)

        ####### to compute
        self.name = None
        self.display_name = None  # optional pretty name for plot labels
        self.rdf = None
        self.TotEvents = None
        self.nEvents = None
        self.rate_err = None
        self.WP = None

        # flags
        self.isPreprocessed = False
        self.isScaled = False
        self.isWPApplied = False
        self.isComputed = False

    def set_scaling(self, scaling_function):
        """Install an online-to-offline scaling, offline = scaling_function(online).

        scaling_function maps a sympy symbol to a sympy expression. The pt_bins
        in GlobalConf are the ONLINE pT binning of the obj (the WP thresholds);
        the rate scan runs on the corresponding offline bins, so that obj and
        ref are compared at equal offline thresholds.
        """
        print(
            "Warning: with a scaling function pt_bins in GlobalConf is the online "
            "pT binning; the rate comparison runs at equal OFFLINE pT thresholds."
        )
        print(
            "If both scaling_function and rate are provided, the rate is not re-computed."
        )
        x = sp.Symbol("$OnlinePt")
        y = sp.Symbol("$OfflinePt")
        expr = scaling_function(x)
        inverse_scaling = sp.solve(sp.Eq(y, expr), x)[-1]
        self.scaling_function = scaling_function
        self.scaling = sp.ccode(expr).replace("$OnlinePt", self.pt_branch)
        self.forward_scaling = sp.lambdify(x, expr, "numpy")
        self.inverse_scaling = sp.lambdify(y, inverse_scaling, "numpy")

    def derive_scaling(self, savefolder=None, label=None, **kwargs):
        """Derive the online-to-offline scaling from the eff_config turn-on.

        Fits the efficiency turn-on vs gen pT for a set of online pT cuts,
        extracts the 95% WP of each and linear-fits them (turn-on fit + linear
        map, as in the TkEleRegression rate studies). The fitted turn-ons and
        the linear fit are saved in savefolder. Extra kwargs (pt_cuts,
        target_eff, bins) go to CutFinder.scaling.derive_scaling_params.
        """
        if self.eff_config is None:
            raise ValueError(f"{self.name}: derive_scaling requires an eff_config.")
        data = self.eff_config.get_gen_data()
        self.scaling_params = derive_scaling_params(
            data["reco_pt"],
            data["pt"],
            data["matched"],
            savefolder=savefolder,
            label=label or self.display_name or self.name,
            **kwargs,
        )
        a, b = self.scaling_params
        self.set_scaling(lambda x, a=a, b=b: a * x + b)
        return self.scaling_params

    def clone(self, cls=None, **kwargs):
        if cls is None:
            cls = self.__class__
        init = cls.__init__
        params = init.__code__.co_varnames
        base_kwargs = {
            "samples_path": self.samples_path,
            "pt_branch": self.pt_branch,
            "score_branch": self.score_branch,
            "preprocess_function": self.func,
            "scaling_function": self.scaling_function,
            "rate": self.rate,
            "tree": self.tree,
            "eff_config": self.eff_config,
        }
        base_kwargs.update(kwargs)
        filtered = {k: v for k, v in base_kwargs.items() if k in params}
        return cls(**filtered)

    def clone_reusing_rdf(self, cls=None, **kwargs):
        """Like clone(), but reuses the already-processed rdf (with its
        in-memory column cache) instead of reopening the samples and
        re-running the preprocessing from scratch.

        Only valid if the clone keeps the same samples_path and
        preprocess_function, so the processed columns are identical.
        """
        new = self.clone(cls, **kwargs)
        if not self.isPreprocessed or self.rdf is None:
            raise RuntimeError(
                "clone_reusing_rdf() requires the source Config to be already computed."
            )
        if new.samples_path != self.samples_path or new.func != self.func:
            raise ValueError(
                "clone_reusing_rdf() requires the same samples_path and "
                "preprocess_function; use clone() to change them."
            )
        new.rdf = self.rdf
        new.TotEvents = self.TotEvents
        new.isPreprocessed = True
        # the shared rdf is already scaled when the source was computed with a
        # scaling function: the clone must not scale it a second time
        new.isScaled = self.isScaled
        return new

    def loadRDF(self):
        if self.rdf is None:
            self.rdf = open_rdf(self.samples_path, self.tree)
            self.TotEvents = self.rdf.Count().GetValue()

    def runPreprocess(self):
        if self.func is not None and not self.isPreprocessed:
            self.rdf = self.func(self.rdf)
            self.isPreprocessed = True

    def getEntries(self):
        if self.nEvents is None:
            self.nEvents = self.rdf.Count().GetValue()

    def scale(self):
        if self.scaling is not None and not self.isScaled:
            # snapshot the current column in an alias: a chained
            # self-referential Redefine is a silent no-op when the preprocess
            # already redefined the branch (e.g. tight_barrel), leaving the
            # branch unscaled and the rate thresholds on the wrong scale
            snapshot = f"{self.pt_branch}_online"
            self.rdf = self.rdf.Alias(snapshot, self.pt_branch)
            self.rdf = self.rdf.Redefine(
                self.pt_branch, self.scaling.replace(self.pt_branch, snapshot)
            )
            self.isScaled = True

    def makeRate(self, bins, maxRate, overwrite=False):
        if self.rate is None or overwrite:
            if "MAXPT" not in self.rdf.GetDefinedColumnNames():
                self.rdf = self.rdf.Define("MAXPT", f"Max({self.pt_branch})")
            else:
                self.rdf = self.rdf.Redefine("MAXPT", f"Max({self.pt_branch})")

            if isinstance(bins, tuple):
                axis = hist.axis.Regular(bins[0], bins[1], bins[2])
                th = self.rdf.Histo1D(("", "", bins[0], bins[1], bins[2]), "MAXPT")
            elif isinstance(bins, np.ndarray):
                axis = hist.axis.Variable(bins)
                th = self.rdf.Histo1D(
                    ("", "", len(bins) - 1, array("d", bins)), "MAXPT"
                )
            else:
                raise ValueError("bins must be either a tuple or a numpy array.")

            th = th.GetValue()
            h = hist.Hist(axis, storage=hist.storage.Weight())

            # re-fill the real bins at their CENTER (under/overflow stay out of
            # the bin contents): a fill at the low edge lands in the wrong bin
            # as soon as the edge is not exactly representable, corrupting the
            # cumulative rate and with it the cut targets
            n_bins = th.GetNbinsX()
            for i in range(1, n_bins + 1):
                center = 0.5 * (th.GetBinLowEdge(i) + th.GetBinLowEdge(i + 1))
                h.fill(center, weight=th.GetBinContent(i))

            # Rate vs threshold (one entry per bin edge), same definition as the
            # cut search (count * maxRate / TotEvents). Events above the last
            # edge still pass every threshold: the overflow enters every point
            # of the cumulative (its own point is the rate at the last edge),
            # otherwise the curve dives to zero while the cut search counts
            # those events.
            overflow = th.GetBinContent(n_bins + 1)
            if overflow:
                h.fill(np.inf, weight=overflow)  # visible to the top-edge lookup
            counts = np.append(
                [th.GetBinContent(i) for i in range(1, n_bins + 1)], overflow
            )
            variances = np.append(
                [th.GetBinError(i) ** 2 for i in range(1, n_bins + 1)],
                th.GetBinError(n_bins + 1) ** 2,
            )
            scale = maxRate / self.TotEvents
            self.rate = np.cumsum(counts[::-1])[::-1] * scale
            self.rate_err = np.sqrt(np.cumsum(variances[::-1])[::-1]) * scale
            return h

    def compute(self):
        if not self.isComputed:
            if self.name is not None:
                pprint(
                    f"[bold green]Computing {self.__class__.__name__}:[/bold green]\n\t{self.name}\n"
                )
            self.loadRDF()
            # Materialize the pt (and score) columns in memory once: the
            # per-pt-bin cut search and WP-rate evaluation run several event
            # loops on this rdf, and without a cache each loop would re-run
            # the whole chain from the source files.
            if not self.isPreprocessed:
                self.runPreprocess()
                cache_cols = [self.pt_branch]
                if self.score_branch is not None:
                    cache_cols.append(self.score_branch)
                self.rdf = self.rdf.Cache(cache_cols)
            self.apply_WP()
            self.rdf = self.rdf.Filter(f"{self.pt_branch}.size()>0")
            self.scale()
            self.getEntries()
            self.isComputed = True

            # self.makeRate(pt_bins, maxRate)

    def apply_WP(self):
        if self.WP is not None and not self.isWPApplied:
            pt_bins = self.WP[0]
            if self.isScaled:
                # WP bins are in the online scale; the pt branch already holds
                # the offline pt (clones reusing a scaled rdf)
                pt_bins = self.forward_scaling(np.asarray(pt_bins, dtype=float))
            self.rdf = applyWP(
                self.pt_branch, self.score_branch, pt_bins, self.WP[1], self.rdf
            )
            self.isWPApplied = True


class ConfigEff:
    """Gen-level efficiency evaluation attached to a ConfigObj/ConfigRef.

    Signal sample -> preprocess_function -> (optional) WP applied to the reco
    objects -> gen electrons geometrically matched (dR) to a surviving reco
    object. Efficiency = gen[matched] / gen[all], vs gen pt and vs gen eta.
    """

    def __init__(
        self,
        *,
        samples_path: str,
        reco: dict,  # {"pt": ..., "eta": ..., "phi": ..., "score": ... (needed for WP)}
        gen: dict,  # {"pt": ..., "eta": ..., "phi": ...}
        preprocess_function: Optional[Callable] = None,
        tree: str = "Events",
        deltaR: float = 0.1,
        gen_acceptance: Optional[str] = None,  # e.g. "abs(GenEl_caloeta) < 1.479"
        pt_bins: Optional[np.ndarray] = None,
        eta_bins: Optional[np.ndarray] = None,
    ):
        self.samples_path = samples_path
        self.reco = reco
        self.gen = gen
        self.func = preprocess_function
        self.tree = tree
        self.deltaR = deltaR
        self.gen_acceptance = gen_acceptance
        self.pt_bins = np.arange(0, 102, 2) if pt_bins is None else np.asarray(pt_bins)
        self.eta_bins = (
            np.arange(-1.6, 1.61, 0.1) if eta_bins is None else np.asarray(eta_bins)
        )

    def get_gen_data(self, wp=None):
        """Run one event loop; return flat arrays {"pt", "eta", "matched", "reco_pt"}.

        wp: optional (pt_bins, score_cuts) applied to the reco objects before
        the gen matching. Requires reco["score"].
        """
        rdf = open_rdf(self.samples_path, self.tree)
        if self.func is not None:
            rdf = self.func(rdf)
        reco_eta, reco_phi, reco_pt = self.reco["eta"], self.reco["phi"], self.reco["pt"]
        if wp is not None:
            # Mask the reco objects per event; NO event Filter: events where no
            # reco object passes the WP still contribute their gen electrons to
            # the denominator as unmatched.
            if "score" not in self.reco:
                raise ValueError("Applying a WP requires reco['score'] in ConfigEff")
            pt_bins, score_cuts = wp_mask_strings(wp[0], wp[1])
            rdf = rdf.Define(
                "gen_eff_WPmask",
                f"WP_mask({self.reco['pt']}, {self.reco['score']}, {pt_bins}, {score_cuts})",
            )
            rdf = rdf.Define("gen_eff_recoEta", f"{reco_eta}[gen_eff_WPmask]")
            rdf = rdf.Define("gen_eff_recoPhi", f"{reco_phi}[gen_eff_WPmask]")
            rdf = rdf.Define("gen_eff_recoPt", f"{reco_pt}[gen_eff_WPmask]")
            reco_eta, reco_phi = "gen_eff_recoEta", "gen_eff_recoPhi"
            reco_pt = "gen_eff_recoPt"
        rdf = rdf.Define(
            "gen_matched",
            f"gen_reco_match({self.gen['eta']}, {self.gen['phi']}, "
            f"{reco_eta}, {reco_phi}, {self.deltaR})",
        )
        rdf = rdf.Define(
            "gen_matched_pt",
            f"gen_reco_match_pt({self.gen['eta']}, {self.gen['phi']}, "
            f"{reco_pt}, {reco_eta}, {reco_phi}, {self.deltaR})",
        )
        gen_pt, gen_eta, matched = self.gen["pt"], self.gen["eta"], "gen_matched"
        matched_pt = "gen_matched_pt"
        if self.gen_acceptance is not None:
            # Restrict the gen denominator to the geometric acceptance (e.g.
            # barrel), so that gen electrons no reco collection can ever match
            # do not pollute the efficiency.
            rdf = rdf.Define("gen_eff_acc", self.gen_acceptance)
            rdf = (rdf.Define("gen_eff_pt", f"{gen_pt}[gen_eff_acc]")
                   .Define("gen_eff_eta", f"{gen_eta}[gen_eff_acc]")
                   .Define("gen_eff_matched", f"gen_matched[gen_eff_acc]")
                   .Define("gen_eff_matched_pt", f"gen_matched_pt[gen_eff_acc]"))
            gen_pt, gen_eta, matched = "gen_eff_pt", "gen_eff_eta", "gen_eff_matched"
            matched_pt = "gen_eff_matched_pt"
        cols = [gen_pt, gen_eta, matched, matched_pt]
        data = extract_columns(rdf, cols)
        return {
            "pt": data[gen_pt],
            "eta": data[gen_eta],
            "matched": data[matched] != 0,
            "reco_pt": data[matched_pt],
        }


class ConfigRef(Config):
    def __init__(
        self,
        *,
        samples_path: Optional[str] = None,
        pt_branch: str,
        preprocess_function: Optional[Callable] = None,
        scaling_function: Optional[Callable] = None,
        rate: Optional[np.array] = None,
        tree: str = "Events",
        score_branch: Optional[str] = None,  # Needed only for applying WP
        WP: Optional[tuple[Iterable, Iterable]] = None,
        eff_config: Optional[ConfigEff] = None,
    ):
        super().__init__(
            samples_path=samples_path,
            pt_branch=pt_branch,
            score_branch=score_branch,
            preprocess_function=preprocess_function,
            scaling_function=scaling_function,
            rate=rate,
            tree=tree,
            eff_config=eff_config,
        )
        self.WP = WP


class ConfigObj(Config):
    def __init__(
        self,
        *,
        samples_path: Optional[str] = None,
        pt_branch: str,
        score_branch: str,
        preprocess_function: Optional[Callable] = None,
        scaling_function: Optional[Callable] = None,
        tree: str = "Events",
        refs: Optional[list[str]] = None,
        eff_config: Optional[ConfigEff] = None,
    ):
        super().__init__(
            samples_path=samples_path,
            pt_branch=pt_branch,
            score_branch=score_branch,
            preprocess_function=preprocess_function,
            scaling_function=scaling_function,
            rate=None,
            tree=tree,
            eff_config=eff_config,
        )

        self.refs = refs
        self.records = dict()

    def addToRecord(self, ref, record_name, bins, cuts, rate, cuts_err=None, chi2=None):
        mask = np.bitwise_and(cuts != -np.inf, cuts > -9999.0)
        bins = bins[mask]
        cuts = cuts[mask]
        if cuts_err is not None:
            cuts_err = cuts_err[mask]
            if isinstance(cuts_err, np.ndarray):
                cuts_err = cuts_err.tolist()
        if isinstance(bins, np.ndarray):
            bins = bins.tolist()
        if isinstance(cuts, np.ndarray):
            cuts = cuts.tolist()
        if isinstance(rate, np.ndarray):
            rate = rate.tolist()
        if isinstance(ref.rate, np.ndarray):
            ref_rate = ref.rate.tolist()
        else:
            ref_rate = ref.rate

        if ref.name not in self.records:
            self.records[ref.name] = {
                "ref_rate": ref_rate,
                record_name: {"bins": bins, "cuts": cuts, "rate": rate},
            }
        else:
            self.records[ref.name][record_name] = {
                "bins": bins,
                "cuts": cuts,
                "rate": rate,
            }
        if chi2 is not None:
            self.records[ref.name][record_name]["chi2"] = chi2
        if cuts_err is not None:
            self.records[ref.name][record_name]["cuts_err"] = cuts_err
