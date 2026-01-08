from CutFinder.algorithms import iterative_bin_cutter
from CutFinder.functions import applyWP
from CutFinder.regressors import bayesian_blocks_gaussian

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
    ):  # (PU200)
        self.pt_bins = pt_bins
        self.maxRate = maxRate
        self.algo = algo
        self.regressor = regressor
        self.fitrange = fitrange
        self.algo_kwargs = algo_kwargs
        self.regressor_kwargs = regressor_kwargs


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
    ):
        if samples_path is None:
            assert rate is not None

        if scaling_function is not None:
            print(
                "Warning: if scaling_function is provided, pt_bins in GlobalConf will be treated like Offline pT bins. The final thresholds will be inverse-mapped to Online pT."
            )
            print(
                "If both scaling_function and rate are provided, scaling_function will only inverse-map the final thresholds, rate will not be re-computed."
            )
            x = sp.Symbol("$OnlinePt")
            y = sp.Symbol("$OfflinePt")
            scaling_function = scaling_function(x)
            inverse_scaling = sp.solve(sp.Eq(y, scaling_function), x)[-1]
            self.scaling = sp.ccode(scaling_function).replace(
                "$OnlinePt", self.pt_branch
            )
            self.inverse_scaling = sp.lambdify(y, inverse_scaling, "numpy")
        else:
            self.scaling = None
            self.inverse_scaling = None

        self.samples_path = samples_path
        self.pt_branch = pt_branch
        self.score_branch = score_branch

        self.func = preprocess_function
        self.rate = rate
        self.tree = tree

        ####### to compute
        self.name = None
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
            "scaling_function": self.scaling,
            "rate": self.rate,
            "tree": self.tree,
        }
        base_kwargs.update(kwargs)
        filtered = {k: v for k, v in base_kwargs.items() if k in params}
        return cls(**filtered)

    def loadRDF(self):
        if self.rdf is None:
            chain = ROOT.TChain(self.tree)
            if "{" in self.samples_path or "}" in self.samples_path:
                range_part = self.samples_path.split("{")[1].split("}")[0]
                start, end = map(int, range_part.split(".."))
                for i in range(start, end + 1):
                    chain.Add(self.samples_path.replace("{" + range_part + "}", str(i)))
            else:
                chain.Add(self.samples_path)
            self.rdf = ROOT.RDataFrame(chain)
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
            self.rdf = self.rdf.Redefine(self.pt_branch, self.scaling)
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
                    ("", "", len(bins) - 1, array("f", bins)), "MAXPT"
                )
            else:
                raise ValueError("bins must be either a tuple or a numpy array.")

            th = th.GetValue()
            h = hist.Hist(axis, storage=hist.storage.Weight())

            for i in range(0, th.GetNbinsX() + 2):
                h.fill(th.GetBinLowEdge(i), weight=th.GetBinContent(i))

            h_scaled = (self.nEvents / self.TotEvents) * maxRate * h / h.integrate(0, 0).value

            rate = np.array([])
            rate_err = np.array([])
            for idx, _ in enumerate(bins):
                integral = h_scaled.integrate(0, idx)
                rate = np.append(rate, integral.value)
                rate_err = np.append(rate_err, integral.variance**0.5)
            self.rate = rate
            self.rate_err = rate_err
            return h

    def compute(self):
        if not self.isComputed:
            if self.name is not None:
                pprint(
                    f"[bold green]Computing {self.__class__.__name__}:[/bold green]\n\t{self.name}\n"
                )
            self.loadRDF()
            self.runPreprocess()
            self.apply_WP()
            self.rdf = self.rdf.Filter(f"{self.pt_branch}.size()>0")
            self.scale()
            self.getEntries()
            self.isComputed = True

            # self.makeRate(pt_bins, maxRate)

    def apply_WP(self):
        if self.WP is not None and not self.isWPApplied:
            self.rdf = applyWP(
                self.pt_branch, self.score_branch, self.WP[0], self.WP[1], self.rdf
            )
            self.isWPApplied = True


class ConfigEff:
    pass
    # implement at the end to plot also efficiencies with the WP that were found


class ConfigRef(Config):
    def __init__(
        self,
        *,
        samples_path: Optional[str] = None,
        pt_branch: str,
        preprocess_function: Optional[Callable] = None,
        scaling_function: Optional[str] = None,
        rate: Optional[np.array] = None,
        tree: str = "Events",
        score_branch: Optional[str] = None,  # Needed only for applying WP
        WP: Optional[tuple[Iterable, Iterable]] = None,
    ):
        super().__init__(
            samples_path=samples_path,
            pt_branch=pt_branch,
            score_branch=score_branch,
            preprocess_function=preprocess_function,
            scaling_function=scaling_function,
            rate=rate,
            tree=tree,
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
        scaling_function: Optional[str] = None,
        tree: str = "Events",
        refs: Optional[list[ConfigRef]] = None,
    ):
        super().__init__(
            samples_path=samples_path,
            pt_branch=pt_branch,
            score_branch=score_branch,
            preprocess_function=preprocess_function,
            scaling_function=scaling_function,
            rate=None,
            tree=tree,
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
