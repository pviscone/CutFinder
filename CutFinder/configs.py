from CutFinder.functions import applyWP

from typing import Optional
from collections.abc import Callable, Iterable

from array import array
import hist
import numpy as np
import ROOT
import sympy as sp


class GlobalConf:
    def __init__(self, pt_bins, maxRate = 31038.96):  #(PU200)
        self.pt_bins = pt_bins
        self.maxRate = maxRate


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

        self.to_compute = True if rate is None else False
        self.to_preprocess = False if preprocess_function is None else True

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
        chain = ROOT.TChain(self.tree)
        chain.Add(self.samples_path)
        self.rdf = ROOT.RDataFrame(chain)
        self.TotEvents = self.rdf.Count().GetValue()

    def runPreprocess(self):
        if self.func is not None:
            self.rdf = self.func(self.rdf)
            self.nEvents = self.rdf.Count().GetValue()

    def getEntries(self):
        self.nEvents = self.rdf.Count().GetValue()

    def scale(self):
        if self.scaling is not None:
            self.rdf = self.rdf.Redefine(self.pt_branch, self.scaling)

    def makeRate(self, bins, maxRate):
        self.rdf = self.rdf.Filter(f"{self.pt_branch}.size()>0").Redefine(
            self.pt_branch, f"Max({self.pt_branch})"
        )
        self.scale()
        self.getEntries()

        if isinstance(bins, tuple):
            axis = hist.axis.Regular(bins[0], bins[1], bins[2])
            th = self.rdf.Histo1D(("", "", bins[0], bins[1], bins[2]), self.pt_branch)
        elif isinstance(bins, np.ndarray):
            axis = hist.axis.Variable(bins)
            th = self.rdf.Histo1D(
                ("", "", len(bins) - 1, array("f", bins)), self.pt_branch
            )
        else:
            raise ValueError("bins must be either a tuple or a numpy array.")

        th = th.GetValue()
        h = hist.Hist(axis, storage=hist.storage.Weight())

        for i in range(1, th.GetNbinsX() + 1):
            h.fill(th.GetBinLowEdge(i), weight=th.GetBinContent(i))

        h = (self.nEvents / self.TotEvents) * maxRate * h / h.integrate(0, 0).value

        rate = np.array([])
        rate_err = np.array([])
        for idx, _ in enumerate(bins[:-1]):
            integral = h.integrate(0, idx)
            rate = np.append(rate, integral.value)
            rate_err = np.append(rate_err, integral.variance**0.5)
        self.rate = rate
        self.rate_err = rate_err

    def compute(self, pt_bins, maxRate):
        if self.rate is None:
            self.loadRDF()
            self.runPreprocess()
            self.apply_WP()
            self.makeRate(pt_bins, maxRate)

    def apply_WP(self):
        if self.WP is not None:
            self.rdf = applyWP(
                self.pt_branch, self.score_branch, self.WP[0], self.WP[1], self.rdf
            )


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
        self.rate_after_cuts = None
        self.rate_after_cuts_err = None
