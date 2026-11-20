CutFinder finds working points (score thresholds) as a function of pT to match target rates and produces plots and records.

# How it works
The idea of CutFinder is the following
1. You define in a config all the samples for which you want to find a WP (and preprocess them is needed) and all the samples that you want to use as a reference (or provide directly the reference rate)
2. A fine binning is used to compute the base rate (without WP) and to find a score cut in every single rate bin
3. The pt_bins vs score_threshold plot is then fitted with a step function to find a set of WP of just few cuts.

# Quickstart

1. Install dependencies:
You can create a conda/mamba enviroment with
`mamba create -f requirements.txt --prefix .venv`

or, if you have already an enviroment you can install the dependencies with `pip` (Achtung! ROOT is required)

2. Run the CLI with an example config:
```bash
python cutFinder -c configs/example.py -o output
```
Use `--regressor-only` to run only the regression on pre-computed cuts.

## Output

For each `ConfigObj` run the tool produces:
- `output/<obj.name>/records.json` — recorded `full` and `fitted` cuts and rates.
- `output/<obj.name>/rates.png`, `cuts.png` (and PDF variants) — plots saved by [`CutFinder.plots.Plotter`](CutFinder/plots.py).

## Notes

- ROOT is used extensively; the CLI declares the C++ helper via `ROOT.gInterpreter.Declare('#include "include/functions.cpp"')` in [cutFinder](cutFinder).
- If online-to-offline scaling functions are provided they are handled in [`CutFinder.configs.Config.scale`](CutFinder/configs.py) and inverse-mapped when computing thresholds.

# Instructions

- Configuration and global settings: [`CutFinder.configs.GlobalConf`](CutFinder/configs.py), [`CutFinder.configs.ConfigObj`](CutFinder/configs.py), [`CutFinder.configs.ConfigRef`](CutFinder/configs.py)
- Algorithms: [`CutFinder.algorithms`](CutFinder/algorithms.py) — algorithms that find bin-by-bin score cuts.
- Regressor: [`CutFinder.regressors`](CutFinder/regressors.py) — algorithms that fits cuts as piecewise-constant blocks.

## Cli options
The cli command is `cutFinder`

- `-c path/to/config.py` specify the path to the config file
- `-o path/to/outfolder` specify the path in which to save the output folder
- `-j int` Controls how many cores RDataframe have to use (default: ALL OF THEM)
- `--regressor-only` often you need to compute the bin-by-bin cuts only once and then finetune the regressor (unless you need a finer binning). With this command you can use the already computed rates and cuts loading them from the `records.json` located in the previously saved output folder.

## Configs
You can find some examples in the config folder.

The config file is a python file in which you have to specify 3 types of Config objects

### [`CutFinder.configs.GlobalConf`](CutFinder/configs.py)
Here have to specify some global settings
- `pt_bins` the fine pt_binning to use to compute the rate (the online pT thresholds of the obj; see "Offline-to-Online scaling").
- `algo` the algorithm to use to compute the bin-by-bin cuts (defined in CutFinder/algorithms.py) (default = iterative_bin_cutter)
- `regressor` the algorithm that perform the step function fit (defined in CutFinder/regressors.py) (default = bayesian_blocks_gaussian)
- `fitrange` the pT range on which the regressor performs the fit (default = None (full range))
- `algo_kwargs` keyword arguments to pass to the algorithm function
- `regressor_kwargs` keyword arguments to pass to the regressor function

### [`CutFinder.configs.ConfigRef`](CutFinder/configs.py)
These are the reference, i.e. CutFinder will find the WP to reproduce the rate defined in the ConfigRefs

You have 2 options:
1. You can directly define the rate curve to reproduce passing it as an array to the `rate` argument
2. Compute the rate form some sample

In the latter case you have to specify
- `samples_path` where the samples are located. You can use xrootd
- `pt_branch` the name of the pt branch

Additionaly you can specify
- `preprocess_function`: A function that takes a ROOT RDataFrame as input and return another RDF. This is useful to preprocess and manipulate your samples (e.g to add some cuts)
- `WP`. [List, List] A set of WP (pt_bins, score_cuts) to apply to the sample
- `score_branch` name of the score branch. Required only if you use the `WP` argument.

### [`CutFinder.configs.ConfigObj`](CutFinder/configs.py)
These are the objective samples, i.e. the samples that will be used to find the new WP.

The arguments are the same of `ConfigRef` but `score_branch` is mandatory.

You can also add a list of names in the `refs` argument to specify just some reference to be compared with for that configobj.

### Manipulations
Config objects can be cloned and manipulated using the `.clone(new_arg=...)` method (it works like cmssw clone).

You can even switch from a `ConfigObj` to a `ConfigRef` and viceversa with `obj_conf.clone(ConfigRef, new_arg=...)`

A cloned config always starts from scratch: it reopens the samples and re-runs the preprocessing on first use. If the clone keeps the same `samples_path` and `preprocess_function` and the source config was already computed (as in `cutFinder`'s fitted-WP evaluation), use `.clone_reusing_rdf(new_arg=...)` instead: it reuses the already-processed (and cached, see below) RDataFrame, avoiding a second full pass over the samples.

Note: `Config.compute()` materializes the pt (and score) columns with `.Cache()` in memory right after the preprocessing, so all the subsequent event loops (per-pt-bin cut search, WP-rate evaluation) run on cached data instead of re-reading the samples for every loop.

### [`CutFinder.configs.ConfigEff`](CutFinder/configs.py)
Attach an `eff_config=ConfigEff(...)` to a `ConfigObj` (and optionally to the `ConfigRef`) to get gen-level efficiency plots for the WPs that were found. For each obj/ref pair and for each variable (gen pt, gen eta) a plot with 3 curves is produced: reference selection, obj with the per-bin WP ("full") and obj with the regression-fitted WP.

- `samples_path`: the (signal) sample used for the efficiency.
- `reco`: dict with `"pt"`, `"eta"`, `"phi"` branch names (post-preprocessing); `"score"` is needed for the obj WPs.
- `gen`: dict with `"pt"`, `"eta"`, `"phi"` gen branch names.
- `preprocess_function`: selection defining the reco objects that enter the matching. **It must NOT filter events** and must mask the eta/phi branches like the pt branch: events where no reco object survives still contribute their gen electrons to the denominator (as unmatched).
- Optional: `tree` (default `"Events"`), `deltaR` gen–reco matching cone (default 0.1), `gen_acceptance` expression restricting the gen denominator to the geometric acceptance (e.g. `"abs(GenEl_caloeta) < 1.479"` — without it, gen electrons that no reco collection can ever match pollute the efficiency), `pt_bins`, `eta_bins` for the efficiency histograms.
- **CAVEAT**: restrict the eta_bins to the acceptance too when `gen_acceptance` is set, otherwise bins outside are just empty.

The matching is "gen electron matched if any (WP-passing) reco object lies within dR"; the efficiency is `gen[matched] / gen[all]`. Columns are extracted with the ForeachSlot-based `CutFinder.extractor.extract_columns` (vector columns must NOT go through `AsNumpy`, which materializes slow object arrays).

## Offline-to-Online scaling
The rate comparison can be performed on the OFFLINE pT scale: `pt_bins` in `GlobalConf` is the ONLINE pT binning of the obj (the WP thresholds, starting at 0), and the rate scan runs on the corresponding OFFLINE bins (`offline = a * online + b`, starting at the offset of the map), so that obj and ref are compared at equal offline thresholds.

There are two ways to enable it:
1. **Derived**: if and only if both the `ConfigRef` and the `ConfigObj` of every comparison have an `eff_config` (and no explicit `scaling_function`), the online-to-offline scaling of each config is derived from its own efficiency turn-on: for a set of online pT cuts (`scaling_kwargs={"pt_cuts": ...}` in `GlobalConf`, default `np.arange(15, 56, 10)`) the turn-on vs gen pT is fitted, the gen pT at 95% efficiency (the so-called 95% WP) is extracted for every cut and the points are linear-fitted as `offline = a * online + b` (the same recipe as the TkEleRegression rate studies, see `CutFinder/scaling.py`). The fitted turn-ons (`turn_on_<cut>.png`) and the linear fit (`online_to_offline_scaling.png`) are saved in `<output>/scaling/<config name>/`; the derived `(a, b)` are printed and stored in `records.json` under `"scaling"`.
2. **Explicit**: pass an online-to-offline `scaling_function` as a lambda of a sympy symbol to the config objects (e.g. `lambda pt: 1.1 * pt + 2`).

Achtung! All the configs must or must not have all together a scaling function defined.

When the scaling is active the output also contains the same WP on the other pT scale: `rates_vs_online.png` (rate vs the online threshold) and `cuts_vs_offline.png` (cuts vs the offline threshold), next to the usual `rates.png` (offline) and `cuts.png` (online).

## Algorithms
The algorithms are function defined in CutFinder/algorithms.py that have
- Input: (ConfigRef, ConfigObj, GlobConf, **kwargs)
- Output: cuts, cuts_err (np.array, np.array), i.e. the cuts to apply to every single rate pt bin and their uncertainty


### iterative_bin_cutter algorithm
Currently this is the only algorithm defined.

Starting from the last pt bin it finds the cut to apply, apply the cut, remove events which contains objects that already triggered in the processed bins and proceed recursively.

The cuts_err are estimated through Bootstrapping that can be controlled with the `n_bootstrap` argument (default = 1000)

Achtung!: the bootstrap slows the process by a lot, do not use high values

## Regressors
The regressor are function defined in CutFinder/regressors.py that have
- Input: rate ptbins, bin-by-bin cuts, sigma (the cuts_err), fitrange (the range on which to perform the fit), **kwargs
- Output: new pt edges, new cuts, chi2 of the fit

### bayesian_blocks_gaussian regressor
Currently this is the only regressor defined.

It employs the bayesian blocks algotithm using a gaussian likelihood.

The `penalty` parameter controls how much adding a new cut is penalized (higher `penalty`, lower the final number of cuts)

## Output
The output folder contains:
- The python config file
- A directory for each `ConfigObj` defined which contains
  - index.php file to allow plot visualization on webeos
  - records.json that contains all the results (WP, cuts, rates, etc. for every reference)
  - A rate plot that show all the reference plots, all the bin-by-bin computed rated and the rates with the final WP.
  - A pt vs score cut plot that show both the bin-by-bin cuts in function of pt and the fitted ones for every reference


# TODO
- Implement `ConfigEff` and incapsulate it into `ConfigObj` and `ConfigRef` to enable plotting efficiency curve for the reference and on the new objects with the new cuts
- Allow to plot only references with given WP (useful for manual finetuning if needed)