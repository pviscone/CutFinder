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
- `--regressor-only` often you need to compute the bin-by-bin cuts only once and then finetune the regressor (unless you need a finer binning). With this command you can use the already computed rates and cuts loading them from the `records.json` located in the previously saved output folder.

## Configs
You can find some examples in the config folder.

The config file is a python file in which you have to specify 3 types of Config objects

### [`CutFinder.configs.GlobalConf`](CutFinder/configs.py)
Here have to specify some global settings
- `pt_bins` the fine pt_binning to use to compute the rate.
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

## Offline-to-Online scaling
You can pass online-to-offline scaling function as lambda function to config objects.

Achtung! All the configs must or must not have all together a scaling function defined.

The scaling function will be used to scale online pt to offline pt when plotting rates and to inverse scale WP pt edges inverting the function wiht sympy

## Algorithms
The algorithms are function defined in CutFinder/algorithms.py that have
- Input: (ConfigRef, ConfigObj, GlobConf, **kwargs)
- Output: cuts, cuts_err (np.array, np.array), i.e. the cuts to apply to every single rate pt bin and their uncertainty


### iterative_bin_cutter algorithm
Currently this is the only algorithm defined.

Starting from the last pt bin it finds the cut to apply, apply the cut, remove events which contains objects that already triggered in the processed bins and proceed recursively.

The cuts_err are estimated through Bootstrapping that can be controlled with the `n_bootstrap` argument (default = 1000)

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