import numpy as np

# ! NB if scaling is applied, the cuts are computed on the offline pT scale.
#! Inverse scaling will be applied in the regressor to map back to online pT.


def iterative_bin_cutter(ref, obj, glob, n_bootstrap=1000):
    ref.makeRate(glob.pt_bins, glob.maxRate)
    cuts = []
    cuts_err = []

    if ref.rate[-1] == 0.0:
        print(
            f"Warning: target rate in bin -1 ({glob.pt_bins[-1]} GeV) is zero (probably due to low stat). No cut will be applied."
        )
        cuts.append(-np.inf)
        cuts_err.append(0.0)
    else:
        # Select max score in the last pt bin
        scores = (
            obj.rdf.Define(
                "scores", f"{obj.score_branch}[{obj.pt_branch} >= {glob.pt_bins[-1]}]"
            )
            .Filter("scores.size() > 0")
            .Define("max_score", "Max(scores)")
        ).AsNumpy(["max_score"])["max_score"]
        rate_bin = len(scores) * (glob.maxRate / obj.TotEvents)

        # Fraction of events to keep in the last pt bin
        f = ref.rate[-1] / rate_bin

        # check if f>=1 otherwise no cut can be applied
        if f > 1.0:
            print(
                f"Warning: target rate in bin -1 ({glob.pt_bins[-1]} GeV) is higher than current rate. No cut will be applied."
            )
            cuts.append(-np.inf)
        elif len(scores) == 0:
            print(
                f"Warning: no events in the last pt bin {glob.pt_bins[-1]} GeV. No cut will be applied."
            )
            cuts.append(-np.inf)
            cuts_err.append(0.0)
        else:
            # find cut to achieve target rate in the last pt bin
            q = np.quantile(scores, 1 - f)
            if len(scores) < 2 or n_bootstrap <= 0:
                sd = 0.0
            else:
                idx = np.random.randint(0, len(scores), size=(n_bootstrap, len(scores)))
                boot_q = np.quantile(scores[idx], 1 - f, axis=1)
                sd = float(np.std(boot_q, ddof=1))
            cuts.append(q)
            cuts_err.append(sd)

    # Apply cut in the last bin
    if cuts[-1] != -np.inf:
        cut_mask = f"({obj.pt_branch} >= {glob.pt_bins[-1]} && {obj.score_branch} >= {cuts[-1]}) || ({obj.pt_branch} < {glob.pt_bins[-1]})"
        rdf = (
            obj.rdf.Define("cut_mask", cut_mask)
            .Redefine(
                obj.score_branch,
                f"{obj.score_branch}[cut_mask]",
            )
            .Redefine(
                obj.pt_branch,
                f"{obj.pt_branch}[cut_mask]",
            )
            .Filter(f"{obj.pt_branch}.size() > 0")
        )
    else:
        rdf = obj.rdf.Define(
            "cut_mask", f"ROOT::VecOps::RVec<bool>({obj.pt_branch}.size(), true)"
        )  # keep all events

    for i in range(len(glob.pt_bins) - 2, -1, -1):
        print(
            f"Processing pt bin {i}: {glob.pt_bins[i]} - {glob.pt_bins[i + 1]} GeV (Obj: {obj.name}, Ref: {ref.name})"
        )
        # Select max score in the current pt bin (from the last to the first) and discard events which already triggered the previous cuts
        scores = (
            # Discard events which already triggered previous cuts
            rdf.Filter(f"Max({obj.pt_branch}) < {glob.pt_bins[i + 1]}")
            # select objs in current pt bin
            .Redefine(
                obj.score_branch,
                f"{obj.score_branch}[{obj.pt_branch} >= {glob.pt_bins[i]} && {obj.pt_branch} < {glob.pt_bins[i + 1]}]",
            )
            .Redefine(
                obj.pt_branch,
                f"{obj.pt_branch}[{obj.pt_branch} >= {glob.pt_bins[i]} && {obj.pt_branch} < {glob.pt_bins[i + 1]}]",
            )
            .Filter(f"{obj.pt_branch}.size() > 0")
            # Select max score
            .Define("max_score", f"Max({obj.score_branch})")
            .AsNumpy(["max_score"])["max_score"]
        )

        rate_bin = len(scores) * (glob.maxRate / obj.TotEvents) + ref.rate[i + 1]
        f = (ref.rate[i] - ref.rate[i + 1]) / (rate_bin - ref.rate[i + 1])
        if f > 1.0:
            print(
                f"Warning: target rate in bin {i} ({glob.pt_bins[i]} GeV) is higher than current rate. No cut will be applied."
            )
            cuts.append(-np.inf)
            cuts_err.append(0.0)
        elif f < 0.0:
            print(
                f"Warning: target rate in bin {i} ({glob.pt_bins[i]} GeV) is lower than previous rate. No cut will be applied."
            )
            cuts.append(-np.inf)
            cuts_err.append(0.0)
        else:
            # find cut to achieve target rate in the current pt bin
            if len(scores) == 0:
                print(
                    f"Warning: no events in the current pt bin {i} ({glob.pt_bins[i]} GeV). No cut will be applied."
                )
                cuts.append(-np.inf)
                cuts_err.append(0.0)
            else:
                q = np.quantile(scores, 1 - f)
                if len(scores) < 2 or n_bootstrap <= 0:
                    sd = 0.0
                else:
                    idx = np.random.randint(
                        0, len(scores), size=(n_bootstrap, len(scores))
                    )
                    boot_q = np.quantile(scores[idx], 1 - f, axis=1)
                    sd = float(np.std(boot_q, ddof=1))
                cuts.append(q)
                cuts_err.append(sd)

        if (
            cuts[-1] != -np.inf and i > 0
        ):  # it's useless to apply cut in the first bin (last iteration)
            cut_mask = f"({obj.pt_branch} >= {glob.pt_bins[i]} && {obj.pt_branch} < {glob.pt_bins[i + 1]} && {obj.score_branch} >= {cuts[-1]}) || ({obj.pt_branch} < {glob.pt_bins[i]})"
            rdf = (
                rdf.Redefine("cut_mask", cut_mask)
                .Redefine(
                    obj.score_branch,
                    f"{obj.score_branch}[cut_mask]",
                )
                .Redefine(
                    obj.pt_branch,
                    f"{obj.pt_branch}[cut_mask]",
                )
                .Filter(f"{obj.pt_branch}.size() > 0")
            )
        print(f"\tCut found : {cuts[-1]} +- {cuts_err[-1]}\n", flush=True)

    cuts = np.array(cuts[::-1])
    cuts_err = np.array(cuts_err[::-1])
    return cuts, cuts_err
