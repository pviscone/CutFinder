import numpy as np

#!All the function are supposed to be used also with functools.partial to be used as preprocess_function in Configs.
#!I.E. RDF must be the last argument of the function and they must return the modified RDF.


def applyWP(pt_branch, score_branch, pt_bins, score_thresholds, rdf):
    pt_bins = np.array(pt_bins).tolist()
    score_thresholds = np.array(score_thresholds).tolist()
    pt_bins = "{" + ", ".join(map(str, pt_bins)) + "}"
    score_thresholds = "{" + ", ".join(map(str, score_thresholds)) + "}"
    rdf = (
        rdf.Filter(f"{pt_branch}.size()>0")
        .Define(
            "WPmask",
            f"WP_mask({pt_branch}, {score_branch}, {pt_bins}, {score_thresholds})",
        )
        .Redefine(pt_branch, f"{pt_branch}[WPmask]")
        .Redefine(score_branch, f"{score_branch}[WPmask]")
    )
    return rdf
