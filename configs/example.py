from CutFinder.configs import ConfigObj, ConfigRef, GlobalConf

import numpy as np

glob = GlobalConf(
    pt_bins=np.concatenate([np.arange(0, 20, 2), np.arange(20, 65, 5)]),
    fitrange=(0, 60),
    algo_kwargs={
        "n_bootstrap": 500,
    },
    regressor_kwargs={"penalty": 15.0},
)

CompositeID = ConfigObj(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/151X_TkElePtRegrTemp_A0/perfNano_10044360_*.root",
    pt_branch="TkEleL2_pt",
    score_branch="TkEleL2_idScore",
)

EllipticID = ConfigRef(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/140Xv0C1/perfNano_8843496_*.root",
    pt_branch="TkEleL2_pt",
)


def tightID(rdf):
    return (
        rdf.Redefine("TkEleL2_pt", "TkEleL2_pt[(TkEleL2_hwQual & 2)==2]")
        .Redefine("TkEleL2_idScore", "TkEleL2_idScore[(TkEleL2_hwQual & 2)==2]")
        .Filter("TkEleL2_pt.size()>0")
    )


EllipticID_Tight = EllipticID.clone(preprocess_function=tightID)
