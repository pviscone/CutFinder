from CutFinder.configs import ConfigObj, ConfigRef, GlobalConf

import numpy as np

glob = GlobalConf(pt_bins=np.linspace(0, 100, 101))

new_id = ConfigObj(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/151X_TkElePtRegrTemp_A0/*.root",
    pt_branch="TkEleL2_pt",
    score_branch="TkEleL2_idScore",
)

elliptic = ConfigRef(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/140Xv0C1/*.root",
    pt_branch="TkEleL2_pt",
)


def tightID(rdf):
    return rdf.Redefine("TkEleL2_pt", "TkEleL2_pt[(TkEleL2_hwQual & 2)==2]").Redefine(
        "TkEleL2_idScore", "TkEleL2_idScore[(TkEleL2_hwQual & 2)==2]"
    )
elliptic_tight = elliptic.clone(preprocess_function=tightID)
custom = new_id.clone(ConfigRef, WP=[[0], [0]])


#!Create refs argument in configObj ????