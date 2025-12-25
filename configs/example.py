from CutFinder.configs import ConfigObj, ConfigRef, GlobalConf

import numpy as np

glob = GlobalConf(pt_bins=np.arange(0, 65, 5))

new_id = ConfigObj(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/151X_TkElePtRegrTemp_A0/perfNano_10044360_{95..100}.root",
    pt_branch="TkEleL2_pt",
    score_branch="TkEleL2_idScore",
)
custom = new_id.clone(ConfigRef, WP=[[0], [0]])

"""
elliptic = ConfigRef(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/140Xv0C1/perfNano_8843496_{95..100}.root",
    pt_branch="TkEleL2_pt",
)


def tightID(rdf):
    return rdf.Redefine("TkEleL2_pt", "TkEleL2_pt[(TkEleL2_hwQual & 2)==2]").Redefine(
        "TkEleL2_idScore", "TkEleL2_idScore[(TkEleL2_hwQual & 2)==2]"
    )


elliptic_tight = elliptic.clone(preprocess_function=tightID)
"""
#!Add in the example refs argument in configObj
