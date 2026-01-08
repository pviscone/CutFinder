from CutFinder.configs import ConfigObj, ConfigRef, GlobalConf

import numpy as np

glob = GlobalConf(pt_bins=np.arange(0, 65, 5))

new_id = ConfigObj(
    samples_path="root://eoscms.cern.ch//eos/cms/store/cmst3/group/l1tr/pviscone/l1teg/fp_ntuples/NuGunAllEta_PU200/FP/151X_TkElePtRegrTemp_A0/perfNano_10044360_{95..100}.root",
    pt_branch="TkEleL2_pt",
    score_branch="TkEleL2_idScore",
)

# Same sample but preappling a custom WP (idScore>0)
custom_ref = new_id.clone(ConfigRef, WP=[[0], [0]])
