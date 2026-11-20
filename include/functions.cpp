#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include <ROOT/RVec.hxx>
#include <vector>

using namespace ROOT;
using namespace ROOT::VecOps;

RVec<bool> WP_mask(const RVecF &pt, const RVecF &score,
                   std::vector<float> pt_bins, std::vector<float> score_cuts) {
  RVec<bool> mask(pt.size(), false);
  for (size_t i = 0; i < pt.size(); ++i) {
    float p = pt[i];
    float s = score[i];

    //Assume pt_bins sorted ascendingly
    for (size_t j = pt_bins.size(); j-- > 0;) {
      if (p >= pt_bins[j]) {
        if (s >= score_cuts[j]){
            mask[i] = true;
        }
        break;
      }
    }
  }
  return mask;
}

RVec<bool> gen_reco_match(const RVecF &gen_eta, const RVecF &gen_phi,
                          const RVecF &reco_eta, const RVecF &reco_phi,
                          double deltaR) {
  // Per gen-level object: is there at least one reco object within deltaR?
  RVec<bool> matched(gen_eta.size(), false);
  for (size_t g = 0; g < gen_eta.size(); ++g) {
    for (size_t r = 0; r < reco_eta.size(); ++r) {
      if (DeltaR(gen_eta[g], reco_eta[r], gen_phi[g], reco_phi[r]) < deltaR) {
        matched[g] = true;
        break;
      }
    }
  }
  return matched;
}

RVecF gen_reco_match_pt(const RVecF &gen_eta, const RVecF &gen_phi,
                        const RVecF &reco_pt, const RVecF &reco_eta,
                        const RVecF &reco_phi, double deltaR) {
  // Per gen-level object: pT of the highest-pT reco object within deltaR
  // (0 when nothing matches).
  RVecF matched_pt(gen_eta.size(), 0.f);
  for (size_t g = 0; g < gen_eta.size(); ++g) {
    float best = 0.f;
    for (size_t r = 0; r < reco_eta.size(); ++r) {
      if (DeltaR(gen_eta[g], reco_eta[r], gen_phi[g], reco_phi[r]) < deltaR &&
          reco_pt[r] > best) {
        best = reco_pt[r];
      }
    }
    matched_pt[g] = best;
  }
  return matched_pt;
}

#endif // !FUNCTIONS_CPP