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

    for (size_t j = pt_bins.size(); j-- > 0;) {
      if (p >= pt_bins[j] && s >= score_cuts[j]) {
        mask[i] = true;
        break;
      }
    }
  }
  return mask;
}

#endif // !FUNCTIONS_CPP