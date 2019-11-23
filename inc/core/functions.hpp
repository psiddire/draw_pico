#ifndef H_FUNCTIONS
#define H_FUNCTIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
 
namespace Functions{
  
  extern const NamedFunc lowDphiFix;
  extern const NamedFunc boostedRegionIdx;
  extern const NamedFunc ntrub;

  bool IsGoodJet(const Baby &b, std::size_t ijet);

  // these are used for doing systematic variations of weights
  enum class Variation{central, up, down};
  NamedFunc MismeasurementCorrection(const std::string &file_path,
                                     const std::string &mismeas_scenario,
                                     Variation variation = Variation::central);
  NamedFunc MismeasurementWeight(const std::string &file_path,
                                 const std::string &mismeas_scenario);
}

#endif
