#ifndef H_ZG_UTILITIES
#define H_ZG_UTILITIES

#include <iostream>
#include <string>
#include <utility>
#include <set>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <stdlib.h>
#include <regex>

#include "TString.h"
#include "TLorentzVector.h"

#include "core/named_func.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/plot_maker.hpp"
#include "core/table.hpp"
#include "core/table_row.hpp"
#include "core/gamma_params.hpp"
#include "core/utilities.hpp"
#include "core/palette.hpp"
#include "core/baby.hpp"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooGenericPdf.h"

// #include "zgamma/KinZfitter.h"

namespace ZgUtilities {
  TLorentzVector AssignL1(const Baby &b, bool gen = false);
  TLorentzVector AssignL2(const Baby &b, bool gen = false);
  TLorentzVector AssignZ (const Baby &b, bool gen = false);
  TLorentzVector AssignH (const Baby &b, bool gen = false);
  TLorentzVector AssignQ1(const Baby &b, bool gen = false);
  TLorentzVector AssignQ2(const Baby &b, bool gen = false);
  TLorentzVector AssignGamma (const Baby &b, bool gen = false);
  TLorentzVector AssignHGG (const Baby &b);
  double lambdaZ(const Baby &b, bool gen = false);
  double cos_theta(const Baby &b, bool gen = false);
  double cos_Theta(const Baby &b, bool gen = false);
  double Getphi(const Baby &b, bool gen = false);
  double KinRefit(const Baby &b);
  std::vector<TLorentzVector> RefitP4(const Baby &b);
  double AssignL1Error(const Baby &b);
  double AssignL2Error(const Baby &b);
}
#endif
