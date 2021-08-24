#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

NamedFunc NminusOne(vector<NamedFunc> cuts, int i) {
  NamedFunc total("1");
  for(int j(0); j < static_cast<int>(cuts.size()); j++) 
    if(j != i)
      total = total && cuts.at(j);
  return total;
}

bool checkBit(int i, int n) {
  return((i%static_cast<int>(pow(2,n+1)))/static_cast<int>(pow(2,n)));
}

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type data = Process::Type::data;
  // string data_path("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string data_path("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/mc/merged_zgmc_llg/");
  // string data_path("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_data/2017/data/merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_data = Process::MakeShared<Baby_pico>("Data",             data, 
                      kBlack,                      {data_path+"*.root"},  trigs);

  proc_data->SetMarkerSize(1);

  vector<shared_ptr<Process>> procs = {proc_data};

  NamedFunc baseline("nphoton > 0 && nll > 0");
  NamedFunc mass_cuts("ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110 && photon_drmin[0] > 0.4");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10"};
  vector<NamedFunc> cut = {"llphoton_m[0] + ll_m[0] >= 185",
                           "llphoton_m[0] > 100 && llphoton_m[0] < 180",
                           "photon_pt[0]/llphoton_m[0] >= 15./110",
                           "photon_drmin[0] > 0.4"};
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(1077)
          .CanvasWidth(900)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return weight*100;
    return weight;
  });
  PlotMaker pm;
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && mass_cuts;
    pm.Push<Hist1D>(Axis(20, 100, 180, "llphoton_m[0]", "m_{ll#gamma}", {}), selection, procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(10, 0, 10, "nphoton", "nphoton", {}), selection, procs, ops).Weight(wgt);
    pm.Push<Hist1D>(Axis(20, -1, 1, "llphoton_cosTheta[0]", "cos(#Theta)", {}), selection, procs, ops).Weight(wgt).Tag("cosTheta"+to_string(i));
    // pm.Push<Hist1D>(Axis(20, -1, 1, "llphoton_costheta[0]", "cos(#theta)", {}), selection, procs, ops).Weight(wgt).Tag("costheta"+to_string(i));
    // pm.Push<Hist1D>(Axis(30, -3, 3, "llphoton_psi[0]", "#phi", {}), selection, procs, ops).Weight(wgt).Tag("phi"+to_string(i));
  }
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

