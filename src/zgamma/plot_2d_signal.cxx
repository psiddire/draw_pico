#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist2d.hpp"
#include "core/utilities.hpp"
using namespace std;
using namespace PlotOptTypes;

int main() {
  gErrorIgnoreLevel = 6000;
  Palette col("txt/colors.txt","default");

  Process::Type sig  =  Process::Type::background;
  string sig_path("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);

  auto proc_hzg = Process::MakeShared<Baby_pico>("HToZg", sig, kRed, {sig_path+"*.root"}, trigs);
  // proc_hzg->SetMarkerSize(0.25);

  NamedFunc baseline("nphoton > 0 && nll > 0");
  NamedFunc lep("(ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15) || "
		"ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10");
  NamedFunc pho("photon_pt[0] > 15 && photon_drmin[0] > 0.4 && photon_sig[0]");
  NamedFunc mass_cuts("ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110");

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight*100;
  });

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> twodim_plotopts = {style().Title(TitleType::info).Overflow(OverflowType::overflow)};

  // PlotOpt style("txt/plot_styles.txt", "CMSPaper");
  // vector<PlotOpt> twodim_plotopts = {style().Title(TitleType::info).Stack(StackType::data_norm)};
  // vector<PlotOpt> twodim_plotopts = {style().Title(TitleType::info).Bottom(BottomType::off).YAxis(YAxisType::linear).Stack(StackType::data_norm).LegendColumns(3)};

  PlotMaker pm;
  vector<shared_ptr<Process>> procs = {proc_hzg};
  for(int i = 0; i < 1; i++) {
    pm.Push<Hist2D>(Axis(20, 50, 150, "ll_m[0]", "m_{ll} [GeV]", {}),
		    Axis(20, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}),
                    baseline && lep && pho && mass_cuts, procs, twodim_plotopts).Tag("Mllg_Mll").Weight(wgt);
    pm.Push<Hist2D>(Axis(20, -1.0, 1.0, "photon_idmva[0]", "Photon ID MVA", {}),
		    Axis(20, 0.00, 0.20, "photon_pterr[0]/photon_pt[0]", "#sigma_{#gamma}", {}),
                    baseline && lep && pho && mass_cuts, procs, twodim_plotopts).Tag("mva_sigma").Weight(wgt);
  }
  // baseline && lep && pho && mass_cuts

  pm.min_print_ = true;
  pm.MakePlots(41.5);
}
