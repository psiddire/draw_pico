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

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type sig  = Process::Type::signal;
  string bfolder("/net/cms17/cms17r0/pico/");
  string ul_path(bfolder+"NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");

  NamedFunc mass_err("mass_err", [](const Baby &b) -> NamedFunc::ScalarType{
      double dml1 = b.llphoton_dml1()->at(0);                                                                                                                             
      double dml2 = b.llphoton_dml2()->at(0);                                                                                                                             
      double dmph = b.llphoton_dmph()->at(0);                                                                                                                             
      return sqrt(dml1*dml1 + dml2*dml2 + dmph*dmph);                                                                                                                     
    });
  NamedFunc mass_res("mass_res", [](const Baby &b) -> NamedFunc::ScalarType{                                                                                              
      double dml1 = b.llphoton_dml1()->at(0);                                                                                                                             
      double dml2 = b.llphoton_dml2()->at(0);                                                                                                                             
      double dmph = b.llphoton_dmph()->at(0);                                                                                                                             
      double mass = b.llphoton_m()->at(0);                                                                                                                                
      return sqrt(dml1*dml1 + dml2*dml2 + dmph*dmph) / mass;                                                                                                              
    });                                                                                                                                                                   

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
      double weight = b.w_lumi();
      return weight*100;
    });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc mass_cuts("ll_m[llphoton_ill[0]] > 50 && "
		      "llphoton_m[0]+ll_m[llphoton_ill[0]]>=185 && "
		      "llphoton_m[0] > 100 && llphoton_m[0] < 180 && "
		      "photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110");
  NamedFunc baseline("nphoton > 0 && nll > 0");
  vector<NamedFunc> lep = {"ll_lepid[llphoton_ill[0]] == 11 && "
			   "el_pt[ll_i1[llphoton_ill[0]]] > 25 && "
			   "el_pt[ll_i2[llphoton_ill[0]]] > 15 && "
			   "el_sig[ll_i1[llphoton_ill[0]]] && "
			   "el_sig[ll_i2[llphoton_ill[0]]]",
                           "ll_lepid[llphoton_ill[0]] == 13 && "
			   "mu_pt[ll_i1[llphoton_ill[0]]] > 20 && "
			   "mu_pt[ll_i2[llphoton_ill[0]]] > 10 && "
			   "mu_sig[ll_i1[llphoton_ill[0]]] && "
			   "mu_sig[ll_i2[llphoton_ill[0]]]"};
  NamedFunc pho("photon_pt[llphoton_iph[0]] > 15 && "
		"photon_drmin[llphoton_iph[0]] > 0.4 && "
		"photon_sig[llphoton_iph[0]]");

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};

  auto proc_UL     = Process::MakeShared<Baby_pico>("HToZ#gamma (UL)", sig, 
						     TColor::GetColor("#ff0000"), {ul_path+"*.root"}, trigs);

  proc_UL->SetLineWidth(3);

  vector<shared_ptr<Process>> procs_UL = {proc_UL};

  PlotMaker pm;
  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && pho && mass_cuts;
    pm.Push<Hist1D>(Axis(100, 0, 5, mass_err, "#sigma_{ll#gamma} [GeV]", {}), selection, procs_UL, ops).Weight(wgt).Tag(lepName[i]+"llg_err");
    pm.Push<Hist1D>(Axis(100, 0.0, 0.05, mass_res, "#sigma_{ll#gamma}/m_{ll#gamma}", {}), selection, procs_UL, ops).Weight(wgt).Tag(lepName[i]+"llg_res");
  }
  pm.min_print_ = true;
  pm.multithreaded_ = false;
  pm.MakePlots(41.5);
}

