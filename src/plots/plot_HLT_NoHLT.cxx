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
  string bfolder("/net/cms29/cms29r0/pico/");
  string sig_path(bfolder+"NanoAODv2/zgamma_signal_l1/2017/signal/merged_zgmc_llg/");

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight*100;
  });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc mass_cuts("ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110");
  NamedFunc baseline("nphoton > 0 && nll > 0");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15 && "
			   "el_sig[ll_i1[0]] && el_sig[ll_i2[0]]",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10 && "
			   "mu_sig[ll_i1[0]] && mu_sig[ll_i2[0]]"};
  NamedFunc pho("photon_pt[0] > 15 && photon_drmin[0] > 0.4 && photon_sig[0]");

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

  auto proc_HLT     = Process::MakeShared<Baby_pico>("HToZ#gamma (HLT)", sig, 
                      TColor::GetColor("#ff0000"), {sig_path+"*.root"}, trigs);
  auto proc_NoHLT  = Process::MakeShared<Baby_pico>("HToZ#gamma (No HLT)", sig, 
                      TColor::GetColor("#0000ff"), {sig_path+"*.root"}, !trigs);

  proc_HLT->SetLineWidth(3);
  proc_NoHLT->SetLineWidth(3);

  vector<shared_ptr<Process>> procs = {proc_HLT, proc_NoHLT};

  PlotMaker pm;
  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && pho && mass_cuts;
    if(i == 0){
      pm.Push<Hist1D>(Axis(15, 0, 150, "el_pt[ll_i1[0]]", "Leading electron p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "el_eta[ll_i1[0]]", "Leading electron #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "el_phi[ll_i1[0]]", "Leading electron #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_phi");
      pm.Push<Hist1D>(Axis(15, 0, 150, "el_pt[ll_i2[0]]", "Trailing electron p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "el_eta[ll_i2[0]]", "Trailing electron #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "el_phi[ll_i2[0]]", "Trailing electron #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_phi");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG10", "L1_SingleEG10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG10");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG15", "L1_SingleEG15", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG15");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG18", "L1_SingleEG18", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG18");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG24", "L1_SingleEG24", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG24");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG26", "L1_SingleEG26", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG26");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG28", "L1_SingleEG28", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG28");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG30", "L1_SingleEG30", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG30");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG32", "L1_SingleEG32", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG32");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG34", "L1_SingleEG34", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG34");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG34er2p1", "L1_SingleEG34er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG34er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG36", "L1_SingleEG36", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG36");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG36er2p1", "L1_SingleEG36er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG36er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG38", "L1_SingleEG38", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG38");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG38er2p1", "L1_SingleEG38er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG38er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG40", "L1_SingleEG40", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG40");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG42", "L1_SingleEG42", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG42");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG45", "L1_SingleEG45", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG45");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleEG50", "L1_SingleEG50", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleEG50");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG18", "L1_SingleIsoEG18", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG18");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG18er2p1", "L1_SingleIsoEG18er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG18er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG20", "L1_SingleIsoEG20", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG20");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG20er2p1", "L1_SingleIsoEG20er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG20er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG22", "L1_SingleIsoEG22", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG22");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG22er2p1", "L1_SingleIsoEG22er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG22er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG24", "L1_SingleIsoEG24", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG24");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG24er2p1", "L1_SingleIsoEG24er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG24er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG26", "L1_SingleIsoEG26", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG26");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG26er2p1", "L1_SingleIsoEG26er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG26er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG28", "L1_SingleIsoEG28", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG28");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG28er2p1", "L1_SingleIsoEG28er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG28er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG30", "L1_SingleIsoEG30", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG30");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG30er2p1", "L1_SingleIsoEG30er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG30er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG32", "L1_SingleIsoEG32", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG32");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG32er2p1", "L1_SingleIsoEG32er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG32er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG33er2p1", "L1_SingleIsoEG33er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG33er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG34", "L1_SingleIsoEG34", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG34");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG34er2p1", "L1_SingleIsoEG34er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG34er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG35", "L1_SingleIsoEG35", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG35");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG35er2p1", "L1_SingleIsoEG35er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG35er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG36", "L1_SingleIsoEG36", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG36");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG36er2p1", "L1_SingleIsoEG36er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG36er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG37", "L1_SingleIsoEG37", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG37");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG38", "L1_SingleIsoEG38", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG38");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG38er2p1", "L1_SingleIsoEG38er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG38er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG40", "L1_SingleIsoEG40", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG40");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_SingleIsoEG40er2p1", "L1_SingleIsoEG40er2p1", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_SingleIsoEG40er2p1");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_15_10", "L1_DoubleEG_15_10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_15_10");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_18_17", "L1_DoubleEG_18_17", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_18_17");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_20_18", "L1_DoubleEG_20_18", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_20_18");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_22_10", "L1_DoubleEG_22_10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_22_10");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_22_12", "L1_DoubleEG_22_12", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_22_12");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_22_15", "L1_DoubleEG_22_15", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_22_15");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_23_10", "L1_DoubleEG_23_10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_23_10");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_24_17", "L1_DoubleEG_24_17", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_24_17");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_25_12", "L1_DoubleEG_25_12", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_25_12");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_25_13", "L1_DoubleEG_25_13", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_25_13");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_25_14", "L1_DoubleEG_25_14", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_25_14");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_LooseIso23_10", "L1_DoubleEG_LooseIso23_10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_LooseIso23_10");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleEG_LooseIso24_10", "L1_DoubleEG_LooseIso24_10", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleEG_LooseIso24_10");
    }
    else if(i == 1){
      pm.Push<Hist1D>(Axis(15, 0, 150, "mu_pt[ll_i1[0]]", "Leading muon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "mu_eta[ll_i1[0]]", "Leading muon #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "mu_phi[ll_i1[0]]", "Leading muon #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_phi");
      pm.Push<Hist1D>(Axis(15, 0, 150, "mu_pt[ll_i2[0]]", "Trailing muon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "mu_eta[ll_i2[0]]", "Trailing muon #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "mu_phi[ll_i2[0]]", "Trailing muon #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_phi");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleMu_12_5", "L1_DoubleMu_12_5", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleMu_12_5");
      pm.Push<Hist1D>(Axis(2, 0, 2, "L1_DoubleMu_15_7", "L1_DoubleMu_15_7", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"L1_DoubleMu_15_7");
    }
  }
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

