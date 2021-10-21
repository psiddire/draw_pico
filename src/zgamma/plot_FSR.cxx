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
  string sig_path(bfolder+"NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");

  NamedFunc llg_mass("llg_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2, g, llg;
    l1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), 0.105);
    l2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), 0.105);
    g.SetPtEtaPhiM(b.photon_pt()->at(b.llphoton_iph()->at(0)), b.photon_eta()->at(b.llphoton_ill()->at(0)), b.photon_phi()->at(b.llphoton_ill()->at(0)), 0.0);
    llg = l1 + l2 + g;
    return llg.M();
  });

  NamedFunc ll_mass("ll_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2, ll;
    l1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), 0.105);
    l2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), 0.105);
    ll = l1 + l2;
    return ll.M();
  });

  NamedFunc llg_corr_mass("llg_corr_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2, g, fsrph, llg;
    l1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), 0.105);
    l2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), 0.105);
    fsrph.SetPtEtaPhiM(b.fsrphoton_pt()->at(0), b.fsrphoton_eta()->at(0), b.fsrphoton_phi()->at(0), 0.0);
    if (b.fsrphoton_muonidx()->at(0)==b.ll_i1()->at(b.llphoton_ill()->at(0)) && l1.DeltaR(fsrph) < 0.4)
      l1 = l1 + fsrph;
    else if (b.fsrphoton_muonidx()->at(0)==b.ll_i2()->at(b.llphoton_ill()->at(0)) && l2.DeltaR(fsrph) < 0.4)
      l2 = l2 + fsrph;
    g.SetPtEtaPhiM(b.photon_pt()->at(b.llphoton_iph()->at(0)), b.photon_eta()->at(b.llphoton_ill()->at(0)), b.photon_phi()->at(b.llphoton_ill()->at(0)), 0.0);
    llg = l1 + l2 + g;
    return llg.M();
  });

  NamedFunc ll_corr_mass("ll_corr_mass",[](const Baby &b) -> NamedFunc::ScalarType{
    TLorentzVector l1, l2, fsrph, ll;
    l1.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i1()->at(b.llphoton_ill()->at(0))), 0.105);
    l2.SetPtEtaPhiM(b.mu_pt()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_eta()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), b.mu_phi()->at(b.ll_i2()->at(b.llphoton_ill()->at(0))), 0.105);
    fsrph.SetPtEtaPhiM(b.fsrphoton_pt()->at(0), b.fsrphoton_eta()->at(0), b.fsrphoton_phi()->at(0), 0.0);
    if (b.fsrphoton_muonidx()->at(0)==b.ll_i1()->at(b.llphoton_ill()->at(0)) && l1.DeltaR(fsrph) < 0.4)
      l1 = l1 + fsrph;
    else if (b.fsrphoton_muonidx()->at(0)==b.ll_i2()->at(b.llphoton_ill()->at(0)) && l2.DeltaR(fsrph) < 0.4)
      l2 = l2 + fsrph;
    ll = l1 + l2;
    return ll.M();
  });

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight*100;
  });

  NamedFunc trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc baseline("nphoton > 0 && nll > 0 && nfsrphoton > 0");
  NamedFunc lep("ll_lepid[llphoton_ill[0]] == 13 && mu_pt[ll_i1[llphoton_ill[0]]] > 20 && mu_pt[ll_i2[llphoton_ill[0]]] > 10");
  NamedFunc pho("photon_pt[llphoton_iph[0]] > 15 && photon_sig[llphoton_iph[0]]");
  NamedFunc mass_cuts("ll_m[llphoton_ill[0]] > 50 && llphoton_m[0]+ll_m[llphoton_ill[0]] >= 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110");

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

  auto proc_UL = Process::MakeShared<Baby_pico>("HToZ#gamma", sig, TColor::GetColor("#ff0000"), {sig_path+"*GluGlu*.root"}, trigs);
  proc_UL->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_UL};

  PlotMaker pm;
  NamedFunc selection = baseline && lep && pho && mass_cuts;
  pm.Push<Hist1D>(Axis(20, 100, 180, llg_mass, "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("llphoton_m");
  pm.Push<Hist1D>(Axis(20, 50, 150, "ll_m[0]", "m_{ll} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ll_m");
  pm.Push<Hist1D>(Axis(20, 100, 180, llg_corr_mass, "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("llphoton_m_corr");
  pm.Push<Hist1D>(Axis(20, 50, 150, ll_corr_mass, "m_{ll} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("ll_m_corr");

  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

