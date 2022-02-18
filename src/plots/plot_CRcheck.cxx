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
  Palette col("txt/colors.txt","default");
  Process::Type back =  Process::Type::background;
  Process::Type sig  =  Process::Type::signal;
  string bfolder("/net/cms29/cms29r0/pico/NanoAODv5/zgamma_channelIslandsv3/2016/");
  string mc_path( bfolder+"mc/merged_zgmc_llg/");
  string sig_path(bfolder+"HToZG/merged_zgmc_llg/");
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_smzg  = Process::MakeShared<Baby_pico>("SM Z#gamma" , back, col("smzg") , {mc_path+"*ZGToLLG*"}, trigs);
  auto proc_ewkzg = Process::MakeShared<Baby_pico>("EWK Z#gamma", back, col("ewkzg"), {mc_path+"*LLAJJ*"},   trigs);
  auto proc_dy    = Process::MakeShared<Baby_pico>("DY     ",     back, col("dy")   , {mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",       back, col("ttg")  , {mc_path+"*TTGJets*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZg(x100)",  sig, col("hzg")  , {sig_path+"*.root"},   trigs);
  proc_smzg->SetLineWidth (0); proc_ttg->SetLineWidth  (0); proc_dy->SetLineWidth   (0); proc_ewkzg->SetLineWidth(0);
  proc_hzg->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_dy, proc_smzg, proc_ttg, proc_ewkzg, proc_hzg};
  NamedFunc els("e-channel",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 25 && b.el_pt()->at(b.ll_i2()->at(0)) > 15;
  });
  NamedFunc mus("#mu-channel",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 20 && b.mu_pt()->at(b.ll_i2()->at(0)) > 10;
  });
  NamedFunc baseline("base",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.llphoton_m()->at(0) + b.ll_m()->at(0) > 185 && b.llphoton_m()->at(0) > 100 && b.llphoton_m()->at(0) < 180 &&
           b.photon_pt()->at(0)/b.llphoton_m()->at(0) >= 15./110 && b.photon_drmin()->at(0) > 0.4;
  });
  NamedFunc baseline_flip("base_mllmllgFlipped",[](const Baby &b) -> NamedFunc::ScalarType{
    return //b.llphoton_m()->at(0) > 100 && 
           b.llphoton_m()->at(0) < 180 &&
           b.ll_m()->at(0) + b.llphoton_m()->at(0) <= 185 &&
           b.photon_pt()->at(0)/b.llphoton_m()->at(0) >= 15./110 && 
           b.photon_drmin()->at(0) > 0.4;
  });
  NamedFunc mll_mllg("mll_mllg", [](const Baby &b) -> NamedFunc::ScalarType{ return b.ll_m()->at(0)+b.llphoton_m()->at(0); });
  NamedFunc p_pte( "photon_pte", [](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pterr()->at(0)/b.photon_pt()->at(0); });
  NamedFunc photon_drj("photon_drj",[](const Baby &b) -> NamedFunc::ScalarType{
    double drjmin(999);
    if(b.njet() == 0) return drjmin;
    TVector3 photon = AssignGamma(b).Vect();
    for(int ij = 0; ij < b.njet(); ij++) {
      TVector3 jet;
      jet.SetPtEtaPhi(b.jet_pt()->at(ij),b.jet_eta()->at(ij),b.jet_phi()->at(ij));
      if(photon.DeltaR(jet) < drjmin) drjmin = photon.DeltaR(jet);
    }
    return drjmin;
  });
  NamedFunc cr1("cr1",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_m()->at(0) > 50 && b.llphoton_m()->at(0) > 100 && b.llphoton_m()->at(0) < 180 &&
           b.photon_pt()->at(b.llphoton_iph()->at(0)) > 15 && 
           b.photon_drmin()->at(0) > 0.4 && 
           (b.llphoton_m()->at(0) + b.ll_m()->at(0) < 185 || b.photon_pt()->at(0)/b.llphoton_m()->at(0) < 15./110);
  });
  vector<NamedFunc> lep = {els, mus};
  vector<string>  tlep =  {"el","mu"};
  PlotOpt lin_stack("txt/plot_styles.txt","CMSPaper");
  lin_stack.Title(TitleType::info).Stack(StackType::signal_overlay)
           .UseCMYK(false) .FileExtensions({"pdf"});
  vector<PlotOpt> ops = {lin_stack};
  vector<PlotOpt> logops = {lin_stack().YAxis(YAxisType::log)};
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.w_lumi();
    if(b.type() >= 200000 && b.type() <= 205000)
      return 100*weight;
    else if(b.type() == 6200)
      return weight/1.664;
    return        weight;
  });
  PlotMaker pm;
  NamedFunc dijet("njet>=2");
  for(int i(0); i < 2; i++) {
    NamedFunc selection = lep.at(i) && baseline && dijet;
    string tag = "base_" + tlep.at(i);
//    pm.Push<Hist1D>(Axis(80,100,180,  "llphoton_m[0]",    "m_{ll#gamma} [GeV]" ,{}),  cr1 && lep.at(i), procs, ops).Weight(wgt).Tag("DYCR");
   pm.Push<Hist1D>(Axis(20,0,1,  "photon_reliso[0]",    "Photon reliso" ,{}), baseline && lep.at(i), procs, logops).Weight(wgt).Tag("DYCR");
   pm.Push<Hist1D>(Axis(20,0,1,  "photon_r9[0]",    "Photon R9" ,{}), baseline && lep.at(i), procs, logops).Weight(wgt).Tag("DYCR");
   pm.Push<Hist1D>(Axis(20,0,0.5,  "photon_hoe[0]",    "Photon H/E" ,{}), baseline && lep.at(i) , procs, ops).Weight(wgt).Tag("DYCR");
   pm.Push<Hist1D>(Axis(16,0.2,1,  "photon_idmva[0]",    "Photon MVA" ,{}), baseline && lep.at(i) , procs, logops).Weight(wgt).Tag("DYCR");
//    pm.Push<Hist1D>(Axis(20,0,4,  photon_drj,    "#Delta R_{min}(#gamma,j)" ,{}), baseline && "njet>0" && lep.at(i), procs, ops).Weight(wgt).Tag("DYCR");
//    pm.Push<Hist1D>(Axis(25,0,1000,  "dijet_m",    "Dijet mass [GeV]" ,{}),  selection, procs, ops).Weight(wgt).Tag("Dijet");
//    pm.Push<Hist1D>(Axis(80,60,140,  "llphoton_m[0]",    "m_{ll#gamma} [GeV]" ,{}),  baseline_flip && lep.at(i), procs, ops).Weight(wgt).Tag("CR_check");
//    pm.Push<Hist1D>(Axis(80,100,180,  "llphoton_m[0]",    "m_{ll#gamma} [GeV]" ,{}),  baseline && lep.at(i), procs, ops).Weight(wgt).Tag("CR_check");
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}
