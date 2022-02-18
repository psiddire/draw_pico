#include <iostream>
#include <string>
#include <vector>
#include <memory>
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
#include "core/hist2d.hpp"
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
  auto proc_ttg   = Process::MakeShared<Baby_pico>("ttbar",       back, col("ttg")  , {mc_path+"*TT_Tune*"}, trigs);
  auto proc_hzg   = Process::MakeShared<Baby_pico>("HToZg(x1000)",back, col("hzg")  , {sig_path+"*.root"},   trigs);
  auto proc_hzgs  = Process::MakeShared<Baby_pico>("HToZg"       ,sig,  kMagenta  , {sig_path+"*.root"},   trigs);
  vector<shared_ptr<Process>> back_procs = {proc_dy, proc_smzg, proc_ttg, proc_ewkzg};
  vector<shared_ptr<Process>> sig_procs = {proc_hzg};
  NamedFunc baseline("baseline",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_m()->at(0) > 50 && b.llphoton_m()->at(0) > 80 && b.llphoton_m()->at(0) < 100 &&
           b.photon_pt()->at(b.llphoton_iph()->at(0)) > 15 && 
           b.photon_drmin()->at(0) > 0.4;
//            b.photon_pt()->at(0)/b.llphoton_m()->at(0) > 15./110;
  });
  NamedFunc els("e-channel",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_lepid()->at(0) == 11 && b.el_pt()->at(b.ll_i1()->at(0)) > 25 && b.el_pt()->at(b.ll_i2()->at(0)) > 15;
  });
  NamedFunc mus("#mu-channel",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.ll_lepid()->at(0) == 13 && b.mu_pt()->at(b.ll_i1()->at(0)) > 20 && b.mu_pt()->at(b.ll_i2()->at(0)) > 10;
  });
  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{ 
    double weight = b.weight();
    if(b.type() >= 200000 && b.type() <= 205000)
//       return 1000000*weight;
      return 1*weight;
    return        weight;
  });
  PlotOpt style("txt/plot_styles.txt","Scatter");
  vector<PlotOpt> bkg_hist = {style().Stack(StackType::data_norm)
                                     .LogMinimum(1)
                                     .CanvasWidth(600)
                                     .LabelSize(0.04)
                                     .YAxis(YAxisType::log)
                                     .Title(TitleType::info)};
  PlotMaker pm;
  vector<vector<shared_ptr<Process>>> procs = {back_procs, sig_procs, {proc_dy}, {proc_smzg}};
  vector<string> tags = {"2D_background","2D_signal",
                         "2D_back_DY", "2D_back_SMZG"};
  for(int i = 0; i < 1; i++) {
//     pm.Push<Hist2D>(Axis(50,50,150, "ll_m[0]",       "m_{ee} [GeV]",{}),
//                     Axis(50,50,300, "llphoton_m[0]", "m_{ee#gamma} [GeV]",{100,80}), 
//                     baseline && els, procs.at(i), bkg_hist).Tag(tags.at(i)).Weight(wgt);
//     pm.Push<Hist2D>(Axis(50,50,150, "ll_m[0]",       "m_{#mu#mu} [GeV]",{}),
//                     Axis(50,50,300, "llphoton_m[0]", "m_{#mu#mu#gamma} [GeV]",{100,80}), 
//                     baseline && mus, procs.at(i), bkg_hist).Tag(tags.at(i)).Weight(wgt);
    pm.Push<Hist2D>(
                    Axis(50,120,320, "ll_m[0]+llphoton_m[0]",      "m_{#mu#mu}+m_{#mu#mu#gamma} [GeV]",{}),
                    Axis(50,0,1,     "photon_pt[0]/llphoton_m[0]", "p_{T,#gamma}/m_{#mu#mu#gamma}",{}), 
                    baseline && mus, procs.at(i), bkg_hist).Tag(tags.at(i)).Weight(wgt);
//     pm.Push<Hist2D>(Axis(50,0,100,  "photon_pt[0]", "photon p_{T} [GeV]",{}),
//                     Axis(50,50,300, "llphoton_m[0]", "m_{#mu#mu#gamma} [GeV]",{100,80}), 
//                     baseline && mus, procs.at(i), bkg_hist).Tag(tags.at(i)).Weight(wgt);
  }
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}
