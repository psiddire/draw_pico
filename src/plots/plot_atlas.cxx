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
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
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
  Process::Type data = Process::Type::data;
  Process::Type sig = Process::Type::signal;

  string folder("/net/cms17/cms17r0/pico/NanoAODv2/");
  string folderdata( folder+"zgamma_data/2017/data/merged_zgmc_llg/");
  string foldersig( folder+"zgamma_signal_ul/2017/signal/merged_zgmc_llg/");

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  auto proc_data = Process::MakeShared<Baby_pico>("Data", data, kBlack, {folderdata+"*.root"}, trigs);
  auto proc_hzg  = Process::MakeShared<Baby_pico>("HToZ#gamma(x20)", sig, kRed, {foldersig+"*.root"}, trigs);

  proc_data->SetMarkerSize(1);
  proc_hzg->SetLineWidth(3);
  vector<shared_ptr<Process> > procs = {proc_data, proc_hzg};

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

  NamedFunc dijet("dijet",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector j1, j2;
      j1.SetPtEtaPhiM(0, 0, 0, 0);
      j2.SetPtEtaPhiM(0, 0, 0, 0);
      int i1(-1), i2(-1);
      for(size_t ij(0); ij < b.jet_pt()->size(); ij++) {
        if(b.jet_isgood()->at(ij) && i1 == -1) {
          i1 = ij;
          j1.SetPtEtaPhiM(b.jet_pt()->at(i1), b.jet_eta()->at(i1), b.jet_phi()->at(i1), b.jet_m()->at(i1));
        } else if(b.jet_isgood()->at(ij) && i2 == -1) {
          i2 = ij;
          j2.SetPtEtaPhiM(b.jet_pt()->at(i2), b.jet_eta()->at(i2), b.jet_phi()->at(i2), b.jet_m()->at(i2));
        }
      }
      return (j1 + j2).M();
    });
  NamedFunc pTt("pTt",[](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 g = AssignGamma(b, false).Vect();
      TVector3 h = AssignH(b, false).Vect();
      TVector3 z = AssignZ(b, false).Vect();
      g.SetZ(0); h.SetZ(0); z.SetZ(0);
      return h.Cross((z-g).Unit()).Mag();
    });
  NamedFunc gen_relpT("gen_relpT",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector h = AssignH(b, true);
      TLorentzVector g = AssignGamma(b, true);
      return g.Pt()/h.M();
    });
  NamedFunc relpT("relpT",[](const Baby &b) -> NamedFunc::ScalarType{
      double iph = b.llphoton_iph()->at(0);
      double pt = b.photon_pt()->at(iph);
      double mass = b.llphoton_m()->at(0);
      return pt/mass;
    });
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
      double weight = b.w_lumi();
      if(b.type() >= 200000 && b.type() <= 205000)
	return 20*weight;
      return weight;
    });

  NamedFunc baseline("nphoton > 0 && nll > 0");
  NamedFunc pho("photon_pt[llphoton_iph[0]] > 15 && photon_drmin[llphoton_iph[0]] > 0.4 && photon_sig[llphoton_iph[0]]");
  vector<NamedFunc> bline = {"ll_m[llphoton_ill[0]] > 50.0",
                             pho,
                             relpT > 15./110.,
                             "llphoton_m[0] > 100 && llphoton_m[0] < 180",
                             "llphoton_m[0] + ll_m[llphoton_ill[0]] >= 185"};
  for(size_t i = 0; i < bline.size(); i++)
    baseline = baseline && bline.at(i);
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

  vector<NamedFunc> cat = { baseline && dijet >= 400 && (lep.at(0) || lep.at(1)),
			    baseline && relpT > 0.4 && dijet < 400 && (lep.at(0) || lep.at(1)),
			    baseline && relpT < 0.40 && pTt > 40 && dijet < 400 && lep.at(0),
			    baseline && relpT < 0.40 && pTt < 40 && dijet < 400 && lep.at(0),
			    baseline && relpT < 0.40 && pTt > 40 && dijet < 400 && lep.at(1),
			    baseline && relpT < 0.40 && pTt < 40 && dijet < 400 && lep.at(1)};

  PlotMaker pm;
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(0), procs, ops).Weight(wgt).Tag("cat1");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(1), procs, ops).Weight(wgt).Tag("cat2");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(2), procs, ops).Weight(wgt).Tag("cat3");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(3), procs, ops).Weight(wgt).Tag("cat4");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(4), procs, ops).Weight(wgt).Tag("cat5");
  pm.Push<Hist1D>(Axis(40, 100, 180, "llphoton_m[0]", "m_{Z#gamma} [GeV]", {}), cat.at(5), procs, ops).Weight(wgt).Tag("cat6");
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

