#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TObject.h"
#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
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
  string folder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_data/2017/");
  string folderdata( folder+"data/merged_zgmc_llg/");

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  auto proc_el_sig = Process::MakeShared<Baby_pico>("DoubleEG",  Process::Type::data, kBlack, {folderdata+"*DoubleEG*.root"}, el_trigs);
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("DoubleMuon",Process::Type::data, kBlack, {folderdata+"*DoubleMuon*.root"}, mu_trigs);

  vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig};

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
      return b.w_lumi();
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

  PlotMaker pm;
  pm.Push<Table>("ATLAS", vector<TableRow>{
      TableRow("VBF enriched",             baseline && dijet >= 400 && (lep.at(0) || lep.at(1)), 0, 0, wgt),
	TableRow("High relative $p_T$",      baseline && relpT > 0.4 && dijet < 400 && (lep.at(0) || lep.at(1)), 0, 0, wgt),
	TableRow("High $p_{Tt}$ ee",         baseline && relpT < 0.40 && pTt > 40 && dijet < 400 && lep.at(0), 0, 0, wgt),
	TableRow("Low $p_{Tt}$ ee",          baseline && relpT < 0.40 && pTt < 40 && dijet < 400 && lep.at(0), 0, 0, wgt),
	TableRow("High $p_{Tt}$ $\\mu\\mu$", baseline && relpT < 0.40 && pTt > 40 && dijet < 400 && lep.at(1), 0, 0, wgt),
	TableRow("Low $p_{Tt}$ $\\mu\\mu$",  baseline && relpT < 0.40 && pTt < 40 && dijet < 400 && lep.at(1), 0, 0, wgt),
	TableRow("Inclusive",                baseline && (lep.at(0) || lep.at(1)), 1, 0, wgt)
  },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}


