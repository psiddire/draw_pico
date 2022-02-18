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

#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

bool getDigit(int num, int place) { // Gives the 2^p place in binary number
  int temp = pow(2,place);
  return((num%(temp*2) - num%temp)/temp);
  }

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("/net/cms29/cms29r0/");
  string foldersig(bfolder+"/pico/NanoAODv5/zgamma_channelIslandsv3/2016/HToZG/unskimmed/*");
  string foldermc(bfolder+"/pico/NanoAODv5/zgamma_channelIslandsv2/2016/mc/unskimmed/*DYJetsToLL*");

  NamedFunc sig_lepid("signal lepton ID",[](const Baby &b) -> NamedFunc::ScalarType{
    int lepid(0);
    for(size_t imc(0); imc < b.mc_id()->size(); imc++) {
      if(b.mc_mom()->at(imc) == 23 && b.mc_mom()->at(b.mc_momidx()->at(imc)) == 25) {
        lepid = abs(b.mc_id()->at(imc));
//         if(lepid == 11 || lepid == 13 || lepid == 15)
          return lepid;
      }
    }
  return lepid;
  });

  NamedFunc dy_stitch("Drell-Yan stitch",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 reco, tru;
    for(size_t iph(0); iph < b.photon_pt()->size(); iph++) {
      reco.SetPtEtaPhi(b.photon_pt()->at(iph),b.photon_eta()->at(iph),b.photon_phi()->at(iph));
      for(size_t imc(0); imc < b.mc_id()->size(); imc++) {
        bitset<15> flags(b.mc_statusflag()->at(imc));
        if(b.mc_id()->at(imc) == 22 && b.mc_status()->at(imc) == 1) {
          if(flags[0] || flags[8]) {
            tru.SetPtEtaPhi(b.mc_pt()->at(imc),b.mc_eta()->at(imc),b.mc_phi()->at(imc));
            if(reco.DeltaR(tru) < 0.1) return false;
          }
        }
      }
    }
    return true;
  });
  NamedFunc lepton_sel("lepton_selection",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.ll_pt()->size() == 0) return false;
    for(size_t ill(0); ill < b.ll_pt()->size(); ill++) {
      if(b.ll_lepid()->at(ill) == 11) {
        for(size_t iel1(0); iel1 < b.el_pt()->size(); iel1++)
          for(size_t iel2(iel1+1); iel2 < b.el_pt()->size(); iel2++)
            if(b.el_pt()->at(iel1)*b.el_sig()->at(iel1) > 25 &&
               b.el_pt()->at(iel2)*b.el_sig()->at(iel2) > 15 &&
               b.el_charge()->at(iel1)+b.el_charge()->at(iel2) == 0)
               return true;
      }
      else if(b.ll_lepid()->at(ill) == 13) {
        for(size_t imu1(0); imu1 < b.mu_pt()->size(); imu1++)
          for(size_t imu2(imu1+1); imu2 < b.mu_pt()->size(); imu2++)
            if(b.mu_pt()->at(imu1)*b.mu_sig()->at(imu1) > 20 &&
               b.mu_pt()->at(imu2)*b.mu_sig()->at(imu2) > 10 &&
               b.mu_charge()->at(imu1)+b.mu_charge()->at(imu2) == 0)
               return true;
      }
    }
    return false;
  });
  NamedFunc nels("Number of electrons with WP98",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int nel(0);
    for(size_t iel(0); iel < b.el_pt()->size(); iel++) {
      double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                         { 0.604775,  0.628743,  0.896462 }};
      double pt(b.el_pt()->at(iel)), eta(abs(b.el_eta()->at(iel)));
      int ipt(1), ieta(-1);
      if(pt>10) ipt = 0;
      if(eta < 0.8) ieta = 0;
      else if(eta < 1.479) ieta = 1;
      else ieta = 2;
      if(b.el_idmva()->at(iel) > wp[ipt][ieta]) nel++;
    }
    return nel;
  });
  NamedFunc dilep_mass("m_{ll}",[](const Baby &b) -> NamedFunc::ScalarType{
    double mass(-1);
    for(int ill(0); ill < b.nll(); ill++) 
      if(mass == -1 || abs(b.ll_m()->at(ill) - 91.1876) < abs(mass - 91.1876))
        mass = b.ll_m()->at(ill);
    return mass;
  });
  NamedFunc photon_cut("photon_cut",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.photon_pt()->size() == 0) return false;
    for(size_t iph(0); iph < b.photon_pt()->size(); iph++) 
      if(b.photon_pt()->at(iph) > 15 && b.photon_sig()->at(iph)) return true;
    return false;
  });
//   NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
//     return abs(b.w_lumi());
//   });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ");
  auto proc_el_sig = Process::MakeShared<Baby_pico>("Sig e",  Process::Type::signal,     kMagenta-4, {foldersig+"*GluGluH*.root"}, sig_lepid == 11);
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("Sig #mu",Process::Type::signal,     kAzure,     {foldersig+"*GluGluH*.root"}, sig_lepid == 13);
  auto proc_el_dy  = Process::MakeShared<Baby_pico>("DY e",   Process::Type::background, kMagenta-4, {foldermc +"*.root"}, "ntruel>=2" && dy_stitch);
  auto proc_mu_dy  = Process::MakeShared<Baby_pico>("DY #mu", Process::Type::background, kAzure,     {foldermc +"*.root"}, "ntrumu>=2" && dy_stitch);
												
  NamedFunc llphoton_cuts("three body invariant mass cuts",[](const Baby &b) -> NamedFunc::ScalarType{
    bool llp_cuts(false);
    double mllp(0), mll(0), gpt(0);
    int best_ill(0);
    double min_dm(999);
    for(size_t ill(0); ill < b.ll_pt()->size(); ill++)
      if(fabs(b.ll_m()->at(ill) - 91.2) < min_dm)
        best_ill = ill;
    for(size_t illp(0); illp < b.llphoton_pt()->size(); illp++) {
      mllp = b.llphoton_m()->at(illp);
      mll  = b.ll_m()->at(b.llphoton_ill()->at(illp));
      gpt  = b.photon_pt()->at(b.llphoton_iph()->at(illp));
      if(b.llphoton_ill()->at(illp) == best_ill)
        if(mllp > 100 && mllp < 180)
          if(mllp + mll >= 185)
            if(gpt/mllp >= 15./110)
              llp_cuts = true;
    }
    return llp_cuts;
  });

	vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig};
// 	vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig, proc_el_dy, proc_mu_dy};

	vector<NamedFunc> cuts = {"1",
                            el_trigs || mu_trigs,
                            "nel>1 || nmu>1",
                            dilep_mass > 50,
                            "nphoton > 0",
                            llphoton_cuts};
	vector<NamedFunc> cutflow;
//   NamedFunc wgt("w_lumi*w_lep*w_photon");
  NamedFunc wgt("w_lumi");
	NamedFunc cut("1");
	for(size_t i = 0; i < cuts.size(); i++) {
	  cut = cut && cuts.at(i);
		cutflow.push_back(cut);
	}
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::overflow)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> ops = {log_lumi, lin_lumi};
  PlotMaker pm;
//   pm.Push<Hist1D>(Axis(2,0,2, !(sig_lepid == 11 || sig_lepid == 13), "ID != lep", {}),"1", samples,ops);
  pm.Push<Hist1D>(Axis(30,-0.5,29.5, sig_lepid, "signal lepton ID", {}),"1", samples,ops).Weight(wgt);
  pm.Push<Table>("AN_sig_cutflow", vector<TableRow>{
      TableRow("Total number of events",                cutflow.at(0),0,0,wgt),
	    TableRow("High level trigger",                    cutflow.at(1),0,0,wgt),
	    TableRow("lepton selections",                     cutflow.at(2),0,0,wgt),
	    TableRow("$m_{ll} > $ 50 GeV",                    cutflow.at(3),0,0,wgt),
	    TableRow("photon ID",                             cutflow.at(4),0,0,wgt),
	    TableRow("three body invariant mass related cut", cutflow.at(5),0,0,wgt),
// 	    TableRow("After including lep+photon SFs"       , cutflow.at(5),0,0,"w_lumi*w_lep*w_photon")
	    },samples,false);
  pm.Push<Table>("eff_sig_cutflow", vector<TableRow>{
      TableRow("Total number of events",                cutflow.at(0),0,0,wgt),
	    TableRow("High level trigger",                    cutflow.at(1),0,0,wgt),
	    TableRow("lepton selections",                     cutflow.at(2),0,0,wgt),
	    TableRow("$m_{ll} > $ 50 GeV",                    cutflow.at(3),0,0,wgt),
	    TableRow("photon ID",                             cutflow.at(4),0,0,wgt),
	    TableRow("three body invariant mass related cut", cutflow.at(5),0,0,wgt)
	    },samples,false, true, false, false, true);
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}
