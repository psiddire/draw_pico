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
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/");
  string foldersig(bfolder+"zgamma_signal_ul/2017/signal/unskimmed/*");

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

  NamedFunc lepton_sel("lepton_selection",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.ll_pt()->size() == 0) return false;
    for(size_t ill(0); ill < b.ll_pt()->size(); ill++) {
      if(b.ll_lepid()->at(ill) == 11) {
        for(size_t iel1(0); iel1 < b.el_pt()->size(); iel1++)
          for(size_t iel2(iel1+1); iel2 < b.el_pt()->size(); iel2++)
            if(b.el_pt()->at(iel1)*b.el_sig()->at(iel1) > 25 &&
               b.el_pt()->at(iel2)*b.el_sig()->at(iel2) > 15 &&
               b.el_charge()->at(iel1)*b.el_charge()->at(iel2) == -1)
               return true;
      }
      else if(b.ll_lepid()->at(ill) == 13) {
        for(size_t imu1(0); imu1 < b.mu_pt()->size(); imu1++)
          for(size_t imu2(imu1+1); imu2 < b.mu_pt()->size(); imu2++)
            if(b.mu_pt()->at(imu1)*b.mu_sig()->at(imu1) > 20 &&
               b.mu_pt()->at(imu2)*b.mu_sig()->at(imu2) > 10 &&
               b.mu_charge()->at(imu1)*b.mu_charge()->at(imu2) == -1)
               return true;
      }
    }
    return false;
  });
  NamedFunc nels("Number of electrons passing the ID",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int nel(0);
    for(size_t iel(0); iel < b.el_pt()->size(); iel++) {
      if(b.el_id()->at(iel)) nel++;
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
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.w_lumi();
  });
  
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  auto proc_el_sig = Process::MakeShared<Baby_pico>("Sig e",  Process::Type::signal,     kMagenta-4, {foldersig+"*.root"}, sig_lepid == 11);
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("Sig #mu",Process::Type::signal,     kAzure,     {foldersig+"*.root"}, sig_lepid == 13);
												
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

  vector<NamedFunc> cuts = {"1",
                            el_trigs || mu_trigs,
                            "nel>1 || nmu>1",
                            dilep_mass > 50,
                            "nphoton > 0",
                            llphoton_cuts};
  vector<NamedFunc> cutflow;

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
  pm.Push<Table>("AN_sig_cutflow", vector<TableRow>{
      TableRow("Total number of events",                cutflow.at(0),0,0,wgt),
	    TableRow("High level trigger",                    cutflow.at(1),0,0,wgt),
	    TableRow("lepton selections",                     cutflow.at(2),0,0,wgt),
	    TableRow("$m_{ll} > $ 50 GeV",                    cutflow.at(3),0,0,wgt),
	    TableRow("photon ID",                             cutflow.at(4),0,0,wgt),
	    TableRow("three body invariant mass related cut", cutflow.at(5),0,0,wgt),
	    },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}
