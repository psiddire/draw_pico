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

int main(){
  gErrorIgnoreLevel = 6000;
  // string bfolder("/net/cms17/cms17r0/pico/NanoAODv7/");
  // string foldersig(bfolder+"zgamma_signal/2017/signal/unskimmed/*");
  string bfolder("/net/cms29/cms29r0/pico/NanoAODv2/");
  string foldersig(bfolder+"zgamma_signal_ul/2017/signal/unskimmed/*");

  NamedFunc sig_lepid("signal lepton ID",[](const Baby &b) -> NamedFunc::ScalarType{
    int lepid(0);
    for(size_t imc(0); imc < b.mc_id()->size(); imc++) {
      if(b.mc_mom()->at(imc) == 23 && b.mc_mom()->at(b.mc_momidx()->at(imc)) == 25) {
        lepid = abs(b.mc_id()->at(imc));
      }
    }
  return lepid;
  });
  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.w_lumi();
  });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc baseline("nphoton > 0 && nll > 0");
  NamedFunc lepton_cut("(ll_lepid[llphoton_ill[0]] == 13 && mu_pt[ll_i1[llphoton_ill[0]]] > 20 && mu_pt[ll_i2[llphoton_ill[0]]] > 10) || (ll_lepid[llphoton_ill[0]] == 11 && el_pt[ll_i1[llphoton_ill[0]]] > 25 && el_pt[ll_i2[llphoton_ill[0]]] > 15)");
  NamedFunc dilep_mass("ll_m[llphoton_ill[0]] > 50");
  NamedFunc photon_cut("photon_pt[llphoton_iph[0]] > 15 && photon_sig[llphoton_iph[0]]");
  NamedFunc llphoton_cuts("llphoton_m[0]+ll_m[llphoton_ill[0]] >= 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110");

  auto proc_el_sig = Process::MakeShared<Baby_pico>("Sig e",   Process::Type::signal, kMagenta-4, {foldersig+"*.root"}, sig_lepid == 11);
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("Sig #mu", Process::Type::signal, kAzure,     {foldersig+"*.root"}, sig_lepid == 13);
  vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig};

  vector<NamedFunc> cuts = {"1",
                            trigs,
			    baseline && lepton_cut,
                            dilep_mass,
			    photon_cut,
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
      TableRow("Base selections",                       cutflow.at(2),0,0,wgt),
      TableRow("lepton selections and $m_{ll} > $ 50 GeV", cutflow.at(3),0,0,wgt),
      TableRow("photon ID",                             cutflow.at(4),0,0,wgt),
      TableRow("three body invariant mass related cut", cutflow.at(5),0,0,wgt),
      },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}
