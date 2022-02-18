#include "TTree.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "THStack.h"
#include "TStyle.h"
#include <iostream>
#include "TLorentzVector.h"

using std::string;
using std::to_string;
using std::set;
using std::cout;
using std::endl;


int main() {

  TFile *fGG = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_GluGluHToZG_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_5.root");
  TFile *fVBF = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_VBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_3.root");
  TFile *fWM = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_WminusH_HToZG_WToAll_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_4.root");
  TFile *fWP = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_WplusH_HToZG_WToAll_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_3.root");
  TFile *fZH = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_ZH_HToZG_ZToAll_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_6.root");
  TFile *fttH = new TFile("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/merged_pico_llg_ttHToZG_M125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_4.root");

  TTree *tGG = (TTree*)fGG->Get("tree");
  TTree *tVBF = (TTree*)fVBF->Get("tree");
  TTree *tWM = (TTree*)fWM->Get("tree");
  TTree *tWP = (TTree*)fWP->Get("tree");
  TTree *tZH = (TTree*)fZH->Get("tree");
  TTree *tttH = (TTree*)fttH->Get("tree");

  string weight = "w_lumi";

  string lumi = "41.5";

  string title = "";

  string el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  string mu_trigs = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)";
  string trigs = el_trigs+"||"+mu_trigs;
  string mass_cuts = "(ll_m[llphoton_ill[0]] > 50 && "
    "llphoton_m[0]+ll_m[llphoton_ill[0]]>=185 && "
    "llphoton_m[0] > 100 && llphoton_m[0] < 180 && "
    "photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110)";
  string baseline = "(nphoton > 0 && nll > 0)";
  string el_sel = "(ll_lepid[llphoton_ill[0]] == 11 && "
    "el_pt[ll_i1[llphoton_ill[0]]] > 25 && "
    "el_pt[ll_i2[llphoton_ill[0]]] > 15 && "
    "el_sig[ll_i1[llphoton_ill[0]]] && "
    "el_sig[ll_i2[llphoton_ill[0]]])";
  string mu_sel = "(ll_lepid[llphoton_ill[0]] == 13 && "
    "mu_pt[ll_i1[llphoton_ill[0]]] > 20 && "
    "mu_pt[ll_i2[llphoton_ill[0]]] > 10 && "
    "mu_sig[ll_i1[llphoton_ill[0]]] && "
    "mu_sig[ll_i2[llphoton_ill[0]]])";
  string pho = "(photon_pt[llphoton_iph[0]] > 15 && "
    "photon_drmin[llphoton_iph[0]] > 0.4 && "
    "photon_sig[llphoton_iph[0]])";

  string cut_mu = mu_trigs+"&&"+baseline+"&&"+mu_sel+"&&"+pho+"&&"+mass_cuts;
  string cut_el = el_trigs+"&&"+baseline+"&&"+el_sel+"&&"+pho+"&&"+mass_cuts;

  string prefix_GGmu = "llg_GGmu";
  string prefix_GGel = "llg_GGel";
  string prefix_VBFmu = "llg_VBFmu";
  string prefix_VBFel = "llg_VBFel";
  string prefix_ZHmu = "llg_ZHmu";
  string prefix_ZHel = "llg_ZHel";
  string prefix_WMmu = "llg_WMmu";
  string prefix_WMel = "llg_WMel";
  string prefix_WPmu = "llg_WPmu";
  string prefix_WPel = "llg_WPel";
  string prefix_ttHmu = "llg_ttHmu";
  string prefix_ttHel = "llg_ttHel";
  string x_axis = "llphoton_m[0]";

  TFile *f = new TFile("Signal.root", "RECREATE");

  TH1D* GGmu = new TH1D(prefix_GGmu.c_str(), title.c_str(), 80, 100, 180);
  tGG->Draw((x_axis+">>"+prefix_GGmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<GGmu->Integral()<<endl;

  TH1D* GGel = new TH1D(prefix_GGel.c_str(), title.c_str(), 80, 100, 180);
  tGG->Draw((x_axis+">>"+prefix_GGel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<GGel->Integral()<<endl;

  TH1D* VBFmu = new TH1D(prefix_VBFmu.c_str(), title.c_str(), 80, 100, 180);
  tVBF->Draw((x_axis+">>"+prefix_VBFmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<VBFmu->Integral()<<endl;

  TH1D* VBFel = new TH1D(prefix_VBFel.c_str(), title.c_str(), 80, 100, 180);
  tVBF->Draw((x_axis+">>"+prefix_VBFel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<VBFel->Integral()<<endl;

  TH1D* WMmu = new TH1D(prefix_WMmu.c_str(), title.c_str(), 80, 100, 180);
  tWM->Draw((x_axis+">>"+prefix_WMmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<WMmu->Integral()<<endl;

  TH1D* WMel = new TH1D(prefix_WMel.c_str(), title.c_str(), 80, 100, 180);
  tWM->Draw((x_axis+">>"+prefix_WMel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<WMel->Integral()<<endl;

  TH1D* WPmu = new TH1D(prefix_WPmu.c_str(), title.c_str(), 80, 100, 180);
  tWP->Draw((x_axis+">>"+prefix_WPmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<WPmu->Integral()<<endl;

  TH1D* WPel = new TH1D(prefix_WPel.c_str(), title.c_str(), 80, 100, 180);
  tWP->Draw((x_axis+">>"+prefix_WPel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<WPel->Integral()<<endl;

  TH1D* ZHmu = new TH1D(prefix_ZHmu.c_str(), title.c_str(), 80, 100, 180);
  tZH->Draw((x_axis+">>"+prefix_ZHmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<ZHmu->Integral()<<endl;

  TH1D* ZHel = new TH1D(prefix_ZHel.c_str(), title.c_str(), 80, 100, 180);
  tZH->Draw((x_axis+">>"+prefix_ZHel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<ZHel->Integral()<<endl;

  TH1D* ttHmu = new TH1D(prefix_ttHmu.c_str(), title.c_str(), 80, 100, 180);
  tttH->Draw((x_axis+">>"+prefix_ttHmu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<ttHmu->Integral()<<endl;

  TH1D* ttHel = new TH1D(prefix_ttHel.c_str(), title.c_str(), 80, 100, 180);
  tttH->Draw((x_axis+">>"+prefix_ttHel).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<ttHel->Integral()<<endl;

  GGmu->Write();
  GGel->Write();
  VBFmu->Write();
  VBFel->Write();
  WMmu->Write();
  WMel->Write();
  WPmu->Write();
  WPel->Write();
  ZHmu->Write();
  ZHel->Write();
  ttHmu->Write();
  ttHel->Write();
  f->Close();

}

