#include "TChain.h"
#include "TList.h"
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

  TChain * tree = new TChain("tree");
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_tt/2017/tt/unskimmed/*.root");

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

  TFile *f = new TFile("HTT.root", "RECREATE");

  string cut_mu = mu_trigs+"&&"+baseline+"&&"+mu_sel+"&&"+pho+"&&"+mass_cuts;
  string cut_el = el_trigs+"&&"+baseline+"&&"+el_sel+"&&"+pho+"&&"+mass_cuts;

  string x_axis = "llphoton_m[0]";

  string prefix_mu = "llg_mu";
  TH1D* analysis_mu = new TH1D(prefix_mu.c_str(), title.c_str(), 80, 100, 180);
  tree->Draw((x_axis+">>"+prefix_mu).c_str(), ("("+cut_mu+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<analysis_mu->Integral()<<endl;

  string prefix_el = "llg_el";
  TH1D* analysis_el = new TH1D(prefix_el.c_str(), title.c_str(), 80, 100, 180);
  tree->Draw((x_axis+">>"+prefix_el).c_str(), ("("+cut_el+")*("+weight+")*("+lumi+")").c_str(), "goff");
  cout<<"number analysis: "<<analysis_el->Integral()<<endl;

  analysis_mu->Write();
  analysis_el->Write();
  f->ls();
}

