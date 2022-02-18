#include "TChain.h"
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

  TChain *tHZG = new TChain("tree");
  tHZG->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/*.root");

  TChain *tHMM = new TChain("tree");
  tHMM->Add("/net/cms17/cms17r0/pico/NanoAODv7/zgamma_mm/2017/mm/merged_zgmc_llg/*.root");

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

  string cut = mu_trigs+"&&"+baseline+"&&"+mu_sel+"&&"+pho+"&&"+mass_cuts;

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
  pad1->cd();
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(10);
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetLeftMargin(0.10);
  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.10);
  pad1->SetBottomMargin(0.10);
  pad1->SetFrameFillStyle(0);
  pad1->SetFrameLineStyle(0);
  pad1->SetFrameLineWidth(3);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderSize(10);

  string prefix_HZG = "llg_HZG";
  string prefix_HMM = "llg_HMM";
  string x_axis = "photon_drmin[llphoton_iph[0]]"; // llphoton_m[0] // ll_m[llphoton_ill[0]] // photon_pt[llphoton_iph[0]] // mu_pt[ll_i2[llphoton_ill[0]]] // photon_drmin[llphoton_iph[0]]

  TH1D* aHZG = new TH1D(prefix_HZG.c_str(), title.c_str(), 100, 0.0, 4.0); // 80, 100, 180 // 100, 50, 150 // 90, 10, 100 // 120, 10, 130 // 100, 0.0, 4.0
  tHZG->Draw((x_axis+">>"+prefix_HZG).c_str(), ("("+cut+")*("+weight+")*("+lumi+")").c_str(), "goff");
  aHZG->Scale(1./aHZG->Integral());
  aHZG->SetLineColor(kRed);
  aHZG->SetLineWidth(3);
  // aHZG->GetXaxis()->SetTitle("min #Delta R(l, #gamma)"); // m_{ll#gamma} [GeV] // m_{ll} [GeV] // p_{T}^{#gamma} [GeV] // p_{T}^{#mu, 2} [GeV] // min #Delta R(l, #gamma)
  aHZG->GetXaxis()->SetLabelSize(0.04);
  aHZG->GetYaxis()->SetLabelSize(0.04);
  cout<<"number analysis: "<<aHZG->Integral()<<endl;

  TH1D* aHMM = new TH1D(prefix_HMM.c_str(), title.c_str(), 100, 0.0, 4.0); // 80, 100, 180 // 100, 50, 150 // 90, 10, 100 // 120, 10, 130 // 100, 0.0, 4.0
  tHMM->Draw((x_axis+">>"+prefix_HMM).c_str(), ("("+cut+")*("+weight+")*("+lumi+")").c_str(), "goff");
  aHMM->Scale(1./aHMM->Integral());
  aHMM->GetXaxis()->SetTitle("min #Delta R(l, #gamma)"); // m_{ll#gamma} [GeV] // m_{ll} [GeV] // p_{T}^{#gamma} [GeV] // p_{T}^{#mu, 2} [GeV] // min #Delta R(l, #gamma)
  aHMM->SetLineColor(kBlue);
  aHMM->SetLineWidth(3);
  cout<<"number analysis: "<<aHMM->Integral()<<endl;

  TLegend *l = new TLegend(0.5, 0.65, 0.9, 0.85);
  l->AddEntry(aHZG, "H #rightarrow Z#gamma", "l");
  l->AddEntry(aHMM, "H #rightarrow #mu#mu", "l");
  
  aHMM->Draw("HIST");
  aHZG->Draw("HISTSAME");
  l->Draw();
  c1->cd();
  pad1->Draw();
  c1->Update();
  c1->SaveAs("dr.png"); // llg.png // ll.png // ptm2.png // dr.png
}

