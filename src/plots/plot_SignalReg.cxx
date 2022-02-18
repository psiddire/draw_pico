#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TStyle.h"
#include <iostream>
#include "TLorentzVector.h"

using std::string;
using std::to_string;
using std::set;
using std::cout;
using std::endl;


TCanvas * newCanvas(string const & name = "", int size = 500)
{
  TSeqCollection * canvases = gROOT->GetListOfCanvases();
  double iCanvas = canvases->GetEntries();
  string canvasName;
  if (name == "") canvasName = "c_g_" + to_string(iCanvas++);
  else canvasName = name;
  return new TCanvas(canvasName.c_str(), canvasName.c_str(), size, size);
}


int main() {

  TChain * tree = new TChain("tree");
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/*.root");

  string el_trigs = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
  string mu_trigs = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)";
  string baseline = "(nphoton > 0 && nll > 0)";
  string electron = "(ll_lepid[llphoton_ill[0]] == 11 && "
    "el_pt[ll_i1[llphoton_ill[0]]] > 25 && "
    "el_pt[ll_i2[llphoton_ill[0]]] > 15 && "
    "el_sig[ll_i1[llphoton_ill[0]]] && "
    "el_sig[ll_i2[llphoton_ill[0]]])";
  string muon = "(ll_lepid[llphoton_ill[0]] == 13 && "
    "mu_pt[ll_i1[llphoton_ill[0]]] > 20 && "
    "mu_pt[ll_i2[llphoton_ill[0]]] > 10 && "
    "mu_sig[ll_i1[llphoton_ill[0]]] && "
    "mu_sig[ll_i2[llphoton_ill[0]]])";
  string photon = "(photon_pt[llphoton_iph[0]] > 15 && "
    "photon_drmin[llphoton_iph[0]] > 0.4 && "
    "photon_sig[llphoton_iph[0]])";
  string mass_cuts = "(ll_m[llphoton_ill[0]] > 50 && "
    "llphoton_m[0]+ll_m[llphoton_ill[0]]>=185 && "
    "llphoton_m[0] > 100 && llphoton_m[0] < 180 && "
    "photon_pt[llphoton_iph[0]]/llphoton_m[0] >= 15./110)";
  string el_analysis = el_trigs+"&&"+baseline+"&&"+electron+"&&"+photon+"&&"+mass_cuts;
  string mu_analysis = mu_trigs+"&&"+baseline+"&&"+muon+"&&"+photon+"&&"+mass_cuts;
  string weight = "1"; // w_lumi

  string x_axis = "sqrt(llphoton_dml1[0]^2 + llphoton_dml2[0]^2 + llphoton_dmph[0]^2) / llphoton_m[0]";
  string x_name = "Mass resolution";

  int nbins = 100;
  float xmin = 0.0;
  float xmax = 0.05;
  float lumi = 41.5;
  string folder_name = "plots/Res/";

  gStyle->SetOptStat(0);
  TCanvas * t_c = 0;
  t_c = newCanvas();

  // string prefix_b = "llg_res_b";
  // string prefix_e = "llg_res_e";
  string title = "";
  // string mu1b = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 0.9)";
  // string mu1o = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) > 0.9 && abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 1.8)";
  // string mu1e = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) > 1.8 && abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 2.4)";
  // string mu2b = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 0.9)";
  // string mu2o = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) > 0.9 && abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 1.8)";
  // string mu2e = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) > 1.8 && abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 2.4)";

  string el1b = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) < 0.8)";
  string el1o = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) > 0.8 && abs(el_eta[ll_i1[llphoton_ill[0]]]) < 1.5)";
  string el1e = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) > 1.5 && abs(el_eta[ll_i1[llphoton_ill[0]]]) < 2.5)";
  string el2b = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) < 0.8)";
  string el2o = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) > 0.8 && abs(el_eta[ll_i2[llphoton_ill[0]]]) < 1.5)";
  string el2e = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) > 1.5 && abs(el_eta[ll_i2[llphoton_ill[0]]]) < 2.5)";

  // string mu_analysis_bb = mu_analysis+"&&"+mu1b+"&&"+mu2b;
  // string mu_analysis_bo = mu_analysis+"&&"+mu1b+"&&"+mu2o;
  // string mu_analysis_be = mu_analysis+"&&"+mu1b+"&&"+mu2e;
  // string mu_analysis_ob = mu_analysis+"&&"+mu1o+"&&"+mu2b;
  // string mu_analysis_oo = mu_analysis+"&&"+mu1o+"&&"+mu2o;
  // string mu_analysis_oe = mu_analysis+"&&"+mu1o+"&&"+mu2e;
  // string mu_analysis_eb = mu_analysis+"&&"+mu1e+"&&"+mu2b;
  // string mu_analysis_eo = mu_analysis+"&&"+mu1e+"&&"+mu2o;
  // string mu_analysis_ee = mu_analysis+"&&"+mu1e+"&&"+mu2e;

  string el_analysis_bb = el_analysis+"&&"+el1b+"&&"+el2b;
  string el_analysis_bo = el_analysis+"&&"+el1b+"&&"+el2o;
  string el_analysis_be = el_analysis+"&&"+el1b+"&&"+el2e;
  string el_analysis_ob = el_analysis+"&&"+el1o+"&&"+el2b;
  string el_analysis_oo = el_analysis+"&&"+el1o+"&&"+el2o;
  string el_analysis_oe = el_analysis+"&&"+el1o+"&&"+el2e;
  string el_analysis_eb = el_analysis+"&&"+el1e+"&&"+el2b;
  string el_analysis_eo = el_analysis+"&&"+el1e+"&&"+el2o;
  string el_analysis_ee = el_analysis+"&&"+el1e+"&&"+el2e;

  string prefix[9] = {"llg_res_bb", "llg_res_bo", "llg_res_be", "llg_res_ob", "llg_res_oo", "llg_res_oe", "llg_res_eb", "llg_res_eo", "llg_res_ee"};
  // string cut[9] = {mu_analysis_bb, mu_analysis_bo, mu_analysis_be, mu_analysis_ob, mu_analysis_oo, mu_analysis_oe, mu_analysis_eb, mu_analysis_eo, mu_analysis_ee};
  string cut[9] = {el_analysis_bb, el_analysis_bo, el_analysis_be, el_analysis_ob, el_analysis_oo, el_analysis_oe, el_analysis_eb, el_analysis_eo, el_analysis_ee};
  int color[9] = {1, 632, 600, 416, 880, 800, 432, 900, 616};
  string lName[9] = {"BB", "BO", "BE", "OB", "OO", "OE", "EB", "EO", "EE"};

  TLegend *legend = new TLegend(0.5, 0.5, 0.88, 0.9);

  for (int i = 0; i < 9; i++) {
    TH1D* analysis = new TH1D(prefix[i].c_str(), title.c_str(), nbins, xmin, xmax);
    analysis->GetXaxis()->SetTitle(x_name.c_str());
    analysis->GetYaxis()->SetRangeUser(0.0, 0.2);
    tree->Draw((x_axis+">>"+prefix[i]).c_str(), ("("+cut[i]+")*("+weight+")").c_str(), "goff");
    cout<<"number analysis: "<<analysis->Integral()*lumi<<endl;
    analysis->Scale(1./analysis->GetEntries());
    analysis->SetLineColor(color[i]);
    analysis->Draw("HIST SAME");
    legend->AddEntry(analysis, lName[i].c_str(), "l");
  }
  legend->Draw();
  t_c->SaveAs((folder_name+"llg_res_El.pdf").c_str());

}

