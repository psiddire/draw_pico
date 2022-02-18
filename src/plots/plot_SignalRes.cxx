#include "TChain.h"
#include "TCanvas.h"
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


void draw_2D(TChain * tree, int nbins, float xmin, float xmax, float ymin, float ymax, float lumi, string x_axis, string x_name, string y_axis, string y_name, string cut, string weight, string folder_name, string prefix) {

  // TH2D analysis(prefix.c_str(), prefix.c_str(), nbins, xmin, xmax, nbins, ymin, ymax);
  TH2D* analysis = new TH2D(prefix.c_str(), prefix.c_str(), nbins, xmin, xmax, nbins, ymin, ymax);
  tree->Draw((y_axis+":"+x_axis+">>"+prefix).c_str(), ("("+cut+")*("+weight+")").c_str(), "goff");

  cout<<"number analysis:  "<<analysis->Integral()*lumi<<endl;

  gStyle->SetOptStat(0);

  TCanvas * t_c = 0;
  t_c = newCanvas();
  analysis->Draw("colz");
  analysis->GetXaxis()->SetTitle(x_name.c_str());
  analysis->GetYaxis()->SetTitle(y_name.c_str());
  analysis->GetYaxis()->SetTitleOffset(1.3);
  t_c->SaveAs((folder_name+prefix+".pdf").c_str());
}


int main() {

  TChain * tree = new TChain("tree");
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/*.root");

  int nbins = 100;
  float xmin = -2.5;
  float xmax = 2.5;
  float ymin = 0.0;
  float ymax = 0.5;
  float ymax_mu = 0.1;
  float ymax_llg = 5;
  float lumi = 41.5;
  string folder_name = "plots/";

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

  string el1_res = "(el_energyErr[ll_i1[llphoton_ill[0]]] / el_pt[ll_i1[llphoton_ill[0]]])";
  string el2_res = "(el_energyErr[ll_i2[llphoton_ill[0]]] / el_pt[ll_i2[llphoton_ill[0]]])";
  string mu1_res = "(mu_ptErr[ll_i1[llphoton_ill[0]]] / mu_pt[ll_i1[llphoton_ill[0]]])";
  string mu2_res = "(mu_ptErr[ll_i2[llphoton_ill[0]]] / mu_pt[ll_i2[llphoton_ill[0]]])";
  string pho_res = "(photon_pterr[llphoton_iph[0]] / photon_pt[llphoton_iph[0]])";
  string llg_err = "sqrt(llphoton_dml1[0]^2 + llphoton_dml2[0]^2 + llphoton_dmph[0]^2)";
  string llg_res = "sqrt(llphoton_dml1[0]^2 + llphoton_dml2[0]^2 + llphoton_dmph[0]^2) / llphoton_m[0]";

  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    if(i == 0){
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax, lumi, "el_eta[ll_i1[llphoton_ill[0]]]", "Electron 1 #eta", el1_res, "Electron 1 #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_eta_el1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax, lumi, "el_eta[ll_i2[llphoton_ill[0]]]", "Electron 2 #eta", el2_res, "Electron 2 #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_eta_el2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", pho_res, "Photon #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_eta_elph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "el_eta[ll_i1[llphoton_ill[0]]]", "Electron 1 #eta", llg_err, "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llg_el1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "el_eta[ll_i2[llphoton_ill[0]]]", "Electron 2 #eta", llg_err, "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llg_el2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", llg_err, "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llg_elph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "el_eta[ll_i1[llphoton_ill[0]]]", "Electron 1 #eta", llg_res, "ll#gamma #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_llgres_el1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "el_eta[ll_i2[llphoton_ill[0]]]", "Electron 2 #eta", llg_res, "ll#gamma #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_llgres_el2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", llg_res, "ll#gamma #sigma/p_{T}", el_analysis, weight, folder_name+lepName[i], "res_llgres_elph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "el_eta[ll_i1[llphoton_ill[0]]]", "Electron 1 #eta", "llphoton_dml1[0]", "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llgerr_el1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "el_eta[ll_i2[llphoton_ill[0]]]", "Electron 2 #eta", "llphoton_dml2[0]", "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llgerr_el2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", "llphoton_dmph[0]", "ll#gamma #sigma [GeV]", el_analysis, weight, folder_name+lepName[i], "res_llgerr_elph");
    }
    else if(i == 1){
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "mu_eta[ll_i1[llphoton_ill[0]]]", "Muon 1 #eta", mu1_res, "Muon 1 #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_eta_mu1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "mu_eta[ll_i2[llphoton_ill[0]]]", "Muon 2 #eta", mu2_res, "Muon 2 #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_eta_mu2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", pho_res, "Photon #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_eta_muph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "mu_eta[ll_i1[llphoton_ill[0]]]", "Muon 1 #eta", llg_err, "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llg_mu1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "mu_eta[ll_i2[llphoton_ill[0]]]", "Muon 2 #eta", llg_err, "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llg_mu2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", llg_err, "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llg_muph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "mu_eta[ll_i1[llphoton_ill[0]]]", "Muon 1 #eta", llg_res, "ll#gamma #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_llgres_mu1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "mu_eta[ll_i2[llphoton_ill[0]]]", "Muon 2 #eta", llg_res, "ll#gamma #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_llgres_mu2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_mu, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", llg_res, "ll#gamma #sigma/p_{T}", mu_analysis, weight, folder_name+lepName[i], "res_llgres_muph");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "mu_eta[ll_i1[llphoton_ill[0]]]", "Muon 1 #eta", "llphoton_dml1[0]", "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llgerr_mu1");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "mu_eta[ll_i2[llphoton_ill[0]]]", "Muon 2 #eta", "llphoton_dml2[0]", "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llgerr_mu2");
      draw_2D(tree, nbins, xmin, xmax, ymin, ymax_llg, lumi, "photon_eta[llphoton_iph[0]]", "Photon #eta", "llphoton_dmph[0]", "ll#gamma #sigma [GeV]", mu_analysis, weight, folder_name+lepName[i], "res_llgerr_muph");
    }
  }
}

