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
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/*.root");

  string weight = "1"; // w_lumi

  float lumi = 41.5;
  string folder_name = "plots/Res/";

  string title = "";
  string mu1b = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 0.8)";
  string mu1o = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) > 0.8 && abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 1.5)";
  string mu1e = "(abs(mu_eta[ll_i1[llphoton_ill[0]]]) > 1.5 && abs(mu_eta[ll_i1[llphoton_ill[0]]]) < 2.5)";
  string mu2b = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 0.8)";
  string mu2o = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) > 0.8 && abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 1.5)";
  string mu2e = "(abs(mu_eta[ll_i2[llphoton_ill[0]]]) > 1.5 && abs(mu_eta[ll_i2[llphoton_ill[0]]]) < 2.5)";
  string phb  = "(abs(photon_eta[llphoton_iph[0]]) < 1.5)";
  string phe  = "(abs(photon_eta[llphoton_iph[0]]) > 1.5)";

  string bbb = mu1b+"&&"+mu2b+"&&"+phb;
  string bob = mu1b+"&&"+mu2o+"&&"+phb;
  string beb = mu1b+"&&"+mu2e+"&&"+phb;
  string obb = mu1o+"&&"+mu2b+"&&"+phb;
  string oob = mu1o+"&&"+mu2o+"&&"+phb;
  string oeb = mu1o+"&&"+mu2e+"&&"+phb;
  string ebb = mu1e+"&&"+mu2b+"&&"+phb;
  string eob = mu1e+"&&"+mu2o+"&&"+phb;
  string eeb = mu1e+"&&"+mu2e+"&&"+phb;
  string bbe = mu1b+"&&"+mu2b+"&&"+phe;
  string boe = mu1b+"&&"+mu2o+"&&"+phe;
  string bee = mu1b+"&&"+mu2e+"&&"+phe;
  string obe = mu1o+"&&"+mu2b+"&&"+phe;
  string ooe = mu1o+"&&"+mu2o+"&&"+phe;
  string oee = mu1o+"&&"+mu2e+"&&"+phe;
  string ebe = mu1e+"&&"+mu2b+"&&"+phe;
  string eoe = mu1e+"&&"+mu2o+"&&"+phe;
  string eee = mu1e+"&&"+mu2e+"&&"+phe;

  TFile *f = new TFile("MuMuSig.root", "RECREATE");
  TList *l = new TList();

  string cut[18] = {bbb, bob, beb, obb, oob, oeb, ebb, eob, eeb, bbe, boe, bee, obe, ooe, oee, ebe, eoe, eee};

  string prefix_res[18] = {"llg_res_bbb", "llg_res_bob", "llg_res_beb", "llg_res_obb", "llg_res_oob", "llg_res_oeb", "llg_res_ebb", "llg_res_eob", "llg_res_eeb",
			   "llg_res_bbe", "llg_res_boe", "llg_res_bee", "llg_res_obe", "llg_res_ooe", "llg_res_oee", "llg_res_ebe", "llg_res_eoe", "llg_res_eee"};
  string x_axis_res = "sqrt(llphoton_dml1[0]^2 + llphoton_dml2[0]^2 + llphoton_dmph[0]^2)";
  for (int i = 0; i < 18; i++) {
    TH1D* analysis = new TH1D(prefix_res[i].c_str(), title.c_str(), 100, 0., 5.0);
    tree->Draw((x_axis_res+">>"+prefix_res[i]).c_str(), ("("+cut[i]+")*("+weight+")").c_str(), "goff");
    cout<<"number analysis: "<<analysis->Integral()*lumi<<endl;
    l->Add(analysis);
  }

  string prefix_m[18] = {"llg_m_bbb", "llg_m_bob", "llg_m_beb", "llg_m_obb", "llg_m_oob", "llg_m_oeb", "llg_m_ebb", "llg_m_eob", "llg_m_eeb",
			 "llg_m_bbe", "llg_m_boe", "llg_m_bee", "llg_m_obe", "llg_m_ooe", "llg_m_oee", "llg_m_ebe", "llg_m_eoe", "llg_m_eee"};
  string x_axis_m = "llphoton_m[0]";
  for (int i = 0; i < 18; i++) {
    TH1D* analysis = new TH1D(prefix_m[i].c_str(), title.c_str(), 50, 100, 150);
    tree->Draw((x_axis_m+">>"+prefix_m[i]).c_str(), ("("+cut[i]+")*("+weight+")").c_str(), "goff");
    cout<<"number analysis: "<<analysis->Integral()*lumi<<endl;
    l->Add(analysis);
  }

  l->Write("histlist", TObject::kSingleKey);
  f->ls();
}

