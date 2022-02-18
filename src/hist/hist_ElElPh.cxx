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
  string el1b = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) < 0.8)";
  string el1o = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) > 0.8 && abs(el_eta[ll_i1[llphoton_ill[0]]]) < 1.5)";
  string el1e = "(abs(el_eta[ll_i1[llphoton_ill[0]]]) > 1.5 && abs(el_eta[ll_i1[llphoton_ill[0]]]) < 2.5)";
  string el2b = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) < 0.8)";
  string el2o = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) > 0.8 && abs(el_eta[ll_i2[llphoton_ill[0]]]) < 1.5)";
  string el2e = "(abs(el_eta[ll_i2[llphoton_ill[0]]]) > 1.5 && abs(el_eta[ll_i2[llphoton_ill[0]]]) < 2.5)";
  string phb  = "(abs(photon_eta[llphoton_iph[0]]) < 1.5)";
  string phe  = "(abs(photon_eta[llphoton_iph[0]]) > 1.5)";

  string bbb = el1b+"&&"+el2b+"&&"+phb;
  string bob = el1b+"&&"+el2o+"&&"+phb;
  string beb = el1b+"&&"+el2e+"&&"+phb;
  string obb = el1o+"&&"+el2b+"&&"+phb;
  string oob = el1o+"&&"+el2o+"&&"+phb;
  string oeb = el1o+"&&"+el2e+"&&"+phb;
  string ebb = el1e+"&&"+el2b+"&&"+phb;
  string eob = el1e+"&&"+el2o+"&&"+phb;
  string eeb = el1e+"&&"+el2e+"&&"+phb;
  string bbe = el1b+"&&"+el2b+"&&"+phe;
  string boe = el1b+"&&"+el2o+"&&"+phe;
  string bee = el1b+"&&"+el2e+"&&"+phe;
  string obe = el1o+"&&"+el2b+"&&"+phe;
  string ooe = el1o+"&&"+el2o+"&&"+phe;
  string oee = el1o+"&&"+el2e+"&&"+phe;
  string ebe = el1e+"&&"+el2b+"&&"+phe;
  string eoe = el1e+"&&"+el2o+"&&"+phe;
  string eee = el1e+"&&"+el2e+"&&"+phe;

  TFile *f = new TFile("ElElSig.root", "RECREATE");
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

