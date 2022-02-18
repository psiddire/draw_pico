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
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_data/2017/data/merged_zgmc_ll_el/*.root");

  string weight = "1"; // w_lumi

  float lumi = 41.5;
  string folder_name = "plots/Res/";

  string title = "";
  string el1b = "(abs(el_eta[ll_i1[0]]) < 0.8)";
  string el1o = "(abs(el_eta[ll_i1[0]]) > 0.8 && abs(el_eta[ll_i1[0]]) < 1.5)";
  string el1e = "(abs(el_eta[ll_i1[0]]) > 1.5 && abs(el_eta[ll_i1[0]]) < 2.5)";
  string el2b = "(abs(el_eta[ll_i2[0]]) < 0.8)";
  string el2o = "(abs(el_eta[ll_i2[0]]) > 0.8 && abs(el_eta[ll_i2[0]]) < 1.5)";
  string el2e = "(abs(el_eta[ll_i2[0]]) > 1.5 && abs(el_eta[ll_i2[0]]) < 2.5)";

  string bb = el1b+"&&"+el2b;
  string bo = el1b+"&&"+el2o;
  string be = el1b+"&&"+el2e;
  string ob = el1o+"&&"+el2b;
  string oo = el1o+"&&"+el2o;
  string oe = el1o+"&&"+el2e;
  string eb = el1e+"&&"+el2b;
  string eo = el1e+"&&"+el2o;
  string ee = el1e+"&&"+el2e;

  TFile *f = new TFile("ElEl.root", "RECREATE");
  TList *l = new TList();

  string cut[9] = {bb, bo, be, ob, oo, oe, eb, eo, ee};

  string prefix_res[9] = {"ll_res_bb", "ll_res_bo", "ll_res_be", "ll_res_ob", "ll_res_oo", "ll_res_oe", "ll_res_eb", "ll_res_eo", "ll_res_ee"};
  string x_axis_res = "sqrt(ll_dml1[0]^2 + ll_dml2[0]^2)";
  for (int i = 0; i < 9; i++) {
    TH1D* analysis = new TH1D(prefix_res[i].c_str(), title.c_str(), 100, 0., 5.0);
    tree->Draw((x_axis_res+">>"+prefix_res[i]).c_str(), ("("+cut[i]+")*("+weight+")").c_str(), "goff");
    cout<<"number analysis: "<<analysis->Integral()*lumi<<endl;
    l->Add(analysis);
  }

  string prefix_m[9] = {"ll_m_bb", "ll_m_bo", "ll_m_be", "ll_m_ob", "ll_m_oo", "ll_m_oe", "ll_m_eb", "ll_m_eo", "ll_m_ee"};
  string x_axis_m = "ll_m[0]";
  for (int i = 0; i < 9; i++) {
    TH1D* analysis = new TH1D(prefix_m[i].c_str(), title.c_str(), 60, 60, 120);
    tree->Draw((x_axis_m+">>"+prefix_m[i]).c_str(), ("("+cut[i]+")*("+weight+")").c_str(), "goff");
    cout<<"number analysis: "<<analysis->Integral()*lumi<<endl;
    l->Add(analysis);
  }

  l->Write("histlist", TObject::kSingleKey);
  f->ls();
}

