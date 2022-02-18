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
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_data/2017/data/merged_zgmc_ll_mu/*.root");

  string weight = "1"; // w_lumi

  float lumi = 41.5;
  string folder_name = "plots/Res/";

  string title = "";
  string mu1b = "(abs(mu_eta[ll_i1[0]]) < 0.9)";
  string mu1o = "(abs(mu_eta[ll_i1[0]]) > 0.9 && abs(mu_eta[ll_i1[0]]) < 1.8)";
  string mu1e = "(abs(mu_eta[ll_i1[0]]) > 1.8 && abs(mu_eta[ll_i1[0]]) < 2.4)";
  string mu2b = "(abs(mu_eta[ll_i2[0]]) < 0.9)";
  string mu2o = "(abs(mu_eta[ll_i2[0]]) > 0.9 && abs(mu_eta[ll_i2[0]]) < 1.8)";
  string mu2e = "(abs(mu_eta[ll_i2[0]]) > 1.8 && abs(mu_eta[ll_i2[0]]) < 2.4)";

  string bb = mu1b+"&&"+mu2b;
  string bo = mu1b+"&&"+mu2o;
  string be = mu1b+"&&"+mu2e;
  string ob = mu1o+"&&"+mu2b;
  string oo = mu1o+"&&"+mu2o;
  string oe = mu1o+"&&"+mu2e;
  string eb = mu1e+"&&"+mu2b;
  string eo = mu1e+"&&"+mu2o;
  string ee = mu1e+"&&"+mu2e;

  TFile *f = new TFile("MuMu.root", "RECREATE");
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

