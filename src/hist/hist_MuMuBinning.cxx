#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "THStack.h"
#include "TStyle.h"
#include <iostream>
#include "TLorentzVector.h"

using namespace std;

int main() {

  TChain * tree = new TChain("tree");
  tree->Add("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_data/2017/data/merged_zgmc_ll_mu/*.root");

  string weight = "1"; // w_lumi

  float lumi = 41.5;
  string folder_name = "plots/Res/";

  cout << tree->GetEntries() << endl;

  vector<float> *ll_m = new vector<float>();
  tree->SetBranchAddress("ll_m", &ll_m);
  vector<int> *ll_i1 = new vector<int>();
  tree->SetBranchAddress("ll_i1", &ll_i1);
  vector<int> *ll_i2 = new vector<int>();
  tree->SetBranchAddress("ll_i2", &ll_i2);
  vector<float> *mu_eta = new vector<float>();
  tree->SetBranchAddress("mu_eta", &mu_eta);
  vector<float> *ll_dml1 = new vector<float>();
  tree->SetBranchAddress("ll_dml1", &ll_dml1);
  vector<float> *ll_dml2 = new vector<float>();
  tree->SetBranchAddress("ll_dml2", &ll_dml2);

  double mu1eta, mu2eta, ll_dm;

  TFile *f = new TFile("MuMuBinning.root", "RECREATE");

  TH1D* mass_err = new TH1D("mass_err", "", 100, 0., 5.0);
  TH1D* mass = new TH1D("mass", "", 60, 60, 120);

  for (int i = 0; i < tree->GetEntries(); i++ ) {
    tree->GetEntry(i);
    if (i % 1000 == 0)
      cout << i << endl;
    mu1eta = abs((*mu_eta)[(*ll_i1)[0]]);
    mu2eta = abs((*mu_eta)[(*ll_i2)[0]]);
    ll_dm = sqrt(pow((*ll_dml1)[0], 2) + pow((*ll_dml2)[0], 2));
    if (mu1eta < 0.9 && mu2eta < 0.9)
      ll_dm = ll_dm * 1.329;
    else if (mu1eta < 0.9 && mu2eta > 0.9 && mu2eta < 1.8)
      ll_dm = ll_dm * 1.322;
    else if (mu1eta < 0.9 && mu2eta > 1.8 && mu2eta < 2.4)
      ll_dm = ll_dm * 1.266;
    else if (mu1eta > 0.9 && mu1eta < 1.8 && mu2eta < 0.9)
      ll_dm = ll_dm * 1.324;
    else if (mu1eta > 0.9 && mu1eta < 1.8 && mu2eta > 0.9 && mu2eta < 1.8)
      ll_dm = ll_dm * 1.307;
    else if (mu1eta > 0.9 && mu1eta < 1.8 && mu2eta > 1.8 && mu2eta < 2.4)
      ll_dm = ll_dm * 1.269;
    else if (mu1eta > 1.8 && mu1eta < 2.4 && mu2eta < 0.9)
      ll_dm = ll_dm * 1.284;
    else if (mu1eta > 1.8 && mu1eta < 2.4 && mu2eta > 0.9 && mu2eta < 1.8)
      ll_dm = ll_dm * 1.308;
    else if (mu1eta > 1.8 && mu1eta < 2.4 && mu2eta > 1.8 && mu2eta < 2.4)
      ll_dm = ll_dm * 1.264;
    mass_err->Fill(ll_dm);
    mass->Fill((*ll_m)[0]);
  }

  f->Write();
  f->Close();

}

