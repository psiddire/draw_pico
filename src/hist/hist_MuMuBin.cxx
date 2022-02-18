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

  TFile *f = new TFile("MuMuBin.root", "RECREATE");

  TH1D* mass_err1 = new TH1D("mass_err1", "", 100, 0., 5.0);
  TH1D* mass1 = new TH1D("mass1", "", 60, 60, 120);
  TH1D* mass_err2 = new TH1D("mass_err2", "", 100, 0., 5.0);
  TH1D* mass2 = new TH1D("mass2", "", 60, 60, 120);
  TH1D* mass_err3 = new TH1D("mass_err3", "", 100, 0., 5.0);
  TH1D* mass3 = new TH1D("mass3", "", 60, 60, 120);
  TH1D* mass_err4 = new TH1D("mass_err4", "", 100, 0., 5.0);
  TH1D* mass4 = new TH1D("mass4", "", 60, 60, 120);
  TH1D* mass_err5 = new TH1D("mass_err5", "", 100, 0., 5.0);
  TH1D* mass5 = new TH1D("mass5", "", 60, 60, 120);
  TH1D* mass_err6 = new TH1D("mass_err6", "", 100, 0., 5.0);
  TH1D* mass6 = new TH1D("mass6", "", 60, 60, 120);
  TH1D* mass_err7 = new TH1D("mass_err7", "", 100, 0., 5.0);
  TH1D* mass7 = new TH1D("mass7", "", 60, 60, 120);
  TH1D* mass_err8 = new TH1D("mass_err8", "", 100, 0., 5.0);
  TH1D* mass8 = new TH1D("mass8", "", 60, 60, 120);
  TH1D* mass_err9 = new TH1D("mass_err9", "", 100, 0., 5.0);
  TH1D* mass9 = new TH1D("mass9", "", 60, 60, 120);
  TH1D* mass_err10 = new TH1D("mass_err10", "", 100, 0., 5.0);
  TH1D* mass10 = new TH1D("mass10", "", 60, 60, 120);
  TH1D* mass_err11 = new TH1D("mass_err11", "", 100, 0., 5.0);
  TH1D* mass11 = new TH1D("mass11", "", 60, 60, 120);
  TH1D* mass_err12 = new TH1D("mass_err12", "", 100, 0., 5.0);
  TH1D* mass12 = new TH1D("mass12", "", 60, 60, 120);
  TH1D* mass_err13 = new TH1D("mass_err13", "", 100, 0., 5.0);
  TH1D* mass13 = new TH1D("mass13", "", 60, 60, 120);

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
    if (ll_dm > 0.4 && ll_dm < 0.8) {
      mass_err1->Fill(ll_dm);
      mass1->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 0.8 && ll_dm < 0.9) {
      mass_err2->Fill(ll_dm);
      mass2->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 0.9 && ll_dm < 1.0) {
      mass_err3->Fill(ll_dm);
      mass3->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.0 && ll_dm < 1.1) {
      mass_err4->Fill(ll_dm);
      mass4->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.1 && ll_dm < 1.2) {
      mass_err5->Fill(ll_dm);
      mass5->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.2 && ll_dm < 1.3) {
      mass_err6->Fill(ll_dm);
      mass6->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.3 && ll_dm < 1.4) {
      mass_err7->Fill(ll_dm);
      mass7->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.4 && ll_dm < 1.6) {
      mass_err8->Fill(ll_dm);
      mass8->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.6 && ll_dm < 1.8) {
      mass_err9->Fill(ll_dm);
      mass9->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 1.8 && ll_dm < 2.0) {
      mass_err10->Fill(ll_dm);
      mass10->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 2.0 && ll_dm < 2.5) {
      mass_err11->Fill(ll_dm);
      mass11->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 2.5 && ll_dm < 3.0) {
      mass_err12->Fill(ll_dm);
      mass12->Fill((*ll_m)[0]);
    }
    else if (ll_dm > 3.0) {
      mass_err13->Fill(ll_dm);
      mass13->Fill((*ll_m)[0]);
    }
  }

  f->Write();
  f->Close();

}

