#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
using namespace std;

int main() {
  gErrorIgnoreLevel = 6000;
  string input_path("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string output_path("/homes/psiddire/draw_pico/tree/");

  TFile f("merged_pico_llg_GluGluHToZG_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_llg_nfiles_5.root");
  TTree *input_tree = (TTree*)f.Get("tree");

  Bool_t HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  // bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  // bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8;
  // Int_t nphoton;
  // Int_t nll;
  // Float_t llphoton_m[4];
  // Int_t llphoton_ill[4];
  // Int_t llphoton_iph[4];
  // Float_t ll_m[4];
  // Float_t ll_lepid[4];
  // Int_t ll_i1[4];
  // Int_t ll_i2[4];
  // Float_t el_pt[4];
  // bool el_sig[4];
  // Float_t mu_pt[4];
  // bool mu_sig[4];
  // Float_t photon_pt[4];
  // Float_t photon_drmin[4];
  // bool photon_sig[4];
  // Float_t w_lumi;

  input_tree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  // input_tree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8);
  // input_tree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8);
  // input_tree->SetBranchAddress("nphoton",&nphoton);
  // input_tree->SetBranchAddress("nll",&nll);
  // input_tree->SetBranchAddress("llphoton_m",&llphoton_m);
  // input_tree->SetBranchAddress("llphoton_ill",&llphoton_ill);
  // input_tree->SetBranchAddress("llphoton_iph",&llphoton_iph);
  // input_tree->SetBranchAddress("ll_m",&ll_m);
  // input_tree->SetBranchAddress("ll_lepid",&ll_lepid);
  // input_tree->SetBranchAddress("ll_i1",&ll_i1);
  // input_tree->SetBranchAddress("ll_i2",&ll_i2);
  // input_tree->SetBranchAddress("el_pt",&el_pt);
  // input_tree->SetBranchAddress("el_sig",&el_sig);
  // input_tree->SetBranchAddress("mu_pt",&mu_pt);
  // input_tree->SetBranchAddress("mu_sig",&mu_sig);
  // input_tree->SetBranchAddress("photon_pt",&photon_pt);
  // input_tree->SetBranchAddress("photon_drmin",&photon_drmin);
  // input_tree->SetBranchAddress("photon_sig",&photon_sig);
  // input_tree->SetBranchAddress("w_lumi",&w_lumi);

  Int_t idx = 0;
  Int_t pdx = 0;

  for (Int_t i = 0; i < input_tree->GetEntries(); i++) {
    // Get the values for the i`th event and fill all our local variables
    // that were assigned to TBranches
    input_tree->GetEntry(i);
    
    if (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) {
      printf("OK");
      // if (nll > 0 && nphoton > 0) {
      // 	idx = llphoton_ill[0];
      // 	pdx = llphoton_iph[0];
      // 	if (ll_lepid[idx] == 11) {
      // 	  if (el_pt[ll_i1[idx]] > 25 && el_pt[ll_i2[idx]] > 15) {
      // 	    if (el_sig[ll_i1[idx]] && el_sig[ll_i2[idx]]) {
      // 	      if (photon_pt[pdx] > 15 && photon_drmin[pdx] > 0.4 && photon_sig[pdx]) {
      // 		if (ll_m[idx] > 50 && llphoton_m[0]+ll_m[idx] > 185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[pdx]/llphoton_m[0] >= 15./110) {
      // 		  printf("%f\n", llphoton_m[0]);
      // 		}
      // 	      }
      // 	    }
      // 	  }
      // 	}
      // }
    }
  }

}

