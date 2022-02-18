#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <bitset>
#include <unistd.h>
#include <getopt.h>
#include "TError.h"
#include "TColor.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/utilities.hpp"
#include "zgamma/zg_utilities.hpp"
using namespace std;
using namespace PlotOptTypes;
using namespace ZgUtilities;

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type sig  = Process::Type::signal;
  string bfolder("/net/cms17/cms17r0/pico/");
  string ul_path(bfolder+"NanoAODv2/zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string rereco_path(bfolder+"NanoAODv7/zgamma_signal/2017/signal/merged_zgmc_llg/");

  NamedFunc nllphoton(       "nllphoton", [](const Baby &b) -> NamedFunc::ScalarType{ return b.nllphoton(); });
  NamedFunc llg_costheta( "llg_costheta", [](const Baby &b) -> NamedFunc::ScalarType{ return b.llphoton_costheta()->at(0); });
  NamedFunc llg_bigTheta( "llg_bigTheta", [](const Baby &b) -> NamedFunc::ScalarType{ return b.llphoton_cosTheta()->at(0); });
  NamedFunc phi(                   "phi", [](const Baby &b) -> NamedFunc::ScalarType{ return Getphi(b); });
  NamedFunc hpt_hm(           "h_ptdh_m", [](const Baby &b) -> NamedFunc::ScalarType{ return b.llphoton_pt()->at(0)/b.llphoton_m()->at(0); });
  NamedFunc gpt_hm(           "g_ptdh_m", [](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pt()->at(0)/b.llphoton_m()->at(0); });
  NamedFunc mll_mllg(         "mll_mllg", [](const Baby &b) -> NamedFunc::ScalarType{ return b.ll_m()->at(0)+b.llphoton_m()->at(0); });
  NamedFunc p_pte(          "photon_pte", [](const Baby &b) -> NamedFunc::ScalarType{ return b.photon_pterr()->at(0)/b.photon_pt()->at(0); });
  NamedFunc p_eb(            "photon_EB", [](const Baby &b) -> NamedFunc::ScalarType{ return abs(b.photon_eta()->at(0)) < 1.4442; });
  NamedFunc p_ee(            "photon_EE", [](const Baby &b) -> NamedFunc::ScalarType{ return abs(b.photon_eta()->at(0)) > 1.566; });
  NamedFunc g_drmax(      "photon_drmax", [](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 photon = AssignGamma(b).Vect();
    TVector3 l1     = AssignL1(b).Vect();
    TVector3 l2     = AssignL2(b).Vect();
    return max(photon.DeltaR(l1),photon.DeltaR(l2));
  });
  NamedFunc sysbal("sysbal",[](const Baby &b) -> NamedFunc::ScalarType{
      TVector3 z = AssignZ(b).Vect();
      TVector3 g = AssignGamma(b).Vect();
      TVector3 j1, j2;
      int i1(-1), i2(-1);
      for(size_t ij(0); ij < b.jet_pt()->size(); ij++)
        if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
        else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
      j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
      j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
      return (z+g+j1+j2).Pt()/(z.Pt()+g.Pt()+j1.Pt()+j2.Pt());
    });
  NamedFunc jjllg_dphi("jjllg_dphi",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 h = AssignH(b).Vect();
    TVector3 jj;
    TVector3 j1, j2;
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
    j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
    jj = j1 + j2;
    return h.DeltaPhi(jj);
  });
  NamedFunc photon_jdr("photon_jdr",[](const Baby &b) -> NamedFunc::ScalarType{ 
    TVector3 g = AssignGamma(b).Vect();
    TVector3 j1, j2;
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    j1.SetPtEtaPhi(b.jet_pt()->at(i1),b.jet_eta()->at(i1),b.jet_phi()->at(i1));
    j2.SetPtEtaPhi(b.jet_pt()->at(i2),b.jet_eta()->at(i2),b.jet_phi()->at(i2));
    return min(g.DeltaR(j1),g.DeltaR(j2));
  });
  NamedFunc photon_zep("photon_zep",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int i1(-1), i2(-1);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) 
      if(b.jet_isgood()->at(ij) && i1 == -1) i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) i2 = ij;
    return abs(b.photon_eta()->at(0) - (b.jet_eta()->at(i1)+b.jet_eta()->at(i2))/2);
  });
  NamedFunc pTt2("pTt2",[](const Baby &b) -> NamedFunc::ScalarType{
    TVector3 g = AssignGamma(b).Vect();
    TVector3 h = AssignH(b).Vect();
    TVector3 z = AssignZ(b).Vect();
    g.SetZ(0); h.SetZ(0); z.SetZ(0);
    return h.Cross(z-g).Mag()/h.Mag();
  }); 
  NamedFunc jet1_pt("jet1_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++)
      if(b.jet_isgood()->at(ij) && i1 == -1) {
	i1 = ij;
        pt = b.jet_pt()->at(ij);
      }
    return pt;
  });
  NamedFunc jet2_pt("jet2_pt",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1), i2(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++)
      if(b.jet_isgood()->at(ij) && i1 == -1)  i1 = ij;
      else if(b.jet_isgood()->at(ij) && i2 == -1) {
	i2 = ij;
	pt = b.jet_pt()->at(ij);
      }
    return pt;
  });
  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight*100;
  });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc mass_cuts("ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110");
  NamedFunc baseline("nphoton > 0 && nll > 0");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15 && "
			   "el_sig[ll_i1[0]] && el_sig[ll_i2[0]]",
                           "ll_lepid[0] == 13 && mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10 && "
			   "mu_sig[ll_i1[0]] && mu_sig[ll_i2[0]]"};
  NamedFunc pho("photon_pt[0] > 15 && photon_drmin[0] > 0.4 && photon_sig[0]");

  // NamedFunc fullsel_el("ll_lepid[llphoton_ill[j]] == 11 && el_pt[ll_i1[llphoton_ill[j]]] > 25 && el_pt[ll_i2[llphoton_ill[j]]] > 15 && llphoton_m[j]+ll_m[j]>=185 && llphoton_m[j] > 100 && llphoton_m[j] < 180 && photon_pt[llphoton_iph[j]]/llphoton_m[j] >= 15./110");
  // NamedFunc fullsel_mu("ll_lepid[llphoton_ill[j]] == 13 && mu_pt[ll_i1[llphoton_ill[j]]] > 20 && mu_pt[ll_i2[llphoton_ill[j]]] > 10 && llphoton_m[j]+ll_m[j]>=185 && llphoton_m[j] > 100 && llphoton_m[j] < 180 && photon_pt[llphoton_iph[j]]/llphoton_m[j] >= 15./110");

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::none)
          .YTitleOffset(1.75)
          .AutoYAxis(false)
          .UseCMYK(false)
          .LeftMargin(0.17)
          .LegendColumns(1)
          .CanvasWidth(800)
          .CanvasHeight(800)
          .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};

  auto proc_UL     = Process::MakeShared<Baby_pico>("HToZ#gamma (UL)", sig, 
                      TColor::GetColor("#ff0000"), {ul_path+"*.root"}, trigs);
  auto proc_Rereco = Process::MakeShared<Baby_pico>("HToZ#gamma (Rereco)", sig, 
                      TColor::GetColor("#0000ff"), {rereco_path+"*.root"}, trigs);

  proc_UL->SetLineWidth(3);
  proc_Rereco->SetLineWidth(3);

  vector<shared_ptr<Process>> procs = {proc_UL, proc_Rereco};

  PlotMaker pm;
  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && pho && mass_cuts;
    pm.Push<Hist1D>(Axis(20, 100, 180, "llphoton_m[0]", "m_{ll#gamma} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"llphoton_m");
    pm.Push<Hist1D>(Axis(20, 50, 150, "ll_m[0]", "m_{ll} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"ll_m");
    pm.Push<Hist1D>(Axis(40, 0, 2, "llphoton_pt[0]/llphoton_m[0]", "p_{T,ll#gamma}/m_{ll#gamma}",{}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"pToverMass");
    // Kinematic MVA training variables
    pm.Push<Hist1D>(Axis(80, -4, 4, phi, "#phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"phi");
    pm.Push<Hist1D>(Axis(80, -4, 4, "llphoton_psi[0]", "#phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"llphoton_phi");
    pm.Push<Hist1D>(Axis(20, -1, 1, llg_bigTheta, "cos(#Theta)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"llphoton_cosTheta");
    pm.Push<Hist1D>(Axis(20, -1, 1, llg_costheta, "cos(#theta)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"llphoton_theta");
    if(i == 0){
      pm.Push<Hist1D>(Axis(15, 0, 150, "el_pt[ll_i1[0]]", "Leading electron p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "el_eta[ll_i1[0]]", "Leading electron #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "el_phi[ll_i1[0]]", "Leading electron #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el1_phi");
      pm.Push<Hist1D>(Axis(15, 0, 150, "el_pt[ll_i2[0]]", "Trailing electron p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "el_eta[ll_i2[0]]", "Trailing electron #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "el_phi[ll_i2[0]]", "Trailing electron #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"el2_phi");
    }
    else if(i == 1){
      pm.Push<Hist1D>(Axis(15, 0, 150, "mu_pt[ll_i1[0]]", "Leading muon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "mu_eta[ll_i1[0]]", "Leading muon #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "mu_phi[ll_i1[0]]", "Leading muon #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu1_phi");
      pm.Push<Hist1D>(Axis(15, 0, 150, "mu_pt[ll_i2[0]]", "Trailing muon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_pt");
      pm.Push<Hist1D>(Axis(50, -2.5, 2.5, "mu_eta[ll_i2[0]]", "Trailing muon #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_eta");
      pm.Push<Hist1D>(Axis(80, -4, 4, "mu_phi[ll_i2[0]]", "Trailing muon #phi", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"mu2_phi");
    }
    // Photon variables
    pm.Push<Hist1D>(Axis(20, -1, 1, "photon_idmva[0]", "Photon ID MVA", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_idmva");
    pm.Push<Hist1D>(Axis(20, 0, 2, "photon_r9[0]", "Photon R9", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_r9");
    pm.Push<Hist1D>(Axis(15, 0, 150, "photon_pt[0]", "Photon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_pt");
    pm.Push<Hist1D>(Axis(20, -2.5, 2.5, "photon_eta[0]", "Photon #eta", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_eta");
    pm.Push<Hist1D>(Axis(35, -0.4, 1, "photon_idmva[0]", "Photon IDMVA (EB)", {}), selection && p_eb, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_idmva_eb");
    pm.Push<Hist1D>(Axis(20, -0.6, 1, "photon_idmva[0]", "Photon IDMVA (EE)", {}), selection && p_ee, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_idmva_ee");
    pm.Push<Hist1D>(Axis(20, 0.4, 5.4, "photon_drmin[0]", "Min #DeltaR(#gamma,l)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_drmin");
    pm.Push<Hist1D>(Axis(20, 0.4, 5.4, g_drmax, "Max #DeltaR(#gamma,l)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_drmax");
    pm.Push<Hist1D>(Axis(40, 0.00, 0.20, "photon_pterr[0]/photon_pt[0]", "#sigma_{#gamma}", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"photon_pTerroverpT");
    // Dijet MVA training variables
    selection = "njet >= 2" && baseline && lep.at(i) && pho && mass_cuts;
    pm.Push<Hist1D>(Axis(12, 30, 210, jet1_pt, "Leading jet p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_jet1_pt");
    pm.Push<Hist1D>(Axis(12, 30, 210, jet2_pt, "Subleading jet p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_jet2_pt");
    pm.Push<Hist1D>(Axis(10, 0, 1, sysbal, "System balance",{}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_sysbal");
    pm.Push<Hist1D>(Axis(12, 0, 9, "dijet_deta", "#Delta#eta(j_{1},j_{2})", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_deta");
    pm.Push<Hist1D>(Axis(10, 0, 3, "dijet_dphi", "#Delta#phi(j_{1},j_{2})", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_dphi");
    pm.Push<Hist1D>(Axis(10, 0, 3, jjllg_dphi, "#Delta#phi(Z#gamma,jj)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_jjllg_dphi");
    pm.Push<Hist1D>(Axis(15, 0.4, 4.4, photon_jdr, "#DeltaR(#gamma,j)", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_photon_jdr");
    pm.Push<Hist1D>(Axis(20, 0, 140, pTt2, "p_{Tt} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_pTt2");
    pm.Push<Hist1D>(Axis(8, 0, 6, photon_zep, "#gamma zeppenfeld",{}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_photon_zeppenfeld");
  }
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

