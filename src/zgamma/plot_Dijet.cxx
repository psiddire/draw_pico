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

const long double PI = acos(-1.L);

long double AddInQuadrature(long double x, long double y){
  if(fabs(y)>fabs(x)){
    const long double temp = y;
    y=x;
    x=temp;
  }
  if(x==0.) return y;
  const long double rat=y/x;
  return fabs(x)*sqrt(1.0L+rat*rat);
}

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*PI);
  return dphi>PI ? 2.L*PI-dphi : dphi;
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

int main() {
  gErrorIgnoreLevel = 6000;
  Palette colors("txt/colors.txt","default");
  Process::Type back = Process::Type::background;
  Process::Type sig  = Process::Type::signal;
  Process::Type data = Process::Type::data;
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/zgamma_mc_ul/2017/");
  string mc_path(bfolder+"mc/merged_zgmc_llg/");
  string sig_path(bfolder+"../../zgamma_signal_ul/2017/signal/merged_zgmc_llg/");
  string data_path(bfolder+"../../zgamma_data/2017/data/merged_zgmc_llg/");


  NamedFunc jet_pt1("jet_pt1",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) {
      float jeta = b.jet_eta()->at(ij); 
      float jphi = b.jet_phi()->at(ij); 
      float geta = b.photon_eta()->at(0);
      float gphi = b.photon_phi()->at(0);
      if (b.ll_lepid()->at(0) == 11) {
	float eeta1 = b.el_eta()->at(b.ll_i1()->at(0));
	float ephi1 = b.el_phi()->at(b.ll_i1()->at(0));
	float eeta2 = b.el_eta()->at(b.ll_i2()->at(0));
	float ephi2 = b.el_phi()->at(b.ll_i2()->at(0));
	if (dR(jeta, eeta1, jphi, ephi1) > 0.4 && dR(jeta, eeta2, jphi, ephi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i1 == -1) {
	  i1 = ij;
	  pt = b.jet_pt()->at(ij);
	}
      }
      else if (b.ll_lepid()->at(0) == 13) {
        float meta1 = b.mu_eta()->at(b.ll_i1()->at(0));
        float mphi1 = b.mu_phi()->at(b.ll_i1()->at(0));
        float meta2 = b.mu_eta()->at(b.ll_i2()->at(0));
        float mphi2 = b.mu_phi()->at(b.ll_i2()->at(0));
	if (dR(jeta, meta1, jphi, mphi1) > 0.4 && dR(jeta, meta2, jphi, mphi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i1 == -1) {
	  i1 = ij;
	  pt = b.jet_pt()->at(ij);
	}
      }
    }
    return pt;
  });

  NamedFunc jet_pt2("jet_pt2",[](const Baby &b) -> NamedFunc::ScalarType{
    int i1(-1), i2(-1);
    double pt(0);
    for(size_t ij(0); ij < b.jet_pt()->size(); ij++) {
      float jeta = b.jet_eta()->at(ij); 
      float jphi = b.jet_phi()->at(ij); 
      float geta = b.photon_eta()->at(0);
      float gphi = b.photon_phi()->at(0);
      if (b.ll_lepid()->at(0) == 11) {
	float eeta1 = b.el_eta()->at(b.ll_i1()->at(0));
	float ephi1 = b.el_phi()->at(b.ll_i1()->at(0));
	float eeta2 = b.el_eta()->at(b.ll_i2()->at(0));
	float ephi2 = b.el_phi()->at(b.ll_i2()->at(0));
	if (dR(jeta, eeta1, jphi, ephi1) > 0.4 && dR(jeta, eeta2, jphi, ephi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i1 == -1) i1 = ij;
	else if (dR(jeta, eeta1, jphi, ephi1) > 0.4 && dR(jeta, eeta2, jphi, ephi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i2 == -1) {
	  i2 = ij;
	  pt = b.jet_pt()->at(ij);
	}
      }
      else if (b.ll_lepid()->at(0) == 13) {
        float meta1 = b.mu_eta()->at(b.ll_i1()->at(0));
        float mphi1 = b.mu_phi()->at(b.ll_i1()->at(0));
        float meta2 = b.mu_eta()->at(b.ll_i2()->at(0));
        float mphi2 = b.mu_phi()->at(b.ll_i2()->at(0));
	if (dR(jeta, meta1, jphi, mphi1) > 0.4 && dR(jeta, meta2, jphi, mphi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i1 == -1) i1 = ij;
	else if (dR(jeta, meta1, jphi, mphi1) > 0.4 && dR(jeta, meta2, jphi, mphi2) > 0.4 && dR(jeta, geta, jphi, gphi) > 0.4 && i2 == -1) {
          i2 = ij;
          pt = b.jet_pt()->at(ij);
        }
      }
    }
    return pt;
  });

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    if(b.type() >= 1000) {
      weight = weight*b.w_photon();
      if(b.type() >= 200000 && b.type() <= 205000)
	return weight*100;
    }
    return weight;
  });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  NamedFunc trigs(el_trigs || mu_trigs);
  NamedFunc mass_cuts("ll_m[0] > 50 && llphoton_m[0]+ll_m[0]>=185 && llphoton_m[0] > 100 && llphoton_m[0] < 180 && photon_pt[0]/llphoton_m[0] >= 15./110");
  NamedFunc baseline("nphoton > 0 && nll > 0");
  vector<NamedFunc> lep = {"ll_lepid[0] == 11 && ll_dr[0] > 0.4 && "
			   "el_pt[ll_i1[0]] > 25 && el_pt[ll_i2[0]] > 15 && "
			   "el_id[ll_i1[0]] && el_id[ll_i2[0]] && "
			   "el_reliso[ll_i1[0]] < 0.35 && el_reliso[ll_i2[0]] < 0.35 && "
			   "el_sip3d[ll_i1[0]] < 4 && el_sip3d[ll_i2[0]] < 4",
                           "ll_lepid[0] == 13 && ll_dr[0] > 0.4 && "
			   "mu_pt[ll_i1[0]] > 20 && mu_pt[ll_i2[0]] > 10 && "
			   "mu_sig[ll_i1[0]] && mu_sig[ll_i2[0]]"};
  NamedFunc pho("photon_pt[0] > 15 && photon_drmin[0] > 0.4 && photon_sig[0]");

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
          .Bottom(BottomType::ratio)
          .CanvasWidth(900)
          .CanvasHeight(900)
          .FileExtensions({"pdf"});

  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt lin_stack = lin_lumi().Stack(StackType::signal_overlay);
  PlotOpt log_stack = log_lumi().Stack(StackType::signal_overlay);
  vector<PlotOpt> ops = {lin_stack};

  auto proc_smzg = Process::MakeShared<Baby_pico>("SM Z#gamma",       back, 
	 	      TColor::GetColor("#16bac5"), {mc_path+"*ZGToLLG*"}, trigs);
  auto proc_dy   = Process::MakeShared<Baby_pico>("DY",               back, 
		      TColor::GetColor("#ffb400"), {mc_path+"*DYJets*"},  trigs && "stitch_dy");
  auto proc_hzg  = Process::MakeShared<Baby_pico>("HToZ#gamma(x100)", sig, 
                      TColor::GetColor("#ff0000"), {sig_path+"*.root"},   trigs);
  auto proc_data = Process::MakeShared<Baby_pico>("Data",             data, 
                      kBlack,                      {data_path+"*.root"},  trigs);

  proc_smzg->SetLineWidth(1);
  proc_dy->SetLineWidth(1);
  proc_hzg->SetLineWidth(3);
  proc_data->SetMarkerSize(1);

  vector<shared_ptr<Process>> procs = {proc_data, proc_dy, proc_smzg, proc_hzg};

  PlotMaker pm;
  string lepName[2] = {"el/", "mu/"};
  for(int i(0); i < 2; i++) {
    NamedFunc selection = baseline && lep.at(i) && pho && mass_cuts;
    pm.Push<Hist1D>(Axis(12, 30, 210, jet_pt1, "Leading jet p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_jet1_pt");
    pm.Push<Hist1D>(Axis(12, 30, 210, jet_pt2, "Subleading jet p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag(lepName[i]+"2jet_jet2_pt");
  }
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

