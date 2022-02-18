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

using namespace std;
using namespace PlotOptTypes;

const long double PI = acos(-1.L);

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*PI);
  return dphi>PI ? 2.L*PI-dphi : dphi;
}

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

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv7/");
  string foldersig(bfolder+"zgamma_mm/2017/mm/unskimmed/*");

  NamedFunc llphoton_cuts("three body invariant mass cuts",[](const Baby &b) -> NamedFunc::VectorType{
    vector<double> llphoton_cuts_;
    llphoton_cuts_.push_back(-1.0); llphoton_cuts_.push_back(-1.0); llphoton_cuts_.push_back(-1.0); llphoton_cuts_.push_back(-1.0);
    TLorentzVector l1, l2, dilep, ph, llg, selZ;
    double mass(-1);
    int l1idx(-1), l2idx(-1);
    if(b.nel() >= 2) {
      for(size_t iel1(0); iel1 < b.el_pt()->size(); iel1++)
	for(size_t iel2(iel1+1); iel2 < b.el_pt()->size(); iel2++)
	  if(b.el_pt()->at(iel1)*b.el_sig()->at(iel1) > 25 &&
	     b.el_pt()->at(iel2)*b.el_sig()->at(iel2) > 15 &&
	     b.el_charge()->at(iel1)*b.el_charge()->at(iel2) == -1) {
	    l1.SetPtEtaPhiM(b.el_pt()->at(iel1), b.el_eta()->at(iel1), b.el_phi()->at(iel1), 0.000511);
	    l2.SetPtEtaPhiM(b.el_pt()->at(iel2), b.el_eta()->at(iel2), b.el_phi()->at(iel2), 0.000511);
	    dilep = l1 + l2;
	    if(mass == -1 || abs(dilep.M() - 91.1876) < abs(mass - 91.1876)) {
	      mass = dilep.M();
	      selZ = dilep;
	      l1idx = iel1;
	      l2idx = iel2;
	    }
	  }
      if (l1idx!=-1 && l2idx!=-1) {
	llphoton_cuts_[0] = 1.0;
	if(mass > 50) {
	  llphoton_cuts_[1] = 1.0;
	  if(b.photon_pt()->size() > 0) {
	    for(size_t iph(0); iph < b.photon_pt()->size(); iph++) {
	      if(b.photon_pt()->at(iph) > 15 && b.photon_sig()->at(iph)) {
		ph.SetPtEtaPhiM(b.photon_pt()->at(iph), b.photon_eta()->at(iph), b.photon_phi()->at(iph), 0.0);
		llphoton_cuts_[2] = 1.0;
		llg = selZ + ph;
		if(llg.M() > 100 && llg.M() < 180)
		  if(llg.M() + mass >= 185)
		    if(ph.Pt()/llg.M() >= 15./110)
		      llphoton_cuts_[3] = 1.0;
	      }
	    }
	  }
	}
      }
    }
    else if(b.nmu() >= 2) {
      for(size_t imu1(0); imu1 < b.mu_pt()->size(); imu1++)
	for(size_t imu2(imu1+1); imu2 < b.mu_pt()->size(); imu2++)
	  if(b.mu_pt()->at(imu1)*b.mu_sig()->at(imu1) > 20 &&
	     b.mu_pt()->at(imu2)*b.mu_sig()->at(imu2) > 10 &&
	     b.mu_charge()->at(imu1)*b.mu_charge()->at(imu2) == -1) {
	    l1.SetPtEtaPhiM(b.mu_pt()->at(imu1), b.mu_eta()->at(imu1), b.mu_phi()->at(imu1), 0.10566);
	    l2.SetPtEtaPhiM(b.mu_pt()->at(imu2), b.mu_eta()->at(imu2), b.mu_phi()->at(imu2), 0.10566);
	    dilep = l1 + l2;
            if(mass == -1 || abs(dilep.M() - 91.1876) < abs(mass - 91.1876)) {
              mass = dilep.M();
	      selZ = dilep;
              l1idx = imu1;
              l2idx = imu2;
            }
          }
      if (l1idx!=-1 && l2idx!=-1) {
	llphoton_cuts_[0] = 1.0;
	if(mass > 50) {
	  llphoton_cuts_[1] = 1.0;
	  if(b.photon_pt()->size() > 0) {
	    for(size_t iph(0); iph < b.photon_pt()->size(); iph++) {
	      if(b.photon_pt()->at(iph) > 15 && b.photon_sig()->at(iph)) {
		ph.SetPtEtaPhiM(b.photon_pt()->at(iph), b.photon_eta()->at(iph), b.photon_phi()->at(iph), 0.0);
		llphoton_cuts_[2] = 1.0;
		llg = selZ + ph;
		if(llg.M() > 100 && llg.M() < 180)
		  if(llg.M() + mass >= 185)
		    if(ph.Pt()/llg.M() >= 15./110)
		      llphoton_cuts_[3] = 1.0;
	      }
	    }
	  }
	}
      }
    }
    return llphoton_cuts_;
  });

  NamedFunc wgt("w_lumi",[](const Baby &b) -> NamedFunc::ScalarType{
    return b.w_lumi();
  });

  // Trigger
  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  auto proc_el_sig = Process::MakeShared<Baby_pico>("Sig e",  Process::Type::signal, kBlack, {foldersig+"*.root"}, el_trigs);
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("#rm{Sig}#; #mu",Process::Type::signal, kBlack, {foldersig+"*.root"}, mu_trigs);

  vector<shared_ptr<Process> > samples = {proc_el_sig, proc_mu_sig};

  vector<NamedFunc> cutflow;
  cutflow.push_back("1");
  cutflow.push_back((el_trigs || mu_trigs));
  cutflow.push_back((el_trigs || mu_trigs) && llphoton_cuts[0.] > 0.0);
  cutflow.push_back((el_trigs || mu_trigs) && llphoton_cuts[1] > 0.0);
  cutflow.push_back((el_trigs || mu_trigs) && llphoton_cuts[2] > 0.0);
  cutflow.push_back((el_trigs || mu_trigs) && llphoton_cuts[3] > 0.0);

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::overflow)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> ops = {log_lumi, lin_lumi};
  PlotMaker pm;
  pm.Push<Table>("AN_hmm_cuts_cutflow", vector<TableRow>{
      TableRow("Total number of events",                cutflow.at(0),0,0,wgt),
      TableRow("High level trigger",                    cutflow.at(1),0,0,wgt),
      TableRow("lepton selections",                     cutflow.at(2),0,0,wgt),
      TableRow("$m_{ll} > $ 50 GeV",                    cutflow.at(3),0,0,wgt),
      TableRow("photon ID",                             cutflow.at(4),0,0,wgt),
      TableRow("three body invariant mass related cut", cutflow.at(5),0,0,wgt),
      },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}
