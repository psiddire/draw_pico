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
  string bfolder("/net/cms17/cms17r0/pico/NanoAODv2/");
  string foldersig(bfolder+"zgamma_gg/2017/gg/unskimmed/*");

  NamedFunc sig_lepid("signal lepton ID",[](const Baby &b) -> NamedFunc::ScalarType{
    int l1idx(1);
    if(b.nel() > 0)
      l1idx = 2;
    else if(b.nmu() > 0)
      l1idx = 3;
    return l1idx;
  });

  NamedFunc llphoton_cuts("three body invariant mass cuts",[](const Baby &b) -> NamedFunc::VectorType{
    vector<double> llphoton_cuts_;
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0);
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0);
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0); 
    llphoton_cuts_.push_back(-1.0);
    TLorentzVector l1, l2, dilep, ph, llg;
    double mass(-1);
    if (b.nel() >= 2) {
      if (b.el_pt()->at(0) > 25 && b.el_pt()->at(1) > 15) {
	llphoton_cuts_[0] = 1.0;
	if (b.el_charge()->at(0)*b.el_charge()->at(1) == -1) {
	  llphoton_cuts_[1] = 1.0;
	  if (b.el_dxy()->at(0) < 0.5 && b.el_dz()->at(0) < 1.0 && b.el_dxy()->at(1) < 0.5 && b.el_dz()->at(1) < 1.0) {
	    llphoton_cuts_[2] = 1.0;
	    if (b.el_sig()->at(0) && b.el_sig()->at(1)) {
	      llphoton_cuts_[3] = 1.0;
	      l1.SetPtEtaPhiM(b.el_pt()->at(0), b.el_eta()->at(0), b.el_phi()->at(0), 0.000511);
	      l2.SetPtEtaPhiM(b.el_pt()->at(1), b.el_eta()->at(1), b.el_phi()->at(1), 0.000511);
	      dilep = l1 + l2;
	      mass = dilep.M();
	      if (mass > 50) {
		llphoton_cuts_[4] = 1.0;
		if (b.photon_pt()->size() > 0) {
		  if (b.photon_pt()->at(0) > 15) {
		    llphoton_cuts_[5] = 1.0;
		    if (b.photon_elveto()->at(0)) {
		      llphoton_cuts_[6] = 1.0;
		      if (b.photon_drmin()->at(0) > 0.4) {
			llphoton_cuts_[7] = 1.0;
			double eta = b.photon_eta()->at(0);
			double mva = b.photon_idmva()->at(0);
			if ((fabs(eta) < 1.4442 && mva > -0.4) || (fabs(eta) > 1.566 && fabs(eta) < 2.5 && mva > -0.58)) {
			  llphoton_cuts_[8] = 1.0;
			  ph.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0.0);
			  llg = dilep + ph;
			  if(llg.M() > 100 && llg.M() < 180)
			    if(llg.M() + mass >= 185)
			      if(ph.Pt()/llg.M() >= 15./110)
				llphoton_cuts_[9] = 1.0;
			}
		      }
		    }
		  }		
		}	      
	      }
	    }
	  }
	}
      }
    }
    else if(b.nmu() >= 2) {
      if (b.mu_pt()->at(0) > 20 && b.mu_pt()->at(1) > 10) {
	llphoton_cuts_[0] = 1.0;
	if (b.mu_charge()->at(0)*b.mu_charge()->at(1) == -1) {
	  llphoton_cuts_[1] = 1.0;
	  if (b.mu_dxy()->at(0) < 0.5 && b.mu_dz()->at(0) < 1.0 && b.mu_dxy()->at(1) < 0.5 && b.mu_dz()->at(1) < 1.0) {
	    llphoton_cuts_[2] = 1.0;
	    if (b.mu_sig()->at(0) && b.mu_sig()->at(1)) {
	      llphoton_cuts_[3] = 1.0;
	      l1.SetPtEtaPhiM(b.mu_pt()->at(0), b.mu_eta()->at(0), b.mu_phi()->at(0), 0.105);
	      l2.SetPtEtaPhiM(b.mu_pt()->at(1), b.mu_eta()->at(1), b.mu_phi()->at(1), 0.105);
	      dilep = l1 + l2;
	      mass = dilep.M();
              if (mass > 50) {
                llphoton_cuts_[4] = 1.0;
		if (b.photon_pt()->size() > 0) {
		  if (b.photon_pt()->at(0) > 15) {
		    llphoton_cuts_[5] = 1.0;
		    if (b.photon_elveto()->at(0)) {
		      llphoton_cuts_[6] = 1.0;
		      if (b.photon_drmin()->at(0) > 0.4) {
			llphoton_cuts_[7] = 1.0;
			double eta = b.photon_eta()->at(0);
			double mva = b.photon_idmva()->at(0);
			if ((fabs(eta) < 1.4442 && mva > -0.4) || (fabs(eta) > 1.566 && fabs(eta) < 2.5 && mva > -0.58)) {
			  llphoton_cuts_[8] = 1.0;
			  ph.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0.0);
			  llg = dilep + ph;
			  if(llg.M() > 100 && llg.M() < 180)
			    if(llg.M() + mass >= 185)
			      if(ph.Pt()/llg.M() >= 15./110)
				llphoton_cuts_[9] = 1.0;
			}
		      }
		    }
		  }		
		}	      
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
    // return 1;
  });

  NamedFunc el_trigs("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  NamedFunc mu_trigs("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
  auto proc_el_sig = Process::MakeShared<Baby_pico>("Electron Channel",  Process::Type::signal, kMagenta-4, {foldersig+"*.root"}, (sig_lepid == 1 || sig_lepid == 2));
  auto proc_mu_sig = Process::MakeShared<Baby_pico>("Muon Channel",  Process::Type::signal, kMagenta-4, {foldersig+"*.root"}, (sig_lepid == 1 || sig_lepid == 3));

  vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig};

  // (el_trigs || mu_trigs) && 
  vector<NamedFunc> cutflow;
  cutflow.push_back("1");
  // cutflow.push_back(el_trigs || mu_trigs);
  cutflow.push_back(llphoton_cuts[0.] > 0.0);
  cutflow.push_back(llphoton_cuts[1] > 0.0);
  cutflow.push_back(llphoton_cuts[2] > 0.0);
  cutflow.push_back(llphoton_cuts[3] > 0.0);
  cutflow.push_back(llphoton_cuts[4] > 0.0);
  cutflow.push_back(llphoton_cuts[5] > 0.0);
  cutflow.push_back(llphoton_cuts[6] > 0.0);
  cutflow.push_back(llphoton_cuts[7] > 0.0);
  cutflow.push_back(llphoton_cuts[8] > 0.0);
  cutflow.push_back(llphoton_cuts[9] > 0.0);

  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::overflow)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> ops = {log_lumi, lin_lumi};
  PlotMaker pm;
  pm.Push<Table>("AN_hgg_cutflow", vector<TableRow>{
      TableRow("Total number of events",                cutflow.at(0),0,0,wgt),
      // TableRow("High level trigger",                    cutflow.at(0),0,0,wgt),
      TableRow("Lepton transverse momentum",            cutflow.at(0),0,0,wgt),
      TableRow("Lepton opposite charge",                cutflow.at(1),0,0,wgt),
      TableRow("Lepton impact parameters",              cutflow.at(2),0,0,wgt),
      TableRow("Lepton identification",                 cutflow.at(3),0,0,wgt),
      TableRow("$m_{ll} > $ 50 GeV",                    cutflow.at(4),0,0,wgt),
      TableRow("Photon transverse momentum",            cutflow.at(5),0,0,wgt),
      TableRow("Photon electron veto",                  cutflow.at(6),0,0,wgt),
      TableRow("Photon delta R minimum",                cutflow.at(7),0,0,wgt),
      TableRow("Photon MVA identification",             cutflow.at(8),0,0,wgt),
      TableRow("Three body invariant mass related cut", cutflow.at(9),0,0,wgt),
      },samples,false);
  pm.min_print_ = true;
  pm.MakePlots(41.5);
}
