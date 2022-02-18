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
#include "core/hist2d.hpp"
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
  string sig_path(bfolder+"NanoAODv2/zgamma_gg/2017/gg/unskimmed/");

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight;
  });

  NamedFunc eephoton_cuts("three body invariant mass cuts",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector l1, l2, dilep, ph, llg;
      if (b.nel() >= 2)
	if (b.el_pt()->at(0) > 25 && b.el_pt()->at(1) > 15)
	  if (b.el_charge()->at(0)*b.el_charge()->at(1) == -1)
	    if (b.el_dxy()->at(0) < 0.5 && b.el_dz()->at(0) < 1.0 && b.el_dxy()->at(1) < 0.5 && b.el_dz()->at(1) < 1.0)
	      if (b.el_sig()->at(0) && b.el_sig()->at(1)) {
		l1.SetPtEtaPhiM(b.el_pt()->at(0), b.el_eta()->at(0), b.el_phi()->at(0), 0.000511);
		l2.SetPtEtaPhiM(b.el_pt()->at(1), b.el_eta()->at(1), b.el_phi()->at(1), 0.000511);
		dilep = l1 + l2;
		if (dilep.M() > 50)
		  if (b.photon_pt()->size() > 0)
		    if (b.photon_pt()->at(0) > 15)
		      if (b.photon_elveto()->at(0))
			// if (b.photon_drmin()->at(0) > 0.4) {
			  // double eta = b.photon_eta()->at(0);
			  // double mva = b.photon_idmva()->at(0);
			  // if ((fabs(eta) < 1.4442 && mva > -0.4) || (fabs(eta) > 1.566 && fabs(eta) < 2.5 && mva > -0.58)) {
			    // ph.SetPtEtaPhiM(b.photon_pt()->at(0), b.photon_eta()->at(0), b.photon_phi()->at(0), 0.0);
			    // llg = dilep + ph;
			    // if(llg.M() > 100 && llg.M() < 180)
			      // if(llg.M() + dilep.M() >= 185)
				// if(ph.Pt()/llg.M() >= 15./110)
				  return dilep.M();
			  // }
			// }
	      }
      return 0;
    });

  NamedFunc trigs("1");

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

  auto proc_UL = Process::MakeShared<Baby_pico>("HTo#gamma#gamma MC sample", sig, TColor::GetColor("#ff0000"), {sig_path+"*.root"}, trigs);
  proc_UL->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_UL};

  PlotMaker pm;
  NamedFunc selection = trigs;
  pm.Push<Hist1D>(Axis(100, 50, 150, eephoton_cuts, "M_{ee} (GeV)", {}), selection, procs, ops).Weight(wgt).Tag("mass_Zee");

  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

