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
  string sig_path(bfolder+"NanoAODv2/zgamma_gg/2017/gg/merged_zgmc_gg/");

  NamedFunc mass_hgg("mass_hgg",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector H = AssignHGG(b);
      return H.M();
    });

  NamedFunc pt_higgs("pt_higgs",[](const Baby &b) -> NamedFunc::ScalarType{
      TLorentzVector H = AssignHGG(b);
      return H.Pt();
    });

  NamedFunc wgt("weight",[](const Baby &b) -> NamedFunc::ScalarType{
    double weight = b.w_lumi();
    return weight;
  });

  NamedFunc trigs("1");
  NamedFunc baseline("nphoton > 1");
  NamedFunc pho("photon_pt[0] > 15 && photon_drmin[0] > 0.4 && photon_sig[0] && photon_pt[1] > 15 && photon_drmin[1] > 0.4 && photon_sig[1]");

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

  auto proc_UL = Process::MakeShared<Baby_pico>("HTo#gamma#gamma", sig, TColor::GetColor("#ff0000"), {sig_path+"merged_pico_gg_GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8_zgmc_gg_nfiles_7.root"}, trigs);
  proc_UL->SetLineWidth(3);
  vector<shared_ptr<Process>> procs = {proc_UL};

  PlotMaker pm;
  NamedFunc selection = baseline && pho;
  pm.Push<Hist1D>(Axis(50, 100, 150, mass_hgg, "M_{#gamma#gamma} (GeV)", {}), selection, procs, ops).Weight(wgt).Tag("mass_higgs");
  pm.Push<Hist1D>(Axis(20, 0, 200, pt_higgs, "Higgs p_{T} (GeV)", {}), selection, procs, ops).Weight(wgt).Tag("pt_higgs");
  pm.Push<Hist1D>(Axis(20, 50, 150, "photon_pt[0]", "Leading photon p_{T} [GeV]", {}), selection, procs, ops).Weight(wgt).Tag("pt_photon");
  pm.Push<Hist1D>(Axis(5, 0, 5, "nphoton", "Number of photons", {}), selection, procs, ops).Weight(wgt).Tag("num_photon");

  pm.min_print_ = true;
  pm.MakePlots(41.5);
}

