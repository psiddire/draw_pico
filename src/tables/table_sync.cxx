#include <iostream>
#include <string>
#include <vector>
#include <memory>

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

int main(){
  gErrorIgnoreLevel = 6000;
  string bfolder("/net/cms29/cms29r0/");
  string foldersig(bfolder+"/pico/NanoAODv5/zgamma_channelIslandsv2/2016/HToZG/unskimmed/*");
  string foldermc(bfolder+"/pico/NanoAODv5/zgamma_channelIslandsv2/2016/mc/unskimmed/*DYJetsToLL*");

  NamedFunc sig_lepid("signal lepton ID",[](const Baby &b) -> NamedFunc::ScalarType{
//     int iz(-1), lepid(0);
    int lepid(0);
//     bool checktau(false);
    for(size_t imc(0); imc < b.mc_id()->size(); imc++) {
      lepid = abs(b.mc_id()->at(imc));
      if(lepid == 11 || lepid == 13 || lepid == 15)
        return lepid;
      if(b.mc_mom()->at(imc) == 23 && b.mc_mom()->at(b.mc_momidx()->at(imc)) == 25) {
        lepid = abs(b.mc_id()->at(imc));
        if(lepid == 11 || lepid == 13 || lepid == 15)
          return lepid;
      }
    }
//       if(b.mc_id()->at(imc) == 23 && b.mc_mom()->at(imc) == 25)
//         iz = imc;
//       if(b.mc_momidx()->at(imc) == iz && b.mc_id()->at(imc) == 15)
//         checktau = true;
//       if(checktau) 
//         if(abs(b.mc_mom()->at(imc)) == 15 && (abs(b.mc_id()->at(imc)) == 11 || abs(b.mc_id()->at(imc)) == 13)) {
//           if(lepid + b.mc_id()->at(imc) == 0)
//             return abs(lepid);
//           lepid = b.mc_id()->at(imc);
//         }
//   if(b.nmu() > 1 && b.nel() > 1) 
//     return b.ll_lepid()->at(0);
//   if(b.nmu() > 1) return 13;
//   if(b.nel() > 1) return 11;
  return lepid;
  });

  NamedFunc nels("Number of electrons with WP98",[](const Baby &b) -> NamedFunc::ScalarType{ 
    int nel(0);
    for(size_t iel(0); iel < b.el_pt()->size(); iel++) {
      double wp[2][3] = {{-0.145237, -0.0315746, -0.032173},
                         { 0.604775,  0.628743,  0.896462 }};
      double pt(b.el_pt()->at(iel)), eta(abs(b.el_eta()->at(iel)));
      int ipt(1), ieta(-1);
      if(pt>10) ipt = 0;
      if(eta < 0.8) ieta = 0;
      else if(eta < 1.479) ieta = 1;
      else ieta = 2;
      if(b.el_idmva()->at(iel) > wp[ipt][ieta]) nel++;
    }
    return nel;
  });
;
  NamedFunc dilep_mass("Dilepton mass",[](const Baby &b) -> NamedFunc::ScalarType{
    double mass(-1);
    for(int ill(0); ill < b.nll(); ill++) 
      if(mass == -1 || abs(b.ll_m()->at(ill) - 91.1876) < abs(mass - 91.1876))
        mass = b.ll_m()->at(ill);
    return mass;
  });

  NamedFunc leading_lep_pt("Leading lepton pt cut",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.ll_pt()->size() == 0) return false;
    int best_ill(-1);
    double mass(-1);
    for(int ill(0); ill < b.nll(); ill++) 
      if(mass == -1 || abs(b.ll_m()->at(ill) - 91.1876) < abs(mass - 91.1876)){
        mass = b.ll_m()->at(ill);
        best_ill = ill;
      }
    if(b.ll_lepid()->at(best_ill) == 11){
      if(b.el_pt()->at(b.ll_i1()->at(best_ill)) > 25)
        return true;
    }
    else if(b.ll_lepid()->at(best_ill) == 13){
      if(b.mu_pt()->at(b.ll_i1()->at(best_ill)) > 20)
        return true;
    }
    return false;
  });
  NamedFunc trailing_lep_pt("Trailing lepton pt cut",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.ll_pt()->size() == 0) return false;
    int best_ill(-1);
    double mass(-1);
    for(int ill(0); ill < b.nll(); ill++) 
      if(mass == -1 || abs(b.ll_m()->at(ill) - 91.1876) < abs(mass - 91.1876)){
        mass = b.ll_m()->at(ill);
        best_ill = ill;
      }
    if(b.ll_lepid()->at(best_ill) == 11){
      if(b.el_pt()->at(b.ll_i2()->at(best_ill)) > 15)
        return true;
    }
    else if(b.ll_lepid()->at(best_ill) == 13){
      if(b.mu_pt()->at(b.ll_i2()->at(best_ill)) > 10)
        return true;
    }
    return false;
  });
  NamedFunc l1emtf("L1EMTF cut for muons",[](const Baby &b) -> NamedFunc::ScalarType{
    if(b.ll_pt()->size() == 0) return true;
    for(size_t ill(0); ill < b.ll_pt()->size(); ill++)
      if(b.ll_m()->at(ill) > 50 && b.ll_lepid()->at(ill) == 13)
        if(b.ll_dphi()->at(ill) < 1.222)
          if(abs(b.mu_eta()->at(b.ll_i1()->at(ill))) > 1.2 &&
             abs(b.mu_eta()->at(b.ll_i2()->at(ill))) > 1.2 &&
             b.ll_deta()->at(ill) < 1.2)
             return false;
    return true;
  });
  NamedFunc photon_pre("Photon preselection",[](const Baby &b) -> NamedFunc::ScalarType{
    for(size_t iph(0); iph < b.photon_pt()->size(); iph++)
      if(b.photon_pt()->at(iph) > 15 && abs(b.photon_eta()->at(iph)) < 2.5 &&
         b.photon_id()->at(iph) && b.photon_elveto()->at(iph))
        return true;
    return false;
  });
  NamedFunc valid_photon("Valid Photon",[](const Baby &b) -> NamedFunc::ScalarType{
    int best_ill(-1);
    double mass(-1);
    for(int ill(0); ill < b.nll(); ill++) 
      if(mass == -1 || abs(b.ll_m()->at(ill) - 91.1876) < abs(mass - 91.1876)){
        mass = b.ll_m()->at(ill);
        best_ill = ill;
      }
    int best_iph(-1);
    for(size_t iph(0); iph < b.photon_pt()->size(); iph++)
      if(b.photon_pt()->at(iph) > 15 && abs(b.photon_eta()->at(iph)) < 2.5 &&
         b.photon_id()->at(iph) && b.photon_elveto()->at(iph)) {
          best_iph = iph;
          break;
        }
    for(size_t illg(0); illg < b.llphoton_pt()->size(); illg++)
      if(b.llphoton_ill()->at(illg) == best_ill && b.llphoton_iph()->at(illg) == best_iph)
        if(b.photon_pt()->at(best_iph)/b.llphoton_m()->at(illg) > 15./110)
          if(b.ll_m()->at(best_ill) + b.llphoton_m()->at(illg) > 185)
            if(b.llphoton_m()->at(illg) > 100 && b.llphoton_m()->at(illg) < 180)
              if(b.photon_drmin()->at(best_iph) > 0.4)
                return true;
    return false;
  });
  NamedFunc stitch_dy("stitch_dy");

  auto proc_el_sig        = Process::MakeShared<Baby_pico>("Sig e",           Process::Type::signal, kMagenta-4, {foldersig+"*.root"}, sig_lepid == 11);
  auto proc_mu_sig        = Process::MakeShared<Baby_pico>("Sig #mu",         Process::Type::signal, kAzure,     {foldersig+"*.root"}, sig_lepid == 13);
  auto proc_tau_sig       = Process::MakeShared<Baby_pico>("Sig #tau",        Process::Type::signal, kAzure,     {foldersig+"*.root"}, sig_lepid == 15);
  auto proc_tot_sig       = Process::MakeShared<Baby_pico>("Total Sig",       Process::Type::signal, kViolet,    {foldersig+"*.root"});
  auto proc_el_dy         = Process::MakeShared<Baby_pico>("DY e",            Process::Type::signal, kMagenta-4, {foldermc +"*.root"}, sig_lepid == 11);
  auto proc_mu_dy         = Process::MakeShared<Baby_pico>("DY #mu",          Process::Type::signal, kAzure,     {foldermc +"*.root"}, sig_lepid == 13);
  auto proc_el_dy_stitch  = Process::MakeShared<Baby_pico>("Stitched DY e",   Process::Type::signal, kMagenta-4, {foldermc +"*.root"}, sig_lepid == 11 && stitch_dy);
  auto proc_mu_dy_stitch  = Process::MakeShared<Baby_pico>("Stitched DY #mu", Process::Type::signal, kAzure,     {foldermc +"*.root"}, sig_lepid == 13 && stitch_dy);

// 	vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig, proc_el_dy, proc_mu_dy, proc_el_dy_stitch, proc_mu_dy_stitch};
	vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig};
// 	vector<shared_ptr<Process> > samples  = {proc_el_sig, proc_mu_sig, proc_tau_sig, proc_tot_sig};
// 	vector<shared_ptr<Process> > samples  = {proc_tot_sig};

	vector<NamedFunc> cuts = {"1",
                            "npv_good>0",
                            (sig_lepid == 11 && nels>1) || (sig_lepid == 13 && "nmu>1"),
                            dilep_mass > 50,
                            leading_lep_pt,
                            trailing_lep_pt,
                            l1emtf,
                            photon_pre,
                            valid_photon};
//                             llphoton_cuts};
  vector<string> cut_labels = {"Total Events",
                               "Good primary vertex",
                               "At least 2 leptons passing ID+preselection(el)",
                               "Valid Z pairing",
                               "Leading lepton $p_{T}$ cut",
                               "Trailing lepton $p_{T}$ cut",
                               "L1EMTF cut (#mu only)",
                               "At least 1 photon passing ID+preselection",
                               "Valid photon"};
	vector<NamedFunc> cutflow;
  NamedFunc wgt("w_lumi*w_lep*w_photon");
	NamedFunc cut("1");
	for(size_t i = 0; i < cuts.size(); i++) {
	  cut = cut && cuts.at(i);
		cutflow.push_back(cut);
	}
  PlotOpt log_lumi("txt/plot_styles.txt","CMSPaper");
  log_lumi.Title(TitleType::info)
          .YAxis(YAxisType::log)
          .Stack(StackType::shapes)
          .Overflow(OverflowType::overflow)
          .FileExtensions({"pdf"});
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  vector<PlotOpt> ops = {log_lumi, lin_lumi};
  PlotMaker pm;
//   pm.Push<Hist1D>(Axis(2,0,2, !(sig_lepid == 11 || sig_lepid == 13), "ID != lep", {}),"1", samples,ops);
//   pm.Push<Hist1D>(Axis(30,-0.5,29.5, sig_lepid, "sig_lepid", {}),"1", samples,ops).Weight(wgt);
//   pm.Push<Hist1D>(Axis(30,-0.5,29.5, sig_lepid, "sig_lepid", {}),"1", samples,ops);
  pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nel", "Nel", {}),!(sig_lepid == 11 || sig_lepid == 13 || sig_lepid == 15), samples,ops).Weight(wgt);
  pm.Push<Hist1D>(Axis(5,-0.5,4.5, "nmu", "Nmu", {}),!(sig_lepid == 11 || sig_lepid == 13 || sig_lepid == 15), samples,ops).Weight(wgt);
//   pm.Push<Hist1D>(Axis(5,-0.5,4.5, nels, "Size of electron list", {}),"1", samples,ops).Weight(wgt);
//   pm.Push<Hist1D>(Axis(60,0,120, dilep_mass, "dilep (ee) mass", {}),"1", samples,ops).Weight(wgt);
  vector<TableRow> cutrows;
  vector<TableRow> cutrows_weight;
  for(size_t ic(0); ic < cuts.size(); ic++) {
    cutrows.push_back(TableRow(cut_labels.at(ic),cutflow.at(ic),0,0,wgt));
    cutrows_weight.push_back(TableRow(cut_labels.at(ic),cutflow.at(ic)));
    }
  pm.Push<Table>("full_sync_cutflow_wlumi",      cutrows,        samples,false);
  pm.Push<Table>("full_eff_sync_cutflow_wlumi",  cutrows,        samples,false,true,false,false,true);
  pm.Push<Table>("full_sync_cutflow_weight",     cutrows_weight, samples,false);
  pm.Push<Table>("full_eff_sync_cutflow_weight", cutrows_weight, samples,false,true,false,false,true);
  pm.min_print_ = true;
  pm.MakePlots(35.9);
}
