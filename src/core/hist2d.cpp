#include "core/hist2d.hpp"

#include <string>
#include <vector>
#include <sstream>

#include <sys/stat.h>

#include "TStyle.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"
#include "TArrow.h"
#include "core/named_func.hpp"

using namespace std;
using namespace PlotOptTypes;

Hist2D::SingleHist2D::SingleHist2D(const Hist2D &figure,
                                   const std::shared_ptr<Process> &process,
                                   const TH2D &hist_template):
  FigureComponent(figure, process),
  clusterizer_(hist_template, 10000),
  proc_and_hist_cut_(figure.cut_ && process->cut_),
  cut_vector_(),
  wgt_vector_(),
  xval_vector_(),
  yval_vector_(){
}

void Hist2D::SingleHist2D::RecordEvent(const Baby &baby){
  const Hist2D& hist = static_cast<const Hist2D&>(figure_);
  size_t min_vec_size;
  bool have_vec = false;

  const NamedFunc &cut = proc_and_hist_cut_;
  if(cut.IsScalar()){
    if(!cut.GetScalar(baby)) return;
  }else{
    cut_vector_ = cut.GetVector(baby);
    if(!HavePass(cut_vector_)) return;
    have_vec = true;
    min_vec_size = cut_vector_.size();
  }

  const NamedFunc &wgt = hist.weight_;
  NamedFunc::ScalarType wgt_scalar = 0.;
  if(wgt.IsScalar()){
    wgt_scalar = wgt.GetScalar(baby);
  }else{
    wgt_vector_ = wgt.GetVector(baby);
    if(!have_vec || wgt_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = wgt_vector_.size();
    }
  }

  const NamedFunc &xval = hist.xaxis_.var_;
  NamedFunc::ScalarType xval_scalar = 0.;
  if(xval.IsScalar()){
    xval_scalar = xval.GetScalar(baby);
  }else{
    xval_vector_ = xval.GetVector(baby);
    if(!have_vec || xval_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = xval_vector_.size();
    }
  }

  const NamedFunc &yval = hist.yaxis_.var_;
  NamedFunc::ScalarType yval_scalar = 0.;
  if(yval.IsScalar()){
    yval_scalar = yval.GetScalar(baby);
  }else{
    yval_vector_ = yval.GetVector(baby);
    if(!have_vec || yval_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = yval_vector_.size();
    }
  }

  if(!have_vec){
    clusterizer_.AddPoint(xval_scalar, yval_scalar, wgt_scalar);
  }else{
    for(size_t i = 0; i < min_vec_size; ++i){
      if(cut.IsVector() && !cut_vector_.at(i)) continue;
      clusterizer_.AddPoint(xval.IsScalar() ? xval_scalar : xval_vector_.at(i),
                            yval.IsScalar() ? yval_scalar : yval_vector_.at(i),
                            wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(i));
    }
  }
}

Hist2D::Hist2D(const Axis &xaxis, const Axis &yaxis, const NamedFunc &cut,
               const std::vector<std::shared_ptr<Process> > &processes,
               const std::vector<PlotOpt> &plot_options):
  Figure(),
  xaxis_(xaxis),
  yaxis_(yaxis),
  cut_(cut),
  weight_("weight"),
  tag_(""),
  plot_options_(plot_options),
  backgrounds_(),
  signals_(),
  datas_(),
  this_opt_(PlotOpt()),
  luminosity_(){
  string xtitle = xaxis_.title_;
  if(xaxis.units_ != "") xtitle += " ["+xaxis.units_+"]";
  string ytitle = yaxis_.title_;
  if(yaxis.units_ != "") ytitle += " ["+yaxis.units_+"]";

  TH2D empty("", (";"+xtitle+";"+ytitle).c_str(),
             xaxis_.Nbins(), &xaxis_.Bins().at(0),
             yaxis_.Nbins(), &yaxis_.Bins().at(0));
  empty.SetStats(false);
  empty.Sumw2(true);
  for(const auto &process: processes){
    TH2D hist_template = empty;
    hist_template.SetFillColor(process->GetFillColor());
    hist_template.SetFillStyle(process->GetFillStyle());
    hist_template.SetLineColor(process->GetLineColor());
    hist_template.SetLineStyle(process->GetLineStyle());
    hist_template.SetLineWidth(process->GetLineWidth());
    hist_template.SetMarkerColor(process->GetMarkerColor());
    hist_template.SetMarkerStyle(process->GetMarkerStyle());
    hist_template.SetMarkerSize(process->GetMarkerSize());
    unique_ptr<SingleHist2D> hist(new SingleHist2D(*this, process, hist_template));

    switch(process->type_){
    case Process::Type::data:
      datas_.push_back(move(hist));
      break;
    case Process::Type::background:
      backgrounds_.push_back(move(hist));
      break;
    case Process::Type::signal:
      signals_.push_back(move(hist));
      break;
    default:
      break;
    }
  }
}

void Hist2D::Print(double luminosity,
                   const string &subdir){
  luminosity_ = luminosity;
  for(const auto &opt: plot_options_){
    this_opt_ = opt;
    this_opt_.MakeSane();
    gStyle->SetColorModelPS(this_opt_.UseCMYK());
    gROOT->ForceStyle();

    MakeOnePlot(subdir);
  }
}

void Hist2D::MakeOnePlot(const string &subdir){
  bool bkg_is_hist;
  switch(this_opt_.Stack()){
  default:
  case StackType::signal_overlay:
  case StackType::signal_on_top:
  case StackType::data_norm:
    bkg_is_hist = true;
    break;
  case StackType::lumi_shapes:
  case StackType::shapes:
    bkg_is_hist = false;
    break;
  }

  const Int_t NRGBs = 5;
  const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs] = { 0.71, 0.50, 1.00, 1.00, 1.00 };
  Double_t green[NRGBs] = { 0.80, 1.00, 1.00, 0.60, 0.50 };
  Double_t blue[NRGBs] = { 0.95, 1.00, 0.50, 0.40, 0.50 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  TCanvas c("canvas","canvas", this_opt_.CanvasWidth(), this_opt_.CanvasHeight());
  c.cd();
  c.SetTicks(1,1);
  c.SetFillStyle(4000);
  c.SetMargin(this_opt_.LeftMargin(),
	      bkg_is_hist ? this_opt_.RightMargin() : 0.05,
	      this_opt_.BottomMargin(),
	      this_opt_.TopMargin());
  if (this_opt_.YAxis() == YAxisType::log)
    c.SetLogz();

  vector<TLine> lines = GetLines();
  vector<shared_ptr<TLatex> > labels = GetLabels(bkg_is_hist);
  TLegend legend(bkg_is_hist ? this_opt_.LeftMargin() : this_opt_.LeftMargin()-0.05, 1.-this_opt_.TopMargin(),
                 bkg_is_hist ? 1.-this_opt_.RightMargin() : 1.-0.05, 1.);
  legend.SetNColumns(datas_.size() + signals_.size() + (bkg_is_hist ? 0 : backgrounds_.size()));
  legend.SetBorderSize(0);
  legend.SetTextSize(0.05);
  legend.SetFillStyle(4000);

  TH2D bkg_hist = GetBkgHist(bkg_is_hist);
  vector<TGraph> bkg_graphs;
  if(bkg_is_hist) bkg_graphs.clear();
  else bkg_graphs = GetGraphs(backgrounds_, true);

  vector<TGraph> sig_graphs = GetGraphs(signals_, true);
  vector<TGraph> data_graphs = GetGraphs(datas_, false);

  for(size_t i = 0; i < datas_.size(); ++i){
    AddEntry(legend, *datas_.at(i), data_graphs.at(i));
  }
  if(!bkg_is_hist){
    for(size_t i = 0; i < backgrounds_.size(); ++i){
      AddEntry(legend, *backgrounds_.at(i), bkg_graphs.at(i));
    }
  }
  for(size_t i = 0; i < signals_.size(); ++i){
    AddEntry(legend, *signals_.at(i), sig_graphs.at(i));
  }

  bkg_hist.Draw("axis");
  if(bkg_is_hist){
    bkg_hist.Draw("COLZ"); // "BOX same"
  }else{
    for(auto &g: bkg_graphs){
      g.Draw("COLZ"); // "p"
    }
  }
  for(auto &l: lines){
    l.Draw();
  }
  for(auto &g: data_graphs){
    g.Draw("COLZ"); // "p"
  }
  for(auto &g: sig_graphs){
    g.Draw("COLZ"); // "p"
  }
  // for(auto &l: labels){
  //   l->Draw("same");
  // }
  legend.Draw("same");
  bkg_hist.Draw("axis same");

    double height = 0.125;
    double width = 0.125;
    TPaveText l1(this_opt_.LeftMargin()+0.16,  -0.014,
		 this_opt_.LeftMargin()+width+0.16, -0.014+height, "NDCNB");
    TPaveText l2(1.- this_opt_.RightMargin()+0.1-width,  this_opt_.BottomMargin(),
		 1.- this_opt_.RightMargin()+0.1,  this_opt_.BottomMargin()+height, "NDCNB");
    TPaveText l3(this_opt_.LeftMargin()+0.16, 1.- this_opt_.TopMargin()-2.0*height,
		 this_opt_.LeftMargin()+0.16+width, 1.- this_opt_.TopMargin()-height, "NDCNB");
    TPaveText l4(1.- this_opt_.RightMargin()+0.1-width, 1.- this_opt_.TopMargin()-2.0*height,
		 1.- this_opt_.RightMargin()+0.1, 1.- this_opt_.TopMargin()-height, "NDCNB");

    TArrow arrow;
    arrow.SetLineColor(kGray+2); arrow.SetFillColor(kGray+2);
    arrow.SetArrowSize(0.015); arrow.SetLineWidth(4);

    l1.AddText("R1");
    l2.AddText("R2");
    l3.AddText("R3");
    l4.AddText("R4");

    SetStyle(l1);
    SetStyle(l2);
    SetStyle(l3);
    SetStyle(l4);
    if(!bkg_is_hist){
      l1.Draw("same");
      l2.Draw("same");
      l3.Draw("same");
      l4.Draw("same");
      arrow.DrawArrow(325, -29, 325, -10);
   }


  if(subdir != "") mkdir(("plots/"+subdir).c_str(), 0777);
  string base_name = subdir != ""
    ? "plots/"+subdir+"/"+Name()
    : "plots/"+Name();
  for(const auto &ext: this_opt_.FileExtensions()){
    string full_name = base_name+"__"+this_opt_.TypeString()+'.'+ext;
    if (Contains(tag_,"FixName:")){
      string tagName=tag_;
      ReplaceAll(tagName, "FixName:", "");
      base_name = subdir != ""
        ? "plots/"+subdir+"/"+tagName
        : "plots/"+tagName;
      full_name = base_name+'.'+ext;
    }
    c.Print(full_name.c_str());
    cout << "open " << full_name << endl;
  }
}

TH2D Hist2D::GetBkgHist(bool bkg_is_hist) const{
  string units = (xaxis_.units_ == yaxis_.units_) ? (xaxis_.units_+"^{2}") : (xaxis_.units_+"*"+yaxis_.units_);
  string z_title = bkg_is_hist
    ? ("Simulated Events/(" + ToString(xaxis_.AvgBinWidth()*yaxis_.AvgBinWidth())+" "+units+")")
    : "";
  string title = ";"+xaxis_.Title()+";"+yaxis_.Title()+";"+z_title;
  TH2D h;
  if(bkg_is_hist && backgrounds_.size() > 0){
    h = backgrounds_.front()->clusterizer_.GetHistogram(luminosity_);
    for(size_t i = 1; i < backgrounds_.size(); ++i){
      TH2D to_add = backgrounds_.at(i)->clusterizer_.GetHistogram(luminosity_);
      h.Add(&to_add);
    }
  }else{
    h = TH2D("", title.c_str(),
             xaxis_.Nbins(), &xaxis_.Bins().at(0),
             yaxis_.Nbins(), &yaxis_.Bins().at(0));
  }
  h.SetTitle(title.c_str());
  h.GetZaxis()->SetTitle(z_title.c_str());
  h.SetStats(0);
  h.SetMinimum(this_opt_.LogMinimum());
  h.SetLabelOffset(0.011);
  h.SetTitleOffset(this_opt_.XTitleOffset(), "x");
  h.SetTitleOffset(this_opt_.YTitleOffset(), "y");
  h.SetTitleOffset(this_opt_.ZTitleOffset(), "z");
  h.SetLabelSize(this_opt_.LabelSize(), "xyz");
  h.SetTitleSize(this_opt_.TitleSize(), "xyz");
  return h;
}

TH2D Hist2D::GetDataHist() const{
  string units = (xaxis_.units_ == yaxis_.units_) ? (xaxis_.units_+"^{2}") : (xaxis_.units_+"*"+yaxis_.units_);
  string z_title = ("Data Events/(" + ToString(xaxis_.AvgBinWidth()*yaxis_.AvgBinWidth())+" "+units+")");
  string title = ";"+xaxis_.Title()+";"+yaxis_.Title()+";"+z_title;
  TH2D h;
  //change this to signal
  if(datas_.size() > 0){
    h = datas_.front()->clusterizer_.GetHistogram();
    for(size_t i = 1; i < datas_.size(); ++i){
      TH2D to_add = datas_.at(i)->clusterizer_.GetHistogram();
      h.Add(&to_add);
    }
  }else{
    h = TH2D("", title.c_str(),
             xaxis_.Nbins(), &xaxis_.Bins().at(0),
             yaxis_.Nbins(), &yaxis_.Bins().at(0));
  }
  h.SetTitle(title.c_str());
  h.GetZaxis()->SetTitle(z_title.c_str());
  h.SetStats(0);
  h.SetMinimum(this_opt_.LogMinimum());
  h.SetLabelOffset(0.011);
  h.SetTitleOffset(this_opt_.XTitleOffset(), "x");
  h.SetTitleOffset(this_opt_.YTitleOffset(), "y");
  h.SetTitleOffset(this_opt_.ZTitleOffset(), "z");
  h.SetLabelSize(this_opt_.LabelSize(), "xyz");
  h.SetTitleSize(this_opt_.TitleSize(), "xyz");
  return h;
}

vector<TGraph> Hist2D::GetGraphs(const vector<unique_ptr<SingleHist2D> > &components,
				 bool lumi_weighted) const{
  vector<TGraph> graphs(components.size());
  for(size_t i = 0; i < components.size(); ++i){
    graphs.at(i) = components.at(i)->clusterizer_.GetGraph(lumi_weighted ? luminosity_ : 1.);
    graphs.at(i).SetName(components.at(i)->process_->name_.c_str());
  }
  return graphs;
}

vector<TLine> Hist2D::GetLines() const{
  vector<TLine> lines;
  for(const auto &v: xaxis_.cut_vals_){
    lines.emplace_back(v, yaxis_.Bins().front(), v, yaxis_.Bins().back());
  }
  for(const auto &v: yaxis_.cut_vals_){
    lines.emplace_back(xaxis_.Bins().front(), v, xaxis_.Bins().back(), v);
  }
  for(auto &l: lines){
    l.SetLineStyle(2);
    l.SetLineWidth(4);
    l.SetLineColor(kBlack);
  }
  return lines;
}

vector<shared_ptr<TLatex> > Hist2D::GetLabels(bool bkg_is_hist) const{
  float left = this_opt_.LeftMargin()+0.03;
  float right = bkg_is_hist ? 1.-this_opt_.RightMargin()-0.03 : 1.-0.05-0.03;
  float top = 1.-this_opt_.TopMargin()-0.03;
  vector<shared_ptr<TLatex> > labels;
  string extra;
  switch(this_opt_.Title()){
  case TitleType::preliminary: extra = "Preliminary"; break;
  case TitleType::simulation: extra = "Simulation"; break;
  case TitleType::simulation_preliminary: extra = "#splitline{Simulation}{Preliminary}"; break;
  case TitleType::simulation_supplementary: extra = "#splitline{Simulation}{Supplementary}"; break;
  case TitleType::supplementary: extra = "Supplementary"; break;
  case TitleType::data: extra = ""; break;
  case TitleType::info: extra = ""; break;
  default:
    ERROR("Did not understand title type "+to_string(static_cast<int>(this_opt_.Title())));
  }
  labels.push_back(make_shared<TLatex>(left-0.013, top+0.005, "#font[62]{CMS}"));
  labels.back()->SetNDC();
  labels.back()->SetTextAlign(13);
  labels.back()->SetTextFont(this_opt_.Font());
  if(extra.find("splitline")!=std::string::npos)
    labels.push_back(make_shared<TLatex>(left-0.022, top-0.052,
					 ("#scale[0.69]{#font[52]{ "+extra+"}}").c_str()));

  else labels.push_back(make_shared<TLatex>(left-0.023, top-0.05,
					    ("#scale[0.67]{#font[52]{ "+extra+"}}").c_str()));
  labels.back()->SetNDC();
  labels.back()->SetTextAlign(13);
  labels.back()->SetTextFont(this_opt_.Font());

  if(this_opt_.Title() == TitleType::simulation_supplementary || this_opt_.Title() == TitleType::supplementary ){
    labels.push_back(make_shared<TLatex>(right, top-0.062,
					 "#scale[0.73]{#font[82]{arXiv:xxxx.xxxxx}}"));
    labels.back()->SetNDC();
    labels.back()->SetTextAlign(33);
    labels.back()->SetTextFont(this_opt_.Font());
  }

  ostringstream oss;
  if (luminosity_tag_ != "") oss << luminosity_tag_ << " fb^{-1} (13 TeV)" << flush;
  else if (luminosity_<1.1) oss << "137 fb^{-1} (13 TeV)" << flush;
  else oss << luminosity_ << " fb^{-1} (13 TeV)" << flush;
  labels.push_back(make_shared<TLatex>(right+0.01, top, oss.str().c_str()));
  labels.back()->SetNDC();
  labels.back()->SetTextAlign(33);
  labels.back()->SetTextFont(this_opt_.Font());
  labels.back()->SetTextSize(0.043);

  return labels;
}

string Hist2D::Name() const{
  string cut = "";
  if(cut_.Name() != "1") cut = "__"+cut_.Name();

  string weight = "";
  if(weight_.Name() != "weight") weight = "__"+weight_.Name();

  if(tag_ == ""){
    return CodeToPlainText(yaxis_.var_.Name()+"__"+xaxis_.var_.Name()+cut+weight);
  }else if (Contains(tag_,"ShortName:")){
    string tagName=tag_;
    ReplaceAll(tagName, "ShortName:", "");
    // return CodeToPlainText(tagName+"__"+xaxis_.var_.Name()+weight);
    return CodeToPlainText(tagName);
  }else{
    if (Contains(tag_,"/")) {
      string token = tag_.substr(0, tag_.find("/"));
      mkdir(("plots/"+token).c_str(), 0777);
    }
    return tag_;
  }
}

Hist2D & Hist2D::Weight(const NamedFunc &weight){
  weight_ = weight;
  return *this;
}

Hist2D & Hist2D::Tag(const std::string &tag){
  tag_ = tag;
  return *this;
}

Hist2D & Hist2D::LuminosityTag(const string &tag){
  luminosity_tag_ = tag;
  return *this;
}

void Hist2D::AddEntry(TLegend &l, const SingleHist2D &h, const TGraph &g) const{
  string name = h.process_->name_;
  ostringstream oss;
  oss << name;
  bool print_rho;
  switch(this_opt_.Title()){
  case TitleType::info:
    print_rho = true;
    break;
  case TitleType::preliminary:
    print_rho = false;
    break;
  case TitleType::simulation:
  case TitleType::simulation_preliminary:
    print_rho = true;
    break;
  case TitleType::simulation_supplementary:
    print_rho = true;
    break;
  case TitleType::supplementary:
  case TitleType::data:
  default:
    print_rho = false;
    break;
  }
  if(h.process_->type_ == Process::Type::signal){
    print_rho = false;
  }

  if(print_rho){
    double rho = h.clusterizer_.GetHistogram(1.).GetCorrelationFactor();
    oss.precision(1);
    oss << " [#rho=" << rho << "]" << flush;
  }
  oss << flush;
  l.AddEntry(&g, oss.str().c_str(), "p");
}

set<const Process*> Hist2D::GetProcesses() const{
  set<const Process*> processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_.get());
  }
  return processes;
}

std::string Hist2D::GetTag() const{
  return tag_;
}

Figure::FigureComponent * Hist2D::GetComponent(const Process *process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

const vector<unique_ptr<Hist2D::SingleHist2D> >& Hist2D::GetComponentList(const Process *process){
  switch(process->type_){
  case Process::Type::data:
    return datas_;
  case Process::Type::background:
    return backgrounds_;
  case Process::Type::signal:
    return signals_;
  default:
    ERROR("Did not understand process type "+to_string(static_cast<long>(process->type_))+".");
    return backgrounds_;
  }
}

void SetStyle(TPaveText &pt){
  pt.SetFillColor(0);
  pt.SetFillStyle(4000);
  pt.SetBorderSize(0);
  pt.SetTextColorAlpha(1,0.5);
}
