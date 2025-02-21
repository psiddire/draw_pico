/*! \class Hist1D

  \brief A full 1D plot with stacked/overlayed histograms

  Hist1D contains all the information necessary to produce a single 1D plot
  containing a combination of background MC, signal MC, and data histograms. The
  content and style of the plot are maintained (mostly) independently, so that
  once the histograms have been filled with data, redrawing with multiple styles
  has minimal overhead.

  To produce a plot, the component histograms in Hist1D::backgrounds_,
  Hist1D::signals_, and Hist1D::datas_ must be filled (e.g., by
  PlotMaker). Once the data is ready, a call to Hist1D::Print() will
  generate the formatted plot for each style contained in
  Hist1D::plot_options_.

*/

/*! \class Hist1D::SingleHist1D

  \brief Container for a TH1D associated with a single Process

  Hist1D::SingleHist1D is mostly a "dumb" container used by Hist1D for
  convenience. It contains a pointer to a single Process and a TH1D.

  Until I have a more elegant solution, it also contains a second TH1D which is
  (ab)used by Hist1D to draw the stacked and luminosity scaled histogram
  without disturbing the data in the main TH1D.
*/

#include "core/hist1d.hpp"

#include <cmath>

#include <algorithm>
#include <sstream>

#include <sys/stat.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TBox.h"
#include "TLegendEntry.h"
#include "TFile.h"

#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  class Counter{
  public:
    Counter():
      count_(0){
    }

    string operator()(){
      return to_string(count_++);
    }
  private:
    unsigned long count_;
  } counter;

  /*!\brief Draws all histograms to current canvas, updating draw_opt to contain
    "same" as needed

    \param[in] hists List of histograms to draw

    \param[in,out] draw_opt Option string used by TH1D to draw
    histogram. Changed to "hist same" after first histogram is drawn
  */
  void DrawAll(const vector<unique_ptr<Hist1D::SingleHist1D> > &hists,
               string &draw_opt, bool reversed = false){
    if(!reversed){
      for(auto &hist: hists){
        hist->scaled_hist_.Draw(draw_opt.c_str());
        if(!Contains(draw_opt, "same")){
          draw_opt = draw_opt + " same";
        }
      }
    }else{
      for(auto h = hists.crbegin(); h != hists.crend(); ++h){
        auto &hist = *h;
        hist->scaled_hist_.Draw(draw_opt.c_str());
        if(!Contains(draw_opt, "same")){
          draw_opt = draw_opt + " same";
        }
      }
    }

    // What can I say...
    // TBox box;
    // box.SetFillColorAlpha(kAzure+10,0.1);
    // box.DrawBox(99.5, 243, 131.5, 251.7);
    // box.DrawBox(99.5, 34.1, 131.5, 35.3);
    // box.DrawBox(99.5, 16.81, 131.5, 17.39);
    // box.Draw();
  }

  /*!\brief Erases x-axis title and labels from plot

    \param[in,out] h Histogram to modify
  */
  void StripXLabels(TH1D &h){
    h.GetXaxis()->SetTitle("");
    h.SetLabelSize(0., "x");
    h.SetTitleSize(0., "x");
  }

  /*!\brief Determines which legend column to put an entry in

    \param[in] entry Number of entries already in legend

    \param[in] n_entries Total entries that will go in legend

    \param[in] n_columns Number of columns in legend

    \return Which column to put entry in. 0 is leftmost, n_columns-1 is
    rightmost.
  */
  size_t GetLegendIndex(size_t entry, size_t n_entries, size_t n_columns){
    size_t entries_per_column = n_entries / n_columns;
    size_t cols_with_extra_entry = n_entries % n_columns;
    size_t this_col = -1;
    size_t this_col_end = 0;
    while(this_col_end <= entry){
      ++this_col;
      this_col_end += entries_per_column;
      if(this_col < cols_with_extra_entry){
        ++this_col_end;
      }
    }
    return this_col;
  }
}

TH1D Hist1D::blank_ = TH1D();

/*!\brief Standard constructor

  \param[in] process Process used to fill histogram

  \param[in] hist A fully styled, and typically unfilled histogram to start from
*/
Hist1D::SingleHist1D::SingleHist1D(const Hist1D &figure,
                                   const std::shared_ptr<Process> &process,
                                   const TH1D &hist):
  FigureComponent(figure, process),
  raw_hist_(hist),
  scaled_hist_(),
  proc_and_hist_cut_(figure.cut_ && process->cut_),
  cut_vector_(),
  wgt_vector_(),
  val_vector_(){
  raw_hist_.Sumw2();
  scaled_hist_.Sumw2();
  raw_hist_.SetBinErrorOption(TH1::kPoisson);
  scaled_hist_.SetBinErrorOption(TH1::kPoisson);
}

void Hist1D::SingleHist1D::RecordEvent(const Baby &baby){
  const Hist1D& stack = static_cast<const Hist1D&>(figure_);
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
  const NamedFunc &wgt = stack.weight_;
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

  const NamedFunc &val = stack.xaxis_.var_;
  NamedFunc::ScalarType val_scalar = 0.;
  if(val.IsScalar()){
    val_scalar = val.GetScalar(baby);
  }else{
    val_vector_ = val.GetVector(baby);
    if(!have_vec || val_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = val_vector_.size();
    }
  }

  if(!have_vec){
    raw_hist_.Fill(val_scalar, wgt_scalar);
  }else{
    for(size_t i = 0; i < min_vec_size; ++i){
      if(cut.IsVector() && !cut_vector_.at(i)) continue;
      raw_hist_.Fill(val.IsScalar() ? val_scalar : val_vector_.at(i),
                     wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(i));
    }
  }
}

/*! Get the maximum of the histogram

  \param[in] max_bound Returns the highest bin content c satisfying
  c<max_bound. Usually infinity.

  \param[in] include_error_bar If true, use bin content+error instead of just
  content

  \param[in] include_overflow If true, also check height of underflow and
  overflow bins

  \return Histogram maximum
*/
double Hist1D::SingleHist1D::GetMax(double max_bound,
                                    bool include_error_bar,
                                    bool include_overflow) const{
  int start_bin = include_overflow ? 0 : 1;
  int end_bin = include_overflow ? (scaled_hist_.GetNbinsX()+1) : scaled_hist_.GetNbinsX();
  double the_max = -numeric_limits<double>::infinity();
  for(int bin = start_bin; bin <= end_bin; ++bin){
    double content = scaled_hist_.GetBinContent(bin);
    if(include_error_bar){
      content += scaled_hist_.GetBinErrorUp(bin);
    }
    if(content > the_max && content < max_bound){
      the_max = content;
    }
  }
  return the_max;
}

/*! Get the minimum of the histogram

  \param[in] min_bound Returns the lowest bin content c satisfying
  c>min_bound. Usually zero.

  \param[in] include_error_bar If true, use bin content-error instead of just
  content

  \param[in] include_overflow If true, also check height of underflow and
  overflow bins

  \return Histogram minimum
*/
double Hist1D::SingleHist1D::GetMin(double min_bound,
                                    bool include_error_bar,
                                    bool include_overflow) const{
  int start_bin = include_overflow ? 0 : 1;
  int end_bin = include_overflow ? scaled_hist_.GetNbinsX() : (scaled_hist_.GetNbinsX()+1);
  double the_min = numeric_limits<double>::infinity();
  for(int bin = start_bin; bin <= end_bin; ++bin){
    double content = scaled_hist_.GetBinContent(bin);
    if(include_error_bar){
      content -= fabs(scaled_hist_.GetBinErrorLow(bin));
    }
    if(content < the_min && content > min_bound){
      the_min = content;
    }
  }
  return the_min;
}

/*! \brief Standard constructor

  \param[in] processes List of process for the component histograms

  \param[in] definition Specification of contents (plotted variable, binning,
  etc.)

  \param[in] plot_options Styles with which to draw plot
*/
Hist1D::Hist1D(const Axis &xaxis, const NamedFunc &cut,
               const std::vector<std::shared_ptr<Process> > &processes,
               const std::vector<PlotOpt> &plot_options):
  Figure(),
  xaxis_(xaxis),
  cut_(cut),
  weight_("weight"),
  tag_(""),
  left_label_({}),
  right_label_({}),
  yaxis_zoom_(1.),
  ratio_numerator_(""),
  ratio_denominator_(""),
  plot_options_(plot_options),
  draw_plot_(true),
  backgrounds_(),
  signals_(),
  datas_(),
  this_opt_(PlotOpt()),
  luminosity_(),
  mc_scale_(),
  mc_scale_error_(){
  if(plot_options_.size() > 0) this_opt_ = plot_options_.front();

  string x_title = xaxis_.title_;
  if(xaxis_.units_ != "") x_title += " ["+xaxis_.units_+"]";

  TH1D empty("", (";"+x_title+";").c_str(), xaxis_.Nbins(), &xaxis_.Bins().at(0));
  empty.SetStats(true);
  empty.Sumw2(true);
  for(const auto &process: processes){
    unique_ptr<SingleHist1D> hist(new SingleHist1D(*this, process, empty));
    hist->raw_hist_.SetFillColor(process->GetFillColor());
    hist->raw_hist_.SetFillStyle(process->GetFillStyle());
    hist->raw_hist_.SetLineColor(process->GetLineColor());
    hist->raw_hist_.SetLineStyle(process->GetLineStyle());
    hist->raw_hist_.SetLineWidth(process->GetLineWidth());
    hist->raw_hist_.SetMarkerColor(process->GetMarkerColor());
    hist->raw_hist_.SetMarkerStyle(process->GetMarkerStyle());
    hist->raw_hist_.SetMarkerSize(process->GetMarkerSize());

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

  blank_.SetFillStyle(0);
  blank_.SetFillColor(kWhite);
  blank_.SetLineWidth(0);
  blank_.SetLineColor(kWhite);
  blank_.SetMarkerSize(0);
  blank_.SetMarkerColor(kWhite);
}

/*! \brief Produce and save formatted plots at given luminosity

  \param[in] luminosity The integrated luminosity with which to draw the plot
*/
void Hist1D::Print(double luminosity,
                   const string &subdir){
  if (!draw_plot_) return;
  luminosity_ = luminosity;
  for(const auto &opt: plot_options_){
    this_opt_ = opt;
    this_opt_.MakeSane();
    gStyle->SetColorModelPS(this_opt_.UseCMYK());
    gROOT->ForceStyle();
    RefreshScaledHistos();
    SetRanges();
    ApplyStyles();
    AdjustFillStyles();

    double bot_min, bot_max;
    vector<TH1D> bot_plots = GetBottomPlots(bot_min, bot_max);
    //I don't know why I can't make this in GetBottomPlots...
    TGraphAsymmErrors bottom_background;
    if(this_opt_.Bottom() != BottomType::off){
      bottom_background = TGraphAsymmErrors(&bot_plots.back());
      bottom_background.SetMinimum(this_opt_.RatioMinimum());
      bottom_background.SetMaximum(this_opt_.RatioMaximum());
      bot_plots.pop_back();
    }

    TGraphAsymmErrors bkg_error = GetBackgroundError();

    StripTopPlotLabels();
    TLine horizontal = GetBottomHorizontal();
    vector<TLine> cut_vals = GetCutLines(GetMinDraw(), GetMaxDraw(), true);
    vector<TLine> bot_cuts = GetCutLines(bot_min, bot_max, false);

    unique_ptr<TCanvas> full;
    unique_ptr<TPad> top, bottom;
    GetPads(full, top, bottom);

    if(this_opt_.AutoYAxis()) FixYAxis(bot_plots);

    if(this_opt_.Bottom() != BottomType::off){
      bottom->cd();


      if(bot_plots.size()>0) bot_plots.at(0).Draw("e0p");
      bottom_background.Draw("2 same");
      string draw_opt = "e0p0 same";
      for(auto &h: bot_plots){
        h.Draw(draw_opt.c_str());
      }


      horizontal.Draw("same");

      bottom->RedrawAxis();
      bottom->RedrawAxis("g");
      for(auto &cut: bot_cuts) cut.Draw();
    }

    top->cd();
    if(this_opt_.YAxis() == YAxisType::log) top->SetLogy(true);
    else top->SetLogy(false);

    string draw_opt = "hist";
    DrawAll(backgrounds_, draw_opt);
    if(this_opt_.ShowBackgroundError() && backgrounds_.size()) bkg_error.Draw("2 same");
    DrawAll(signals_, draw_opt, true);

    ReplaceAll(draw_opt, "hist", "e0p");
    // Turn on error on zero
    vector<TH1::EBinErrorOpt> default_bin_error_option;
    if(this_opt_.ErrorOnZeroData()) {
      for(auto &hist: datas_) {
        default_bin_error_option.push_back(hist->scaled_hist_.GetBinErrorOption());
        hist->scaled_hist_.SetBinErrorOption(TH1::kPoisson);
      }
    }
    DrawAll(datas_, draw_opt, true);
    //// Return to default setting -> This doesn't work. log plots also show error on zero
    //if(this_opt_.ErrorOnZeroData()) {
    //  int iPlot = 0;
    //  for(auto &hist: datas_) hist->scaled_hist_.SetBinErrorOption(default_bin_error_option[iPlot++]);
    //}
    for(auto &cut: cut_vals) cut.Draw();

    vector<shared_ptr<TLegend> > legends = GetLegends();
    for(auto &legend: legends){
      legend->Draw();
    }

    top->RedrawAxis();
    top->RedrawAxis("g");

    if(left_label_.size()>0){
      for (unsigned ilabel(0); ilabel<left_label_.size(); ilabel++) {
        TLatex label; 
        label.SetTextFont(this_opt_.Font()+10);label.SetTextSize(this_opt_.ExtraLabelSize());
        label.SetTextAlign(13);
        size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
        if(this_opt_.DisplayLumiEntry()) ++num_plots;
        double legend_height = this_opt_.TrueLegendHeight(num_plots);
        double left_bound = this_opt_.LeftMargin()+0.05;
        double bottom_bound = 1-legend_height-this_opt_.LegendPad()*2-(ilabel+1)*this_opt_.ExtraLabelSize()-0.02;
        label.DrawLatexNDC(left_bound, bottom_bound, left_label_[ilabel].c_str());
      }
    }

    if(right_label_.size()>0){
      for (unsigned ilabel(0); ilabel<right_label_.size(); ilabel++) {
        TLatex label; 
        label.SetTextFont(this_opt_.Font()+10);label.SetTextSize(this_opt_.ExtraLabelSize());
        label.SetTextAlign(33);
        size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
        if(this_opt_.DisplayLumiEntry()) ++num_plots;
        double legend_height = this_opt_.TrueLegendHeight(num_plots);
        double right_bound = 1-this_opt_.RightMargin()-0.03;
        double bottom_bound = 1-legend_height-0.01-this_opt_.LegendPad()*2-(ilabel+1)*this_opt_.ExtraLabelSize()-0.02;
        label.DrawLatexNDC(right_bound, bottom_bound, right_label_[ilabel].c_str());
      }
    }

    vector<shared_ptr<TLatex> > title_text = GetTitleTexts();
    for(auto &x: title_text){
      x->Draw();
    }

    // Printing values to terminal
    if(this_opt_.PrintVals()){
      TH1D *hdata = (datas_.size() ? &(datas_[0]->scaled_hist_) : 0);
      TH1D *hmc = (backgrounds_.size() ? &(backgrounds_[0]->scaled_hist_) : 0);
      TH1D *hbot = (bot_plots.size() ? &(bot_plots[0]) : 0);
      if(hdata==0 || hmc==0 || hbot==0 || hdata->Integral()==0|| hbot->Integral()==0) cout<<"Printing values failed: no histogram or no entries in histogram"<<endl;
      else {
	int digits = floor(log10(max(hdata->GetBinContent(hdata->GetMaximumBin()), 
				     hmc->GetBinContent(hmc->GetMaximumBin())))+1.);
	//// Digits for error are calculated with the sqrt, and added 2 to print with one decimal
	int edigits = floor(log10(sqrt(max(hdata->GetMaximum(), hmc->GetMaximum())))+1.)+2;
	cout<<endl<<"Printing values for "<<Name()<<". Data/MC = "
	    <<RoundNumber(hdata->Integral(), 2,hmc->Integral()) <<endl;
	for(int bin=1; bin<=hdata->GetNbinsX(); bin++){
	  cout<<"Bin "<<setw(5)<<hdata->GetBinLowEdge(bin)<<","<<setw(5)<<hdata->GetBinLowEdge(bin+1)<<": Data = ";
	  cout<<setw(digits)<<RoundNumber(hdata->GetBinContent(bin),1)<<" +- "<<setw(edigits)<<RoundNumber(hdata->GetBinError(bin),1);
	  cout<<", MC = "<<setw(digits+2)<<RoundNumber(hmc->GetBinContent(bin),1)<<" +- "
	      <<setw(edigits)<<RoundNumber(hmc->GetBinError(bin),1);
	  if(this_opt_.Bottom() != BottomType::off)
	    cout<<"   ->    Ratio  = "<<setw(6)<<RoundNumber(hbot->GetBinContent(bin),3)
		<<" +- "<<setw(5)<<RoundNumber(hbot->GetBinError(bin),3);
	  cout<<endl;
	} // Loop over histogram bins
      }
    }

    if(subdir != "") mkdir(("plots/"+subdir).c_str(), 0777);
    string base_name = subdir != ""
      ? "plots/"+subdir+"/"+Name()
      : "plots/"+Name();
    for(const auto &ext: this_opt_.FileExtensions()){
      // string full_name = base_name+"__"+this_opt_.TypeString()+'.'+ext;
      string full_name = base_name+'.'+ext;
      if (Contains(tag_,"FixName:")){
        string tagName=tag_;
        ReplaceAll(tagName, "FixName:", "");
        base_name = subdir != ""
          ? "plots/"+subdir+"/"+tagName
          : "plots/"+tagName;
        full_name = base_name+'.'+ext;
      }
      full->Print(full_name.c_str());
      cout << "open " << full_name << endl;
    }

  }
}

set<const Process*> Hist1D::GetProcesses() const{
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

std::string Hist1D::GetTag() const {
  return tag_;
}

Figure::FigureComponent * Hist1D::GetComponent(const Process *process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

string Hist1D::Name() const{
  string cut = "";
  if(cut_.Name() != "1") cut = "__"+cut_.Name();

  string weight = "";
  if(weight_.Name() != "weight") weight = "__"+weight_.Name();
  
  if(tag_ == ""){
    return CodeToPlainText(xaxis_.var_.Name()+cut+weight);
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

string Hist1D::Title() const{
  bool cut = (cut_.Name() != "" && cut_.Name() != "1");
  bool weight = weight_.Name() != "weight";
  if(cut && weight){
    return CodeToRootTex(cut_.Name())+" (weight="+CodeToRootTex(weight_.Name())+")";
  }else if(cut){
    return CodeToRootTex(cut_.Name());
  }else if(weight){
    return CodeToRootTex("weight="+weight_.Name());
  }else{
    return "";
  }
}

Hist1D & Hist1D::Weight(const NamedFunc &weight){
  weight_ = weight;
  return *this;
}

Hist1D & Hist1D::Tag(const string &tag){
  tag_ = tag;
  return *this;
}

Hist1D & Hist1D::LuminosityTag(const string &tag){
  luminosity_tag_ = tag;
  return *this;
}

Hist1D & Hist1D::LeftLabel(const vector<string> &label){
  left_label_ = label;
  return *this;
}

Hist1D & Hist1D::RightLabel(const vector<string> &label){
  right_label_ = label;
  return *this;
}

Hist1D & Hist1D::YAxisZoom(const double &yaxis_zoom){
  yaxis_zoom_ = yaxis_zoom;
  return *this;
}

Hist1D & Hist1D::RatioTitle(const string &numerator,
                            const string &denominator){
  ratio_numerator_ = numerator;
  ratio_denominator_ = denominator;
  return *this;
}

Hist1D & Hist1D::DrawPlot(const bool &draw_plot) {
  draw_plot_ = draw_plot;
  return *this;
}

/*!\brief Generates stacked and scaled histograms from unstacked and unscaled
  ones

  Sets bin contents for all required Hist1D::SingleHist1D::scaled_hist_ to the
  appropriate values using the Hist1D::SingleHist1D::raw_hist_ containing the
  unstacked contents at 1 fb^{-1}
*/
void Hist1D::RefreshScaledHistos(){
  InitializeHistos();
  MergeOverflow();
  ScaleHistos();
  StackHistos();
  NormalizeHistos();
  FixAsymmErrors();
}

/*!\brief Sets all Hist1D::SingleHist1D::scaled_hist_ to corresponding
  Hist1D::SingleHist1D::raw_hist_
*/
void Hist1D::InitializeHistos() const{
  for(auto &hist: backgrounds_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("bkg_"+hist->process_->name_+"_"+counter()).c_str());
  }
  for(auto &hist: signals_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("sig_"+hist->process_->name_+"_"+counter()).c_str());
  }
  for(auto &hist: datas_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("dat_"+hist->process_->name_).c_str());
//     hist->scaled_hist_.SetName(("dat_"+hist->process_->name_+"_"+counter()).c_str());
  }
}

/*!\brief Moves overflow (underflow) contents into last (first) visible bin
s */
void Hist1D::MergeOverflow() const{
  bool underflow = false, overflow = false;
  switch(this_opt_.Overflow()){
  default:
  case OverflowType::none:
    underflow = false;
    overflow = false;
    break;
  case OverflowType::underflow:
    underflow = true;
    overflow = false;
    break;
  case OverflowType::overflow:
    underflow = false;
    overflow = true;
    break;
  case OverflowType::both:
    underflow = true;
    overflow = true;
    break;
  }

  for(auto &hist: backgrounds_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
  for(auto &hist: signals_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
  for(auto &hist: datas_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
}

/*!\brief Scales histograms to required luminosity
 */
void Hist1D::ScaleHistos() const{
  for(auto &hist: backgrounds_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
    hist->scaled_hist_.Scale(luminosity_);
  }
  for(auto &hist: signals_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
    hist->scaled_hist_.Scale(luminosity_);
  }
  for(auto &hist: datas_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
  }
}

/*!\brief Stacks histograms if necessary for current plot style
 */
void Hist1D::StackHistos() const{
  if(this_opt_.Stack() == StackType::signal_overlay
     || this_opt_.Stack() == StackType::signal_on_top
     || this_opt_.Stack() == StackType::data_norm){
    for(size_t ibkg = backgrounds_.size() - 2; ibkg < backgrounds_.size(); --ibkg){
      backgrounds_.at(ibkg)->scaled_hist_ = backgrounds_.at(ibkg)->scaled_hist_ + backgrounds_.at(ibkg+1)->scaled_hist_;
    }
    if(backgrounds_.size() && this_opt_.Stack() == StackType::signal_on_top){
      for(auto &hist: signals_){
        hist->scaled_hist_ = hist->scaled_hist_ + backgrounds_.front()->scaled_hist_;
      }
    }
  }
}

/*!\brief Normalize histograms to data or 100%*(bin width) if needed for current
  style
*/
void Hist1D::NormalizeHistos() const{
  mc_scale_ = 1.;
  mc_scale_error_ = 1.;
  if(this_opt_.Stack() == StackType::data_norm){
    if(datas_.size() == 0 || backgrounds_.size() == 0) return;
    int nbins = xaxis_.Nbins();
    double data_error, mc_error;
    double data_norm = datas_.front()->scaled_hist_.IntegralAndError(0, nbins+1, data_error, "width");
    double mc_norm = backgrounds_.front()->scaled_hist_.IntegralAndError(0, nbins+1, mc_error, "width");
    mc_scale_ = data_norm/mc_norm;
    if(this_opt_.PrintVals()){
      cout << "MC scale factor: " << mc_scale_ << endl;
    }
    mc_scale_error_ = hypot(data_norm*mc_error, mc_norm*data_error)/(mc_norm*mc_norm);
    for(auto &hist: backgrounds_){
      hist->scaled_hist_.Scale(mc_scale_);
    }
    for(auto h = datas_.begin(); h != datas_.end(); ++h){
      auto &hist = *h;
      double dumb;
      double this_integral = hist->scaled_hist_.IntegralAndError(0, nbins+1, dumb, "width");
      hist->scaled_hist_.Scale(this_integral == 0. ? 1. : data_norm/this_integral);
    }
  }else if(this_opt_.Stack() == StackType::shapes){
    for(auto &hist: backgrounds_){
      Normalize(hist->scaled_hist_, 100., true);
    }
    for(auto &hist: signals_){
      Normalize(hist->scaled_hist_, 100., true);
    }
    for(auto &hist: datas_){
      Normalize(hist->scaled_hist_, 100., true);
    }
  }
}

/*!\brief Restore asymmetric error bars if histogram unweighted"
 */
void Hist1D::FixAsymmErrors() const{
  for(const auto &hist: backgrounds_){
    FixAsymmErrors(hist);
  }
  for(const auto &hist: signals_){
    FixAsymmErrors(hist);
  }
  for(unsigned data_idx = 0; data_idx < datas_.size(); data_idx++) {
    bool error_on_zero = (data_idx == (datas_.size()-1)) && this_opt_.ErrorOnZeroData();
    FixAsymmErrors(datas_[data_idx], error_on_zero);
  }
}

/*!\brief Restore asymmetric error bars if histogram unweighted"
 */
void Hist1D::FixAsymmErrors(const unique_ptr<SingleHist1D> &sh1d, bool error_on_zero_data) const{
  TH1D &h = sh1d->scaled_hist_;
  const TArrayD *sumw2 = h.GetSumw2();

  double tolerance = 4.*numeric_limits<double>::epsilon();

  bool turn_off_sumw2 = true;
  for(int bin = 0; turn_off_sumw2 && bin < h.GetNcells(); ++bin){
    double w2 = sumw2->At(bin);
    double content = h.GetBinContent(bin);
    if(fabs(w2-content) > tolerance*content){
      turn_off_sumw2 = false;
    }
  }

  if(turn_off_sumw2 || error_on_zero_data) {
    h.Sumw2(false);
  }
}

/*!\brief Set y-axis plotting range
 */
void Hist1D::SetRanges() const{
  double the_min = GetMinDraw();
  double the_max = GetMaxDraw()/yaxis_zoom_;

  double ratio = GetLegendRatio();

  double top, bottom;
  switch(this_opt_.YAxis()){
  default:
  case YAxisType::linear:
    bottom = the_min >= 0. ? 0. : the_min;
    top = bottom+ratio*(the_max-bottom);
    break;
  case YAxisType::log:
    bottom = the_min > this_opt_.LogMinimum() ? the_min : this_opt_.LogMinimum();
    top = exp(log(bottom)+ratio*(log(the_max)-log(bottom)));
    break;
  }

  for(auto &hist: backgrounds_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: signals_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: datas_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
}

/*!\brief Set label styles and title for all histograms
 */
void Hist1D::ApplyStyles() const{
  for(auto &hist: backgrounds_){
    StyleHisto(hist->scaled_hist_);
  }
  for(auto &hist: signals_){
    StyleHisto(hist->scaled_hist_);
  }
  for(auto &hist: datas_){
    StyleHisto(hist->scaled_hist_);
  }
}

/*!\brief Set label styles and title for a histogram

  \param[in,out] h Histogram to be adjusted
*/
void Hist1D::StyleHisto(TH1D &h) const{
  h.GetXaxis()->SetTitleOffset(this_opt_.XTitleOffset());
  h.GetYaxis()->SetTitleOffset(this_opt_.YTitleOffset());
  h.SetNdivisions(this_opt_.NDivisions(), "xyz");
  h.SetLabelSize(this_opt_.LabelSize(), "xyz");
  h.SetTitleSize(this_opt_.TitleSize(), "xyz");
  h.SetLabelFont(this_opt_.Font(), "xyz");
  h.SetTitleFont(this_opt_.Font(), "xyz");
  
  double bin_width = xaxis_.AvgBinWidth();

  bool have_vec_ = false; // To be properly filled from the SingleHist1D...
  string Yunit = (have_vec_?"Entries":"Events");
  string yunit = (have_vec_?"entries":"events");
  
  ostringstream title;
  switch(this_opt_.Stack()){
  default:
    DBG("Unrecognized stack option " << static_cast<int>(this_opt_.Stack()) << ".");
    /* FALLTHRU */
  case StackType::signal_overlay:
    /* FALLTHRU */
  case StackType::signal_on_top:
    /* FALLTHRU */
  case StackType::data_norm:
    /* FALLTHRU */
  case StackType::lumi_shapes:
    if(xaxis_.units_ == "" && bin_width == 1){
      title << Yunit;    
      break;
    }
    else{
      title << Yunit<<" / " << bin_width;
      if(xaxis_.units_ != "") title << " " << xaxis_.units_;
      break;
    }
  case StackType::shapes:
    if(this_opt_.Title() == TitleType::simulation_supplementary){
      if(xaxis_.units_ == "" && bin_width == 1){
	title << "% "<<yunit;
	break;
      }
      else{
	title << "% "<<yunit<<" / " << bin_width;
	if(xaxis_.units_ != "") title << " " << xaxis_.units_;
	break;
      }
    }
    else{  
      if(xaxis_.units_ == "" && bin_width == 1){
	title << "% of "<<yunit;
	break;
      }
      else{
	title << "% of "<<yunit<<" / " << bin_width;
	if(xaxis_.units_ != "") title << " " << xaxis_.units_;
	break;
      }
    }
  }
  
  h.GetYaxis()->SetTitle(title.str().c_str());
}

/*!\brief Make histograms a hollow line for unstacked styles
 */
void Hist1D::AdjustFillStyles() const{
  if(this_opt_.BackgroundsStacked()) return;

  for(auto &bkg: backgrounds_){
    TH1D &h = bkg->scaled_hist_;
    h.SetFillStyle(0);
    h.SetLineColor(h.GetFillColor());
    h.SetLineWidth(3);
  }
}

/*!\brief Generated canvas and pads for top and bottom plots

  To make object positioning simpler, the top and bottom pads span the entire
  canvas and only the margins are adjusted. This has the added bonus of making
  the "0" on the top y-axis visible for free.

  \param[out] c Full plot canvas

  \param[out] top Pad for main plot

  \param[out] bottom Pad for bottom plot (e.g. ratio plot)
*/
void Hist1D::GetPads(unique_ptr<TCanvas> &c,
                     unique_ptr<TPad> &top,
                     unique_ptr<TPad> &bottom) const{
  c.reset(new TCanvas("canvas", "canvas", this_opt_.CanvasWidth(), this_opt_.CanvasHeight()));
  c->cd();
  top.reset(new TPad(("top_pad_"+counter()).c_str(), "top_pad", 0., 0., 1., 1.));
  bottom.reset(new TPad(("bottom_pad_"+counter()).c_str(), "bottom_pad", 0., 0., 1., 1.));
  c->SetMargin(0., 0., 0., 0.);
  c->SetTicks(1,1);
  c->SetFillStyle(4000);
  top->SetTicks(1,1);
  top->SetFillStyle(4000);
  bottom->SetTicks(1,1);
  bottom->SetFillStyle(4000);
  if(this_opt_.Bottom() == BottomType::off){
    top->SetMargin(this_opt_.LeftMargin(),
                   this_opt_.RightMargin(),
                   this_opt_.BottomMargin(),
                   this_opt_.TopMargin());
  }else{
    top->SetMargin(this_opt_.LeftMargin(),
                   this_opt_.RightMargin(),
                   this_opt_.BottomHeight(),
                   this_opt_.TopMargin());
    bottom->SetMargin(this_opt_.LeftMargin(),
                      this_opt_.RightMargin(),
                      this_opt_.BottomMargin(),
                      1.-this_opt_.BottomHeight());
  }
  bottom->Draw();
  top->Draw();
}

/*!\brief Adjust y-axis title offset based on y-axis range

  \param[in,out] bottom_plots Ratio or other bottom-half plots whose title
  offsets also need to be adjusted
*/
void Hist1D::FixYAxis(vector<TH1D> &bottom_plots) const{
  double offset = this_opt_.YTitleOffset();
  if(this_opt_.YAxis() == YAxisType::log){
    offset = 1.5;
  }else{
    double the_max = GetMaxDraw()*GetLegendRatio();
    int digits = fabs(floor(log10(the_max))-1)+2;
    digits = max(2, min(6, digits));//For huge axis scale, reverts to scientific notation

    //Scale offset by good empirical numbers
    offset = 0.6+0.25*digits;
  
  }
  for(auto &hist: backgrounds_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: signals_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: datas_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: bottom_plots){
    if(this_opt_.Bottom() != BottomType::ratio) hist.SetTitleOffset(offset, "y");
  }
}

/*!\brief Get text to print at top of plot

  Depending on current plot style, this may be the CMS {Preliminary, Simulation,
  etc...} with luminosity or the cut and weight for the plot

  \return List of text items to be printed in title region of plot
*/
vector<shared_ptr<TLatex> > Hist1D::GetTitleTexts() const{
  vector<shared_ptr<TLatex> > out;
  double left = this_opt_.LeftMargin();
  double right = 1.-this_opt_.RightMargin();
  double bottom = 1.-this_opt_.TopMargin();
  double top = 1.;
  if(this_opt_.Title() == TitleType::info){
    if(Title() != ""){
      out.push_back(make_shared<TLatex>(0.5*(left+right), 0.5*(bottom+top),
                                        Title().c_str()));
      out.back()->SetNDC();
      out.back()->SetTextAlign(22);
      out.back()->SetTextFont(this_opt_.Font());

      //Adjust title to fit in available space
      double max_width, max_height;
      GetTitleSize(max_width, max_height, true);
      UInt_t width, height;
      out.back()->GetBoundingBox(width, height);
      while(width > max_width || height > max_height){
        out.back()->SetTextSize(0.8*out.back()->GetTextSize());
        out.back()->GetBoundingBox(width, height);
      }
      while(width < 0.5*max_width && height < 0.5*max_height){
        out.back()->SetTextSize(1.25*out.back()->GetTextSize());
        out.back()->GetBoundingBox(width, height);
      }
    }
  }else{
    string extra;
    switch(this_opt_.Title()){
    case TitleType::preliminary: extra = "Preliminary"; break;
    case TitleType::simulation: extra = "Simulation"; break;
    case TitleType::simulation_preliminary: extra = "Simulation Preliminary"; break;
    case TitleType::simulation_supplementary: extra = "Simulation Supplementary"; break;
    case TitleType::supplementary: extra = "Supplementary"; break;
    case TitleType::data: extra = ""; break;
    case TitleType::info:
    default:
      ERROR("Did not understand title type "+to_string(static_cast<int>(this_opt_.Title())));
    }
   
    if (this_opt_.TitleInFrame()) {
      double in_frame_bottom = 1.-this_opt_.TopMargin()-this_opt_.LegendPad()-this_opt_.TitleSize();
      double in_frame_left = this_opt_.LeftMargin()+this_opt_.LegendPad();
      out.push_back(make_shared<TLatex>(in_frame_left, in_frame_bottom,
		  		      ("#font[62]{CMS}#scale[0.74]{#font[52]{ "+extra+"}}").c_str()));
    }
    else
      out.push_back(make_shared<TLatex>(left, bottom+0.2*(top-bottom),
		  		      ("#font[62]{CMS}#scale[0.74]{#font[52]{ "+extra+"}}").c_str()));
    
   
    out.back()->SetNDC();
    out.back()->SetTextAlign(11);
    out.back()->SetTextFont(this_opt_.Font());
    out.back()->SetTextSize(this_opt_.TitleSize());

    ostringstream oss;
    if(this_opt_.Stack() != StackType::shapes) {
      if (luminosity_tag_ != "") oss << luminosity_tag_ << " fb^{-1} (13 TeV)" << flush;
      else if (luminosity_<1.1) oss << "137 fb^{-1} (13 TeV)" << setprecision(1) << flush;
      else oss << setprecision(3) << luminosity_ << " fb^{-1} (13 TeV)" << flush;
    } else oss << "13 TeV" << flush;
    out.push_back(make_shared<TLatex>(right, bottom+0.2*(top-bottom),
                                      oss.str().c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(31);
    out.back()->SetTextFont(this_opt_.Font());
    out.back()->SetTextSize(this_opt_.TitleSize());
  }
  return out;
}

/*!\brief Get uncertainty on total background

  \return Graph representing the MC background error band
*/
TGraphAsymmErrors Hist1D::GetBackgroundError() const{
  TGraphAsymmErrors g;
  if(backgrounds_.size() == 0){
    TH1D h("", "", xaxis_.Nbins(), &xaxis_.Bins().at(0));
    g = TGraphAsymmErrors(&h);
  }else{
    g = TGraphAsymmErrors(&(backgrounds_.front()->scaled_hist_));
  // set the color of the error band to the line color, accomodating data-to-data plots
    g.SetFillColorAlpha(backgrounds_.front()->scaled_hist_.GetLineColor(),0.2);
  }
  //g.SetFillStyle(3002);
  g.SetLineWidth(0);
  g.SetMarkerSize(0);
  return g;
}

/*!\brief Get vertical lines at cut values

  \param[in] y_min Lower bound of y-axis

  \param[in] y_max Upper bound of y-axis

  \return Lines at x-coordinate of cut value and y-coordinates running from
  bottom of plot to bottom of legend
*/
vector<TLine> Hist1D::GetCutLines(double y_min, double y_max, bool adjust_bottom) const{
  double bottom = y_min;
  if(adjust_bottom){
    switch(this_opt_.YAxis()){
    default:
      DBG("Bad YAxis type " << static_cast<int>(this_opt_.YAxis()));
      /* FALLTHRU */
    case YAxisType::linear: bottom = y_min >= 0. ? 0. : y_min; break;
    case YAxisType::log:    bottom = y_min > this_opt_.LogMinimum() ? y_min : this_opt_.LogMinimum(); break;
    }
  }
  vector<TLine> out(xaxis_.cut_vals_.size());
  for(double cut: xaxis_.cut_vals_){
    out.emplace_back(cut, bottom, cut, y_max);
    out.back().SetNDC(false);
    out.back().SetLineStyle(2);
    out.back().SetLineColor(kBlack);
    out.back().SetLineWidth(3);
  }
  for(double cut: xaxis_.hard_cut_vals_){
    out.emplace_back(cut, bottom, cut, y_max);
    out.back().SetNDC(false);
    out.back().SetLineStyle(9);
    out.back().SetLineColor(kBlack);
    out.back().SetLineWidth(3);
  }

  return out;
}

/*!\brief Get ratio or other plots drawn on the lower pad

  \param [out] the_min Y-axis minimum across plots for lower pad

  \param [out] the_max Y-axis maximum across plots for lower pad

  \return Set of plots to be drawn on lower pad. These may be ratio plots or
  something else depending on the current plot style
*/
std::vector<TH1D> Hist1D::GetBottomPlots(double &the_min, double &the_max) const{
  if(this_opt_.Bottom() == BottomType::off) return vector<TH1D>();

  TH1D denom;
  vector<TH1D> out;
  if(backgrounds_.size() != 0){
    denom = backgrounds_.front()->scaled_hist_;
  }else if(datas_.size() != 0){
    denom = datas_.front()->scaled_hist_;
  }else if(signals_.size() != 0){
    denom = signals_.front()->scaled_hist_;
  }else{
    ERROR("No histograms available to make bottom plot");
  }
  bool stacked;
  switch(this_opt_.Stack()){
  case StackType::signal_overlay:
  case StackType::signal_on_top:
  case StackType::data_norm:
    stacked = true; break;
  case StackType::lumi_shapes:
  case StackType::shapes:
    stacked = false; break;
  default:
    ERROR("Bad stack type: "+to_string(static_cast<int>(this_opt_.Stack())));
    break;
  }
  if(stacked && backgrounds_.size()){
    out.push_back(backgrounds_.front()->scaled_hist_);
    out.back().SetName(("bot_plot_bkg_"+backgrounds_.front()->process_->name_+"_"+counter()).c_str());
  }else{
    for(const auto &h: backgrounds_){
      out.push_back(h->scaled_hist_);
      out.back().SetName(("bot_plot_bkg_"+h->process_->name_+"_"+counter()).c_str());
    }
  }
  for(const auto &h: datas_){
    out.push_back(h->scaled_hist_);
    out.back().SetName(("bot_plot_data_"+h->process_->name_+"_"+counter()).c_str());
  }
  if(!stacked || this_opt_.Bottom() == BottomType::sorb || this_opt_.Bottom() == BottomType::sorb_cut_upper){
    for(const auto &h: signals_){
      out.push_back(h->scaled_hist_);
      out.back().SetName(("bot_plot_sig_"+h->process_->name_+"_"+counter()).c_str());
    }
  }
  if(!out.size()) return vector<TH1D>();
  TH1D band = out.front();
  for(size_t i = 0; (i+1) < out.size(); ++i){
    out.at(i) = out.at(i+1);
  }
  out.back() = band;
  out.back().SetFillStyle(1001);
  out.back().SetFillColorAlpha(backgrounds_.front()->scaled_hist_.GetLineColor(),0.2);
  out.back().SetLineWidth(0);
  out.back().SetMarkerStyle(0);
  out.back().SetMarkerSize(0);
  out.back().SetName(("bot_plot_band_"+counter()).c_str());

  for(int bin = 0; bin <= denom.GetNbinsX()+1; ++bin){
    denom.SetBinError(bin, 0.);
  }

  switch(this_opt_.Bottom()){
  case BottomType::ratio:
    for(auto &h: out){
      //commented code makes error bars on zero entry bins
      //h.Sumw2(false);
      //h.SetBinErrorOption(TH1::kPoisson);
      h.Divide(&denom);
      //for(int bin = 0; bin <= denom.GetNbinsX()+1; ++bin){
      //  if (h.GetBinContent(bin) < 1.0e-12) {
      //    if (denom.GetBinContent(bin) > 0) {
      //      h.SetBinContent(bin, 0);
      //      h.SetBinError(bin, 1.8/denom.GetBinContent(bin));
      //    }
      //    else
      //      //hack to avoid seeing any plotted points for bins with no signal or background
      //      h.SetBinContent(bin, -1);
      //  }
      //}
    }
    break;
  case BottomType::diff:
    for(auto &h: out){
      h = h - denom;
    }
    break;
  case BottomType::sorb:
    for(auto &h: out){
      double totb = denom.Integral(0,denom.GetNbinsX()+1);
      if (totb < 1.0e-10) totb = 1.0; //avoid divide-by-zeros
      double tots = h.Integral(0,h.GetNbinsX()+1);
      if (tots < 1.0e-10) tots = 1.0; //avoid divide-by-zeros
      for(int bin = 0; bin <= h.GetNbinsX()+1; ++bin) {
        double total_background_yield = denom.Integral(bin,denom.GetNbinsX()+1)/totb;
        if (total_background_yield < 1.0e-10) { //avoid divide-by-zeros and negative yields
          h.SetBinContent(bin, 0);
          h.SetBinError(bin, 0);
        }
        else {
          double sqrtb = sqrt(total_background_yield);
          double herr(0); 
          for(int ebin = bin; ebin <= h.GetNbinsX()+1; ++ebin) 
            herr += pow(h.GetBinError(ebin),2);
          herr = sqrt(herr)/sqrtb;
          h.SetBinContent(bin, h.Integral(bin,h.GetNbinsX()+1)/tots/sqrtb);
          h.SetBinError(bin, herr/tots);
        }
      }
    }
    break;
  case BottomType::sorb_cut_upper:
    for(auto &h: out){
      double totb = denom.Integral(0,denom.GetNbinsX()+1);
      if (totb < 1.0e-10) totb = 1.0; //avoid divide-by-zeros
      double tots = h.Integral(0,h.GetNbinsX()+1);
      if (tots < 1.0e-10) tots = 1.0; //avoid divide-by-zeros
      for(int bin = h.GetNbinsX()+1; bin >= 0; --bin) {
        double total_background_yield = denom.Integral(0,bin)/totb;
        if (total_background_yield < 1.0e-10) { //avoid divide-by-zeros and negative yields
          h.SetBinContent(bin, 0.0);
          h.SetBinError(bin, 0.0);
        }
        else {
          double sqrtb = sqrt(total_background_yield);
          double herr(0); 
          for(int ebin = bin; ebin >= 0; --ebin) 
            herr += pow(h.GetBinError(ebin),2);
          herr = sqrt(herr)/sqrtb;
          h.SetBinContent(bin, h.Integral(0,bin)/tots/sqrtb);
          h.SetBinError(bin, herr/tots);
        }
      }
    }
    break;
  case BottomType::off:
  default:
    ERROR("Bad type for bottom plot: "+to_string(static_cast<int>(this_opt_.Bottom())));
    break;
  }

  the_min = numeric_limits<double>::infinity();
  the_max = -numeric_limits<double>::infinity();
  for(auto &h: out){
    h.SetNdivisions(this_opt_.NDivisionsBottom(), "y");
    for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
      double hi = h.GetBinContent(bin)+h.GetBinErrorUp(bin);
      double lo = h.GetBinContent(bin)-fabs(h.GetBinErrorLow(bin));
      if(hi>the_max) the_max = hi;
      if(lo<the_min) the_min = lo;
    }
  }

  //make bottom plot title
  if(this_opt_.Bottom() == BottomType::ratio){
    the_min = this_opt_.RatioMinimum();
    the_max = this_opt_.RatioMaximum();
    for(auto &h: out){
      string num = ratio_numerator_;
      string den = ratio_denominator_;
      if(datas_.size() != 0){
        if(num == "") num = "Data";
        if(den == "") den = "MC";
      }else{
        if(num == "") num = "MC";
        if(den == "") den = backgrounds_.front()->process_->name_;
      }
      h.GetYaxis()->CenterTitle();
      h.GetYaxis()->SetTitle(("#frac{"+num+"}{"+den+"}").c_str());
      h.SetTitleSize(h.GetTitleSize("y")/1.25,"y");
      h.SetTitleOffset(h.GetTitleOffset("y")/1.25, "y");
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
  }else if(this_opt_.Bottom() == BottomType::diff){
    for(auto &h: out){
      if(datas_.size() != 0) h.GetYaxis()->SetTitle("Data-MC");
      //else h.GetYaxis()->SetTitle(backgrounds_.front()->process_->name_.c_str());
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
      h.GetYaxis()->CenterTitle();
    }
  }else if(this_opt_.Bottom() == BottomType::sorb){
    the_min = 0;
    the_max = 1.2*the_max > 1 ? 1.2*the_max : 1;
    for(auto &h: out){
      if(signals_.size() != 0) h.GetYaxis()->SetTitle("#frac{#varepsilon_{S}}{#sqrt{#varepsilon_{B}}} (lower cut)");
      h.GetYaxis()->CenterTitle();
      h.SetTitleSize(h.GetTitleSize("y")/1.25,"y");
      h.SetTitleOffset(h.GetTitleOffset("y"), "y");
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
  }else if(this_opt_.Bottom() == BottomType::sorb_cut_upper){
    the_min = 0;
    the_max = 1.2*the_max > 1 ? 1.2*the_max : 1;
    for(auto &h: out){
      if(signals_.size() != 0) h.GetYaxis()->SetTitle("#frac{#varepsilon_{S}}{#sqrt{#varepsilon_{B}}} (upper cut)");
      h.GetYaxis()->CenterTitle();
      h.SetTitleSize(h.GetTitleSize("y")/1.25,"y");
      h.SetTitleOffset(h.GetTitleOffset("y"), "y");
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
  }
  return out;
}

/*!\brief Get horizontal line drawn at "agreement" value for bottom plots

  E.g. Line is at 1 for ratio plots, 0 for difference plots, etc.

  \return Line at appropriate height depending on plot style
*/
TLine Hist1D::GetBottomHorizontal() const{
  double left = xaxis_.Bins().front();
  double right = xaxis_.Bins().back();
  double y;
  switch(this_opt_.Bottom()){
  case BottomType::ratio: y = 1.; break;
  case BottomType::diff: y = 0.; break;
  case BottomType::sorb: y = -1000; break;
  case BottomType::sorb_cut_upper: y = -1000; break;
  case BottomType::off: y = 0.; break;
  default:
    y = 0.;
    DBG("Invalid BottomType: " << to_string(static_cast<int>(this_opt_.Bottom())));
  }

  TLine line(left, y, right, y);
  line.SetNDC(false);
  line.SetLineStyle(2);
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  return line;
}

/*!Remove x-axis labels and title from plots in top pad if necessary
 */
void Hist1D::StripTopPlotLabels() const{
  if(this_opt_.Bottom() == BottomType::off) return;
  for(auto &hist: backgrounds_){
    StripXLabels(hist->scaled_hist_);
  }
  for(auto &hist: signals_){
    StripXLabels(hist->scaled_hist_);
  }
  for(auto &hist: datas_){
    StripXLabels(hist->scaled_hist_);
  }
}

/*!\brief Get highest drawn point below max_bound across all component
  histograms

  \param[in] max_bound Only consider points below this value in finding the
  maximum

  \return The highest drawn point below max_bound across all component
  histograms
*/
double Hist1D::GetMaxDraw(double max_bound) const{
  double the_max = -numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_max = hist->GetMax(max_bound, this_opt_.ShowBackgroundError());
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: signals_){
    double this_max = hist->GetMax(max_bound, false);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: datas_){
    double this_max = hist->GetMax(max_bound, true);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  return the_max;
}

/*!\brief Get lowest drawn point above min_bound across all component histograms

  \param[in] min_bound Only consider points above this value in finding the
  minimum

  \return The lowest drawn point greater than min_bound across all component
  histograms
*/
double Hist1D::GetMinDraw(double min_bound) const{
  double the_min = numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_min = hist->GetMin(min_bound, this_opt_.ShowBackgroundError());
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: signals_){
    double this_min = hist->GetMin(min_bound, false);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: datas_){
    double this_min = hist->GetMin(min_bound, true);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  return the_min;
}

/*!\brief Get list of legends emulating single legend with multiple columns

  Legends are filled down-columns first, then across rows. Data samples are
  added first, then signals, then backgrounds, with the order within each group
  preserved from Hist1D::Hist1D()

  \return Legends with all processes and possibly MC normalization if plot style
  requires it
*/
vector<shared_ptr<TLegend> > Hist1D::GetLegends(){
  size_t n_entries = datas_.size() + signals_.size() + backgrounds_.size();
  if(this_opt_.DisplayLumiEntry()) ++n_entries;
  size_t n_columns = min(n_entries, static_cast<size_t>(this_opt_.LegendColumns()));

  double left = this_opt_.LeftMargin()+this_opt_.LegendLeftPad()+this_opt_.LegendPad();
  double top = 1.-this_opt_.TopMargin()-this_opt_.LegendPad();
  double bottom = top-this_opt_.TrueLegendHeight(n_entries);
  if (this_opt_.TitleInFrame()) left += 0.3;

  double delta_x = this_opt_.TrueLegendWidth(n_entries);
  vector<shared_ptr<TLegend> > legends(n_columns);
  for(size_t i = 0; i < n_columns; ++i){
    double left_column_offset = (i == 0) ? this_opt_.LegendLeftColumnOffset() : 0;
    double x = left+i*delta_x+left_column_offset;
    legends.at(i) = make_shared<TLegend>(x, bottom, x+this_opt_.LegendMarkerWidth(), top);
    legends.at(i)->SetFillStyle(0);
    legends.at(i)->SetBorderSize(0);
    legends.at(i)->SetTextSize(this_opt_.TrueLegendEntryHeight(n_entries));
    legends.at(i)->SetTextFont(this_opt_.Font());
  }

  size_t entries_added = 0;
  AddEntries(legends, datas_, "lep", n_entries, entries_added);
  AddEntries(legends, backgrounds_, this_opt_.BackgroundsStacked() ? "f" : "l", n_entries, entries_added);
  AddEntries(legends, signals_, "l", n_entries, entries_added);
  //  AddEntries(legends, backgrounds_, this_opt_.BackgroundsStacked() ? "f" : "l", n_entries, entries_added);

  //Add a dummy legend entry to display MC normalization
  if(this_opt_.DisplayLumiEntry()){
    auto &leg = legends.at(GetLegendIndex(entries_added, n_entries, legends.size()));
    ostringstream label;
    if (luminosity_tag_ != "") label << fixed  << "L=" << setprecision(3) << luminosity_tag_ << " fb^{-1}";
    else if(luminosity_ != 1.0) label << fixed  << "L=" << setprecision(3) << luminosity_ << " fb^{-1}";
    else label << fixed << setprecision(1) << "L=137 fb^{-1}";
    //else label << fixed << setprecision(1) << "L=" << 36.8 << " fb^{-1}";
    if(this_opt_.Stack() == StackType::data_norm && datas_.size() > 0){
      label << ", (" << 100.*mc_scale_ << "#pm" << 100.*mc_scale_error_ << ")%";
    }
    auto entry = leg->AddEntry(&blank_, label.str().c_str(), "f");
    entry->SetFillStyle(0);
    entry->SetFillColor(kWhite);
    entry->SetLineWidth(0);
    entry->SetLineColor(kWhite);
    entry->SetMarkerStyle(0);
    entry->SetMarkerStyle(kWhite);
  }

  return legends;
}

/*!\brief Distribute processes from list of histograms across legends

  \param[in,out] legends Legends to which entries are added

  \param[in] hists Component histograms whose processes go in legends

  \param[in] style Option string specifying legend marker type. Subset of "flep"

  \param[in] n_entries Number of entries that need to fit in legends

  \param[in,out] entries_added Entries already in legend
*/
void Hist1D::AddEntries(vector<shared_ptr<TLegend> > &legends,
                        const vector<unique_ptr<SingleHist1D> > &hists,
                        const string &style,
                        size_t n_entries,
                        size_t &entries_added) const{
  for(auto h = hists.cbegin(); h != hists.cend(); ++h){
    size_t legend_index = GetLegendIndex(entries_added, n_entries, legends.size());
    TLegend &legend = *legends.at(legend_index);
    string label = (*h)->process_->name_.c_str();
    if(this_opt_.Title() == TitleType::info){
      double value;
      switch(this_opt_.Stack()){
      default:
        DBG("Bad stack option: " << static_cast<int>(this_opt_.Stack()));
        /* FALLTHRU */
      case StackType::signal_overlay:
        /* FALLTHRU */
      case StackType::signal_on_top:
        /* FALLTHRU */
      case StackType::data_norm:
        value = GetYield(h);
        if(value>=1.){
          label += " [N=" + FixedDigits(value, 2) + "]";
        }else{
          label += " [N=" + FixedDigits(value, 1) + "]";
        }
        break;
      case StackType::lumi_shapes:
        /* FALLTHRU */
      case StackType::shapes:
        value = GetMean(h);
        label += " [#mu=" + FixedDigits(value,3) + "]";
        break;
      }
    }
    //Shrink text size if label is long
    double fudge_factor = 0.25;//Not sure how TLegend width affects marker width, but this seems to work
    double max_width = (this_opt_.TrueLegendWidth(n_entries)-this_opt_.LegendMarkerWidth()*fudge_factor) * this_opt_.CanvasWidth();
    double max_height = this_opt_.TrueLegendEntryHeight(n_entries) * this_opt_.LegendDensity() * this_opt_.CanvasHeight();
    TLatex latex(0.5, 0.5, label.c_str());
    latex.SetTextSize(legend.GetTextSize());
    latex.SetTextFont(legend.GetTextFont());
    UInt_t width, height;
    latex.GetBoundingBox(width, height);
    while(width > max_width || height > max_height){
      latex.SetTextSize(0.95*latex.GetTextSize());
      for(auto &leg: legends){
        leg->SetTextSize(0.95*leg->GetTextSize());
      }
      latex.GetBoundingBox(width, height);
    }

    legend.AddEntry(&((*h)->scaled_hist_), label.c_str(), style.c_str());
    ++entries_added;
  }
}

/*!\brief Get factor by which to expand y-axis range to fit legend

  \return Factor by which the upper bound of the top plot's y-axis needs to be
  expanded to make room for the legend
*/
double Hist1D::GetLegendRatio() const{
  size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
  if(this_opt_.DisplayLumiEntry()) ++num_plots;
  double legend_height = this_opt_.TrueLegendHeight(num_plots);
  double top_plot_height;
  if(this_opt_.Bottom() == BottomType::off){
    top_plot_height = 1.-this_opt_.TopMargin()-this_opt_.BottomMargin();
  }else{
    top_plot_height = 1.-this_opt_.TopMargin()-this_opt_.BottomHeight();
  }
  return top_plot_height/(top_plot_height-legend_height-2.*this_opt_.LegendPad());
}

/*!\brief Get integrated number of weighted entries in histogram

  Possibly varying bin widths are accounted for.

  \param[in] h Iterator to histogram in one of Hist1D::datas_,
  Hist1D::signals_, Hist1D::backgrounds_ for which total yield will be
  found
*/
double Hist1D::GetYield(std::vector<std::unique_ptr<SingleHist1D> >::const_iterator h) const{
  TH1D hist = (*h)->scaled_hist_;

  //Subtract underlying histogram
  if((*h)->process_->type_ == Process::Type::background
     && h != (--backgrounds_.cend())
     && this_opt_.BackgroundsStacked()){
    hist = hist - (*(++h))->scaled_hist_;
  }

  //Want yield, not area, so divide out average bin width
  //N.B.: Can't just use hist.Integral() in case of varying bin width
  double raw_integral = hist.Integral("width");
  int nbins = hist.GetNbinsX();
  double left = hist.GetBinLowEdge(1);
  double right = hist.GetBinLowEdge(nbins+1);
  double width = (right-left)/nbins;
  return raw_integral/width;
}

/*!\brief Get mean of histogram

  Possibly varying bin widths are accounted for.

  \param[in] h Iterator to histogram in one of Hist1D::datas_,
  Hist1D::signals_, Hist1D::backgrounds_ for which mean will be found
*/
double Hist1D::GetMean(std::vector<std::unique_ptr<SingleHist1D> >::const_iterator h) const{
  TH1D hist = (*h)->scaled_hist_;

  //Subtract underlying histogram
  if((*h)->process_->type_ == Process::Type::background
     && h != (--backgrounds_.cend())
     && this_opt_.BackgroundsStacked()){
    hist = hist - (*(++h))->scaled_hist_;
  }

  return hist.GetMean();
}

/*!\brief Get width and height of title region

  \param[out] width Width of title region

  \param[out] height Height of title region

  \param[in] in_pixels If true, measure dimensions in pixels. If false, measure
  in NDC.
*/
void Hist1D::GetTitleSize(double &width, double &height, bool in_pixels) const{
  width = 1.-this_opt_.LeftMargin()-this_opt_.RightMargin();
  height = this_opt_.TopMargin();
  if(in_pixels){
    width *= this_opt_.CanvasWidth();
    height *= this_opt_.CanvasHeight();
  }
}

const vector<unique_ptr<Hist1D::SingleHist1D> >& Hist1D::GetComponentList(const Process *process){
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
