#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>

#include <TCanvas.h>
#include <TAxis.h>
#include <TColor.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>

#include "string.hh"
#include "algebra.hh"
#include "type_traits.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::get;
using std::tie;
using namespace ivanp;
using namespace ivanp::math;

using h_t = TH1F;
using h_ptr = std::unique_ptr<h_t>;

h_ptr make_band(
  const std::vector<double>& bins,
  const std::vector<double>& height
) {
  h_t *h = new h_t("","",bins.size()-1,bins.data());
  for (unsigned i=0, n=height.size(); i<n; ++i) {
    h->SetBinError(i+1,height[i]);
  }
  h->SetStats(0);
  h->SetMarkerStyle(0);
  h->SetLineWidth(1); // gives legend color boxes outlines

  return h_ptr(h);
}

std::array<h_ptr,2> make_outline(h_t* h) {
  auto* xa = h->GetXaxis();
  const unsigned nbins = h->GetNbinsX();
  std::array<h_ptr,2> hh {
    h_ptr(new h_t("","",nbins,xa->GetXbins()->GetArray())),
    h_ptr(new h_t("","",nbins,xa->GetXbins()->GetArray()))
  };
  for (unsigned i=1; i<=nbins; ++i) {
    const auto x = h->GetBinError(i);
    get<0>(hh)->SetBinContent(i, x);
    get<1>(hh)->SetBinContent(i,-x);
  }
  for (auto& a : hh) {
    a->SetMarkerStyle(0);
    a->SetLineWidth(1);
    a->SetLineColor(1);
    a->SetLineStyle(h->GetLineStyle());
    a->SetLineColor(h->GetLineColor());
  }
  return std::move(hh);
}

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " input.HepData" << endl;
    return 1;
  }
  bool burst = false, corr = false;
  for (int i=2; i<argc; ++i) {
    if (!strcmp(argv[i],"burst")) burst = true;
    else if (!strcmp(argv[i],"corr")) corr = true;
  }

  struct bin {
    double min, max, xsec, stat;
    std::map<std::string,double> unc;
  };
  std::map<std::string,std::vector<bin>> vars;

  { std::ifstream hepdata(argv[1]);
  auto it = vars.end();
  unsigned line_n = 0;
  for (std::string line, tok; std::getline(hepdata,line); ) {
    ++line_n;
    if (it==vars.end()) {
      if (starts_with(line,"*dataset:")) {
        const auto emp = vars.emplace(std::piecewise_construct,
          std::forward_as_tuple(line.substr(line.rfind('/')+1)),
          tie()
        );
        if (!emp.second) {
          cerr << "repeated variable: " << emp.first->first << endl;
          continue;
        }
        it = emp.first;
      }
    } else {
      const bool star = starts_with(line,"*");
      if (it->second.size()==0 && star) continue;
      if (line.size() && !star) {
        it->second.emplace_back();
        bin& b = it->second.back();

        const auto d1 = line.find(';');
        tok = line.substr(0,d1);
        if (starts_with(tok,">=")) tok.erase(0,2);
        std::stringstream ss(tok);
        ss >> b.min >> tok >> b.max;
        if (tok!="TO") b.max = b.min;

        const auto d2 = line.find('(',d1+1);
        ss.clear();
        ss.str(line.substr(d1+1,d2-d1-1));
        ss >> b.xsec >> tok >> b.stat;
        if (tok!="+-") throw std::runtime_error("missing +- in bin line");

        const char* cstr = line.c_str();
        for (size_t i=d2+1; line[i]!=';'; ) {
          // cout << cstr+i << endl;
          if (!starts_with(cstr+i,"DSYS")) throw std::runtime_error(
            "missing DSYS in bin line");
          const auto eq  = line.find('=',i);
          const auto col = line.find(':',eq+1);
          auto end = line.find(',',col+1);
          if (end==std::string::npos) end = line.find(')',col+1);
          const size_t sep = std::find(cstr+eq+1,cstr+col,',')-cstr;

          const auto unc = line.substr(col+1,end-col-1);
          if (!b.unc.emplace(
            unc,
            sep==col // one value
            ? std::stod(line.substr(eq+1,col-eq-1))
            : std::max( std::abs(std::stod(line.substr(eq +1,sep-eq -1))),
                        std::abs(std::stod(line.substr(sep+1,col-sep-1))) )
          ).second) throw std::runtime_error(cat(
            "duplicate uncert source \'",unc,"\' on line ",line_n));

          i = end+1;
        }

      } else it = vars.end();
    }
  }}

  // ================================================================

  static const std::unordered_map<std::string,std::string> tex {
    {"N_j_30", "N_{jets}"},
    {"N_j_50", "N_{jets}^{ #geq50 GeV}"},
    {"pT_yy", "p_{T}^{#gamma#gamma} [GeV]"},
    {"pTt_yy", "p_{Tt}^{#gamma#gamma} [GeV]"},
    {"pT_yyjj_30", "p_{T}^{#gamma#gammajj} [GeV]"},
    {"HT_30", "H_{T} [GeV]"},
    {"yAbs_yy", "|y_{#gamma#gamma}|"},
    {"yAbs_j1_30", "|y_{j1}|"},
    {"yAbs_j2_30", "|y_{j2}|"},
    {"Dphi_j_j_30", "|#Delta#phi_{jj}|"},
    {"Dphi_j_j_30_signed", "#Delta#phi_{jj}"},
    {"Dphi_yy_jj_30", "|#Delta#phi_{#gamma#gamma,jj}|"},
    {"pT_j1_30", "p_{T}^{j1} [GeV]"},
    {"pT_j2_30", "p_{T}^{j2} [GeV]"},
    {"cosTS_yy", "|cos #theta*|"},
    {"m_jj_30", "m_{jj} [GeV]"},
    {"Dy_j_j_30", "|#Deltay_{jj}|"},
    {"Dy_y_y", "|#Deltay_{#gamma#gamma}|"},
    {"maxTau_yyj_30", "max #tau_{#gamma#gammaj} [GeV]"},
    {"sumTau_yyj_30", "sum #tau_{#gamma#gammaj} [GeV]"},
    {"fid_incl", "Inclusive"},
    {"fid_VBF", "VBF enhanced"},
    {"fid_lep1", "N_{lept} #geq 1"}
  };

  TCanvas canv;
  canv.SetBottomMargin(0.13);
  canv.SetRightMargin(0.035);
  canv.SetTopMargin(0.03);
  if (!burst) canv.SaveAs("uncert.pdf[");

  gPad->SetTickx();
  gPad->SetTicky();

  // bool skip = true;
  for (const auto& var : vars) {
    cout << var.first << endl;

    // collect bin edges
    auto edges = var.second | [](const auto& b){ return b.min; };
    edges.push_back(var.second.back().max);

    // collect uncertainties
    auto uncs = var.second | [](const auto& b){
      return std::vector<double> {
        b.unc.at("lumi"),
        [&b]{
          double x = 0;
          for (const auto& unc : b.unc) {
            if (unc.first=="lumi" || unc.first=="fit") continue;
            x += sq(unc.second);
          }
          return std::sqrt(x);
        }(),
        b.unc.at("fit"),
        b.stat
      };
    };

    // for (auto& unc : uncs) {
    //   for (auto x : unc)
    //     cout << ' ' << x;
    //   cout << endl;
    // }

    // partial sums in quadrature
    for (auto& unc : uncs)
      for (unsigned i=1; i<unc.size(); ++i)
        unc[i] = qadd(unc[i],unc[i-1]);

    // divide by cross section
    tie(uncs,var.second) * [](auto& unc, const auto& b){
      for (auto& u : unc) u /= b.xsec;
    };

    // cout << endl;
    // for (auto& unc : uncs) {
    //   for (auto x : unc)
    //     cout << ' ' << x;
    //   cout << endl;
    // }

    transpose_container_t<decltype(uncs)> tuncs(uncs.front().size());
    for (unsigned i=0; i<uncs.front().size(); ++i) {
      tuncs[i].resize(uncs.size());
      for (unsigned j=0; j<uncs.size(); ++j)
        tuncs[i][j] = uncs[j][i];
    }

    static const std::vector<std::array<int,3>> styles {
      {{kAzure-6,1,1}},
      {{kAzure+8,1,3}},
      {{kAzure-8,1,2}},
      {{17,1,1}}
    };

    const auto bands = tie(tuncs,styles) *
      [&edges](const auto& unc, const auto& style){
        auto band = make_band(edges, unc);
        band->SetFillColor(get<0>(style));
        band->SetLineColor(get<1>(style));
        band->SetLineStyle(get<2>(style));
        auto outline = make_outline(band.get());
        return std::array<h_ptr,3> {
          std::move(band),
          std::move(get<0>(outline)),
          std::move(get<1>(outline))
        };
      };

    const auto& total = bands.back();
    get<0>(total)->SetTitle("");
    TAxis *xa = get<0>(total)->GetXaxis(),
          *ya = get<0>(total)->GetYaxis();
    xa->SetTitle(tex.at(var.first).c_str());
    xa->SetTitleOffset(0.95);
    ya->SetTitleOffset(0.75);
    ya->SetTitle("#Delta#sigma_{fid}/#sigma_{fid}");
    xa->SetTitleSize(0.06);
    xa->SetLabelSize(0.05);
    ya->SetTitleSize(0.065);
    ya->SetLabelSize(0.05);

    int max = std::ceil( std::abs( *std::max_element(
      tuncs.back().begin(), tuncs.back().end() ) )
      + ( tuncs.back().size()>1 ? 0.5 : 0. )
    );
    if (max%2 && max!=1) max += 1;
    max = std::min(max,8);
    ya->SetRangeUser(-max,max);
    get<0>(total)->Draw("E2");

    get<1>(total)->Draw("same");
    get<2>(total)->Draw("same");

    if (starts_with(var.first,"N_j_")) {
      for (unsigned i=0, n=edges.size()-1; i<n; ++i) {
        xa->SetBinLabel( i+1, cat(
          n-i>1 ? " = " : " #geq ", std::ceil(edges[i])
        ).c_str() );
      }
      xa->SetLabelSize(0.08);
    } else if (var.first.substr(0,4)=="fid_") {
      xa->SetBinLabel(1,"");
    }

    // draw in oposite order, so that smaller values can be seen
    for (unsigned i=bands.size()-1; i; ) {
      --i;
      for (unsigned j=0, n=bands[i].size(); j<n; ++j)
        bands[i][j]->Draw("E2same" + (j ? 2 : 0));
    }

    gPad->RedrawAxis();

    static const std::vector<const char*> labels {
      "Luminosity",
      "#oplus Correction factor",
      "#oplus Signal extraction",
      "#oplus Statistics"
    };

    TLegend leg(0.12,0.1525,0.72,0.2525);
    leg.SetLineWidth(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.041);
    leg.SetNColumns(2);
    tie(bands,labels) * [&leg](const auto& band, const char* lbl){
      leg.AddEntry(get<0>(band).get(),lbl,"f");
    };
    leg.Draw();

    TLatex l;
    l.SetTextColor(1);
    l.SetNDC();
    l.SetTextFont(72);
    l.DrawLatex(0.135,0.83,"ATLAS");
    l.SetTextFont(42);
    l.DrawLatex(0.255,0.83,"Internal");
    // l.DrawLatex(0.255,0.83,"Preliminary");
    l.SetTextFont(42);
    l.DrawLatex(0.135,0.89,
      "#it{H} #rightarrow #gamma#gamma, "
      "#sqrt{#it{s}} = 13 TeV, 36.1 fb^{-1}, "
      "m_{H} = 125.09 GeV"
    );
    l.SetTextFont(42);

    canv.SaveAs(burst ? (var.first+".pdf").c_str() : "uncert.pdf");
  }
  if (!burst) canv.SaveAs("uncert.pdf]");
}
