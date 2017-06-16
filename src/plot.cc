#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
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
#include "lists.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

template <typename... T> struct bad_type;

using std::cout;
using std::cerr;
using std::endl;
using std::get;
using std::tie;
using namespace ivanp;
using namespace ivanp::math;

template <typename T>
struct less_deref {
  constexpr bool operator()(const T* lhs, const T* rhs) const {
    return (*lhs) < (*rhs);
  }
};

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
        bool geq;
        if ((geq = starts_with(tok,">="))) tok.erase(0,2);
        std::stringstream ss(tok);
        ss >> b.min >> tok >> b.max;
        if (geq) b.max = b.min+1;
        else if (tok!="TO") b.max = b.min;

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
  if (!burst) canv.SaveAs(cat("uncert",corr  ? "_corr" : "",".pdf[").c_str());

  gPad->SetTickx();
  gPad->SetTicky();

  for (const auto& var : vars) {
    cout << var.first << endl;

    std::vector<const std::string*> corr_selected, corr_other;
    if (corr) { // select most significant contributions
      std::vector<std::pair<const std::string*,std::vector<double>>> corr_uncs;
      for (const auto& bin : var.second) {
        for (const auto& unc : bin.unc) {
          const auto& name = unc.first;
          if (name=="lumi" || name=="fit") continue;
          auto it = std::find_if(corr_uncs.begin(),corr_uncs.end(),
            [&name](const auto& x){ return name == *x.first; });
          if (it==corr_uncs.end()) {
            corr_uncs.emplace_back(&name,decltype(it->second){});
            it = --corr_uncs.end();
          }
          it->second.push_back(unc.second/bin.xsec);
        }
      }
      auto corr_uncs_sorted = corr_uncs | [](const auto& x){
        // sum squares of relative unc in each bin
        return std::make_pair(x.first,
          std::accumulate(x.second.begin(),x.second.end(),0.,
            [](auto total, auto next){ return total + sq(next); }));
      };
      std::sort(corr_uncs_sorted.begin(),corr_uncs_sorted.end(),
        [](const auto& a, const auto& b){ return a.second > b.second; });

      const unsigned n = std::min(4u,(unsigned)corr_uncs_sorted.size());
      corr_selected.reserve(n);
      for (unsigned i=n; i; )
        --i, corr_selected.emplace_back(corr_uncs_sorted[i].first);
      if (corr_uncs_sorted.size()>n) {
        corr_other.reserve(corr_uncs_sorted.size()-n);
        for (unsigned i=n; i<corr_uncs_sorted.size(); ++i)
          corr_other.emplace_back(corr_uncs_sorted[i].first);
      }

      // for (const auto* x : corr_selected)
      //   cout <<"  "<< *x << endl;
    }

    // collect bin edges
    const auto edges = ( var.second | [](const auto& b){ return b.min; } )
                     << var.second.back().max;

    // collect uncertainties
    auto uncs = var.second | [&](const auto& b){
      if (!corr) {
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
      } else {
        cout << b.min << endl;
        return ( corr_selected | [&b](const auto& s){
          const auto err = b.unc.at(*s);
          cout <<"  "<< *s <<' '<< err << endl;
          return err;
        } ) << std::sqrt(std::accumulate(
          corr_other.begin(),corr_other.end(),0.,
          [&b](auto total, const auto* s){ return total + sq(b.unc.at(*s)); }
        ));
      }
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

    // transpose_container_t<decltype(uncs)> tuncs(uncs.front().size());
    std::vector<std::vector<double>> tuncs(uncs.front().size());
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
    static const std::vector<std::array<int,3>> styles_corr {
      {{kOrange+9,1,1}},
      {{kOrange+3,1,2}},
      {{kOrange+8,1,1}},
      {{kOrange-2,1,2}},
      {{kOrange-9,1,1}}
    };

    const auto bands = tie(tuncs, corr ? styles_corr : styles) *
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

    const auto max = *std::max_element(tuncs.back().begin(),tuncs.back().end());
    auto range = std::exp2( std::ceil( std::log2(max) ) );
    if (range > 8) range = 8;
    else if (max/range > 0.7) range *= 2;
    ya->SetRangeUser(-range,range);
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

    static const std::vector<std::string> labels {
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
    tie(bands,
        !corr ? labels : (corr_selected | [i=0](auto* s) mutable {
          return (i++ ? "#oplus " : "") + *s;
        }) << "#oplus others"
      ) * [&leg](const auto& band, const std::string& lbl){
        leg.AddEntry(get<0>(band).get(),lbl.c_str(),"f");
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

    canv.SaveAs(cat(
      burst ? var.first : "uncert",
      corr  ? "_corr" : "",
      ".pdf").c_str());
  }
  if (!burst) canv.SaveAs(cat("uncert",corr  ? "_corr" : "",".pdf]").c_str());
}
