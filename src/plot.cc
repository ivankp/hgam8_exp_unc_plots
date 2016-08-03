#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <vector>
#include <array>
#include <unordered_map>
#include <cmath>
#include <cstring>
#include <memory>

#include <TCanvas.h>
#include <TAxis.h>
#include <TColor.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define burst

using namespace std;

using h_t = TH1F;
using h_ptr = std::unique_ptr<h_t>;

h_ptr make_band(
  const vector<double>& bins,
  const vector<double>& height
) {
  const unsigned nbins = bins.size();
  vector<double> bins_width(nbins), cent(nbins,0);

  for (unsigned i=0; i<nbins-1; ++i)
    bins_width[i] = bins[i+1]-bins[i];

  h_t *h = new h_t("","",nbins-1,bins.data());
  for (unsigned i=0, n=height.size(); i<n; ++i) {
    h->SetBinError(i+1,height[i]);
  }
  h->SetStats(0);
  h->SetMarkerStyle(0);
  h->SetLineWidth(1); // gives legend color boxes outlines

  return h_ptr(h);
}

array<h_ptr,2> make_outline(h_t* h) {
  auto* xa = h->GetXaxis();
  const unsigned nbins = h->GetNbinsX();
  array<h_ptr,2> hh {
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

constexpr unsigned ncol = 5;

int main(int argc, char const *argv[]) {
  if (argc>3 || (argc==2 && !strcmp(argv[1],"-h"))) {
    cout << "usage: " << argv[0] << " [input_file]" << endl;
    cout << "Default input_file is bands.dat" << endl;
    return 1;
  }

  vector<pair<string,array<vector<double>,ncol>>> vars;

  bool hit_var = false;
  unsigned col = 0;

  ifstream dat(argc==1 ? "bands.dat" : argv[1]);
  if (dat.is_open()) {
    string str;
    while ( dat >> str ) {
      if (str=="var") {
        hit_var = true;
      } else if (hit_var) {
        vars.emplace_back(move(str), array<vector<double>,ncol>());
        hit_var = false;
        col = 0;
      } else if (col<ncol) {
        vars.back().second[col++].emplace_back(stold(str));
        if (col==ncol) col = 0;
      }
    }
    dat.close();
  } else {
    cout << "Unable to open file " << argv[1] << endl;
    return 1;
  }

  unordered_map<string,string> tex {
    {"Njets", "N_{jets}"},
    {"pT_yy", "p_{T}^{#gamma#gamma} [GeV]"},
    {"yAbs_yy", "|y_{#gamma#gamma}|"},
    {"dphi_jj", "|#Delta#phi_{jj}|"},
    {"pT_j1", "p_{T}^{j1} [GeV]"},
    {"Njets50", "N_{jets}^{ #geq50 GeV}"},
    {"cosTS", "|cos #theta*|"},
    {"m_jj", "m_{jj} [GeV]"},
    {"dy_jj", "|y_{jj}|"},
    {"fid_incl", "Inclusive"},
    {"fid_VBF", "VBF enhanced"},
    {"fid_lep1", "N_{lept} #geq 1"}
  };

  TCanvas canv;
  canv.SetBottomMargin(0.13);
  canv.SetRightMargin(0.035);
  canv.SetTopMargin(0.03);
  #ifndef burst
    canv.SaveAs("uncert.pdf[");
  #endif

  gPad->SetTickx();
  gPad->SetTicky();

  for (const auto& var : vars) {
    cout << var.first << endl;

    auto total = make_band(get<0>(var.second),get<1>(var.second));
    total->SetFillColor(17);
    total->SetTitle("");
    TAxis *xa = total->GetXaxis(),
          *ya = total->GetYaxis();
    xa->SetTitle(tex[var.first].c_str());
    xa->SetTitleOffset(0.95);
    ya->SetTitleOffset(0.75);
    ya->SetTitle("#Delta#sigma_{fid}/#sigma_{fid}");
    xa->SetTitleSize(0.06);
    xa->SetLabelSize(0.05);
    ya->SetTitleSize(0.065);
    ya->SetLabelSize(0.05);

    int max = ceil( abs( *max_element(
      get<1>(var.second).begin(),get<1>(var.second).end() ) )
      + ( total->GetNbinsX()>1 ? 0.5 : 0. )
    );
    if (max%2 && max!=1) max += 1;
    max = std::min(max,8);
    ya->SetRangeUser(-max,max);
    total->SetLineColor(1);
    total->SetLineStyle(1);
    total->Draw("E2");

    auto total_outline = make_outline(total.get());
    get<0>(total_outline)->Draw("same");
    get<1>(total_outline)->Draw("same");

    if (var.first.substr(0,5)=="Njets") {
      for (unsigned i=0, n=get<0>(var.second).size()-1; ; ++i) {
        stringstream ss;
        if (n-i>1) {
          ss << " = " << ceil(get<0>(var.second)[i]);
          xa->SetBinLabel(i+1,ss.str().c_str());
        } else {
          ss << " #geq " << ceil(get<0>(var.second)[i]);
          xa->SetBinLabel(i+1,ss.str().c_str());
          break;
        }
      }
      xa->SetLabelSize(0.08);
    } else if (var.first.substr(0,4)=="fid_") {
      xa->SetBinLabel(1,"");
    }

    auto lce = make_band(get<0>(var.second),get<2>(var.second));
    lce->SetFillColor(kAzure-8);
    lce->SetLineColor(1);
    lce->SetLineStyle(2);
    lce->Draw("sameE2");

    auto lce_outline = make_outline(lce.get());
    get<0>(lce_outline)->Draw("same");
    get<1>(lce_outline)->Draw("same");

    auto lc = make_band(get<0>(var.second),get<3>(var.second));
    lc->SetFillColor(kAzure+8);
    lc->SetLineColor(1);
    lc->SetLineStyle(3);
    lc->Draw("sameE2");

    auto lc_outline = make_outline(lc.get());
    get<0>(lc_outline)->Draw("same");
    get<1>(lc_outline)->Draw("same");

    auto lumi = make_band(get<0>(var.second),get<4>(var.second));
    lumi->SetFillColor(kAzure-6);
    lumi->SetLineColor(1);
    lumi->SetLineStyle(1);
    lumi->Draw("sameE2");

    auto lumi_outline = make_outline(lumi.get());
    get<0>(lumi_outline)->Draw("same");
    get<1>(lumi_outline)->Draw("same");

    gPad->RedrawAxis();

    // TLegend leg(0.58,0.725,0.96,0.925);
    TLegend leg(0.12,0.1525,0.72,0.2525);
    leg.SetLineWidth(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.041);
    leg.SetNColumns(2);
    leg.AddEntry(lumi .get(),"Luminosity","f");
    leg.AddEntry(lc   .get(),"#oplus Correction factor","f");
    leg.AddEntry(lce  .get(),"#oplus Signal extraction","f");
    leg.AddEntry(total.get(),"#oplus Statistics","f");
    leg.Draw();

    TLatex l;
    l.SetTextColor(1);
    l.SetNDC();
    l.SetTextFont(72);
    l.DrawLatex(0.135,0.83,"ATLAS");
    l.SetTextFont(42);
    l.DrawLatex(0.255,0.83,"Internal");
    l.SetTextFont(42);
    l.DrawLatex(0.135,0.89,
      "#it{H} #rightarrow #gamma#gamma, "
      "#sqrt{#it{s}} = 13 TeV, 13.3 fb^{-1}, "
      "m_{H} = 125.09 GeV"
    );
    l.SetTextFont(42);
    // l.DrawLatex(0.255,0.83,"m_{H} = 125.09 GeV");

    #ifdef burst
      canv.SaveAs((var.first+".pdf").c_str());
    #else
      canv.SaveAs("uncert.pdf");
    #endif
  }
  #ifndef burst
    canv.SaveAs("uncert.pdf]");
  #endif

  return 0;
}
