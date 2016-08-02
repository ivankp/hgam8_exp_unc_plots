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
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using namespace std;

std::unique_ptr<TGraphAsymmErrors> make_band(
  const vector<double>& bins,
  const vector<double>& height
) {
  const unsigned nbins = bins.size();
  vector<double> bins_width(nbins), cent(nbins,0);

  for (unsigned i=0; i<nbins-1; ++i)
    bins_width[i] = bins[i+1]-bins[i];

  TH1F *h = new TH1F("","",nbins-1,bins.data());
  h->SetStats(0);

  auto *g = new TGraphAsymmErrors(nbins,bins.data(),cent.data(),
                                  0,bins_width.data(),
                                  height.data(),height.data());

  g->SetHistogram(h);
  g->GetXaxis()->SetRangeUser(bins[0],bins.back()+bins_width.back());
  g->SetLineWidth(1); // gives legend color boxes outlines

  return std::unique_ptr<TGraphAsymmErrors>(g);
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
    {"Njets50", "N_{jets}^{#geq50 GeV}"},
    {"cosTS", "|cos #theta*|"},
    {"m_jj", "m_{jj} [GeV]"},
    {"dy_jj", "|y_{jj}|"},
    {"fid_incl", "Inclusive"},
    {"fid_VBF", "VBF enhanced"},
    {"fid_lep", "N_{lept} #geq 1"}
  };

  TCanvas canv;
  canv.SetBottomMargin(0.12);
  canv.SetRightMargin(0.03);
  canv.SetTopMargin(0.03);
  canv.SaveAs("uncert.pdf[");

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
    xa->SetTitleOffset(1.05);
    ya->SetTitleOffset(1.);
    ya->SetTitle(
      "Fractional uncertainty on cross section, #Delta#sigma_{fid}/#sigma_{fid}");
    xa->SetTitleSize(0.05);
    xa->SetLabelSize(0.04);
    ya->SetTitleSize(0.045);
    ya->SetLabelSize(0.04);

    double max = ceil( abs( *max_element(
      get<1>(var.second).begin(),get<1>(var.second).end() ) ) );
    if (unsigned(max)%2) max += 1;
    max = std::min(max,8.);
    ya->SetRangeUser(-max,max);
    total->Draw("a2");

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
      xa->SetLabelSize(0.06);
    } else if (var.first.substr(0,4)=="fid_") {
      xa->SetBinLabel(1,"");
    }

    auto lce = make_band(get<0>(var.second),get<2>(var.second));
    lce->SetFillColor(kAzure-8);
    // lce->SetLineStyle(2);
    lce->Draw("2");

    auto lc = make_band(get<0>(var.second),get<3>(var.second));
    lc->SetFillColor(kAzure+8);
    // lc->SetLineStyle(3);
    lc->Draw("2");

    auto lumi = make_band(get<0>(var.second),get<4>(var.second));
    lumi->SetFillColor(kAzure-6);
    // lumi->SetLineStyle(4);
    lumi->Draw("2");

    gPad->RedrawAxis();

    TLegend leg(0.58,0.725,0.96,0.925);
    leg.SetLineWidth(0);
    leg.SetFillColorAlpha(0,0);
    leg.AddEntry(lumi .get(),"Luminosity","f");
    leg.AddEntry(lc   .get(),"#oplus Correction factor syst.","f");
    leg.AddEntry(lce  .get(),"#oplus Signal extraction syst.","f");
    leg.AddEntry(total.get(),"#oplus Statistics","f");
    leg.Draw();

    TLatex l;
    l.SetTextColor(1);
    l.SetNDC();
    l.SetTextFont(72);
    l.DrawLatex(0.15,0.87,"ATLAS");
    l.SetTextFont(42);
    l.DrawLatex(0.27,0.87,"Internal");
    l.SetTextFont(42);
    l.DrawLatex(0.15,0.80,
      "#it{H} #rightarrow #gamma#gamma, "
      "#sqrt{#it{s}} = 13 TeV, 13.3 fb^{-1}"
    );

    canv.SaveAs("uncert.pdf");
  }
  canv.SaveAs("uncert.pdf]");

  return 0;
}
