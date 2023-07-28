#include "TH1.h"
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TObjString.h>
#include <TPie.h>
#include <array>
#include <event.h>
#include <memory>
#include <string>
#include <tools.h>
#include <vector>

void make_pie_plot(std::vector<std::pair<std::string, double>> data,
                   std::string filename) {
  constexpr std::array<int, 10> col{kRed,   kBlue, kViolet, kYellow, kOrange,
                                    kGreen, kGray, kTeal,   kPink};
  auto pie = std::make_unique<TPie>("final state", "final state", data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    pie->SetEntryVal(i, data[i].second);
    pie->SetEntryFillColor(i, col[i % col.size()]);
    pie->SetEntryFillStyle(i, 1000 + i / col.size());
    pie->SetEntryLabel(i, data[i].first.c_str());
  }
  auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 700, 700);
  canvas->cd();
  // pie->SetRadius(0.25);
  pie->SetCircle(0.5, 0.5 - .1, 0.35);
  auto leg = pie->MakeLegend(.6, .6, .9, .9);
  pie->SetLabelFormat("%perc");
  pie->SetLabelsOffset(-.2);
  pie->Draw("");
  leg->Draw("SAME");
  canvas->SaveAs(filename.c_str());
  return;
}

event::channel getmode_nuwro(TObjString &code) {
  const int mode = code.GetString().Atoi();
  // cout << "mode = " << mode << endl;
  switch (mode) {
  case 1:
    return event::channel::QE;
    break;
  case 11:
    return event::channel::RES;
    break;
  case 26:
    return event::channel::DIS;
    break;
  case 2:
    return event::channel::MEC;
    break;
  case 16:
  default:
    return event::channel::Other;
    break;
  }
}

event::channel get_mode_genie(const TObjString &code) {
  if (code.GetString().Contains("QES")) {
    return event::channel::QE;
  } else if (code.GetString().Contains("RES")) {
    return event::channel::RES;
  } else if (code.GetString().Contains("DIS")) {
    return event::channel::DIS;
  } else if (code.GetString().Contains("MEC")) {
    return event::channel::MEC;
  }
  return event::channel::Other;
}

void normalize(TH1 *h, double xsec, double total) {
  h->Scale(1 / total * xsec, "width");
}

std::pair<double, double> get_xsec(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  // spline->SaveAs("wrong.root");
  TSpline3 sp("sp", spline);
  TF1 func(
      "spline", [&](double *x, double *) { return sp.Eval(*x); }, 0,
      h_rate->GetXaxis()->GetXmax(), 0);
  for (int ii = 1; ii <= h_rate->GetNbinsX(); ii++) {
    double bin_c = h_rate->GetBinContent(ii);
    double bin_up = h_rate->GetXaxis()->GetBinUpEdge(ii);
    double bin_low = h_rate->GetXaxis()->GetBinLowEdge(ii);
    double bin_width = bin_up - bin_low;
    if (bin_c < 1 || func.Integral(bin_low, bin_up) == 0) {
      continue;
    }
    fluxint += bin_c / func.Integral(bin_low, bin_up) * bin_width;
  }
  double event_rate = h_rate->Integral();
  return {event_rate, event_rate / fluxint};
}

void BinLogX(TAxis *axis) {
  // void XGLUtils::BinLogX(TAxis *axis)
  //
  //  Method for the correct logarithmic binning of histograms
  //  copied and modified from AliTPCcalibBase

  const Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  if (from < 1e-2)
    from = 1e-2;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  const Double_t factor = TMath::Power(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i - 1];
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}