#pragma once
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TObjString.h>
#include <TSpline.h>
#include <event.h>
#include <memory>
#include <plot.h>
#include <string>
#include <vector>

void make_pie_plot(std::vector<std::pair<std::string, double>> data,
                   std::string filename);

template <typename T, typename U> void save(T &&t, U &&u) {
  for (auto &&tt : t) {
    tt->SetDirectory(&*u);
    ResetStyle(tt);
    auto c1 = getCanvas();
    tt->Draw("HIST");
    tt->SetLineWidth(2);
    // c1->SaveAs((tt->GetName() + std::string(".pdf")).c_str());
    // c1->SaveAs((tt->GetName() + std::string(".png")).c_str());
    // c1->SaveAs((tt->GetName() + std::string(".eps")).c_str());
    // c1->SetLogy();
    // c1->SaveAs((tt->GetName() + std::string("_log.pdf")).c_str());
    // c1->SaveAs((tt->GetName() + std::string("_log.png")).c_str());
    // c1->SaveAs((tt->GetName() + std::string("_log.eps")).c_str());
  }
  u->Write();
  for (auto &&tt : t) {
    tt->SetDirectory(0);
  }
}

template <typename T> void draw(T &&tt, const char *title = nullptr, const char *opt = "HIST") {
  ResetStyle(tt);
  auto c1 = getCanvas();
  if (title) {
    tt->SetTitle(title);
  }
  tt->Draw(opt);
  tt->SetLineWidth(2);
  c1->SaveAs((tt->GetName() + std::string(".pdf")).c_str());
  c1->SaveAs((tt->GetName() + std::string(".png")).c_str());
  c1->SaveAs((tt->GetName() + std::string(".eps")).c_str());
  c1->SetLogy();
  c1->SaveAs((tt->GetName() + std::string("_log.pdf")).c_str());
  c1->SaveAs((tt->GetName() + std::string("_log.png")).c_str());
  c1->SaveAs((tt->GetName() + std::string("_log.eps")).c_str());
}

template <typename T, typename U>
void draw_same(T &&x1, U &&x2, std::string title = "", double chi2 = 0.) {
  TLegend leg{0.55, 0.7, 0.85, 0.9};
  ResetStyle(&leg);
  auto max = x2->GetMaximum();
  x1->SetMaximum(max * 1.8);
  // x2->SetMaximum(max * 1.8);
  x1->SetMinimum(0);
  x2->SetMinimum(0);
  ResetStyle(x1);
  ResetStyle(x2);
  x1->SetLineColor(kRed);
  x2->SetLineColor(kBlue);
  x1->SetTitle(title.c_str());
  x2->SetTitle(title.c_str());
  if (chi2) {
    leg.SetHeader(("#Chi^{2} = " + std::to_string(chi2)).c_str());
  }
  leg.AddEntry(&*x1, "NuWro", "l");
  leg.AddEntry(&*x2, "Data", "l");
  auto c1 = getCanvas();
  x1->Draw("hist");
  x2->Draw("E0 same");
  leg.Draw();
  c1->SaveAs((x1->GetName() + std::string(".pdf")).c_str());
  c1->SaveAs((x1->GetName() + std::string(".png")).c_str());
  c1->SaveAs((x1->GetName() + std::string(".eps")).c_str());
}

event::channel getmode_nuwro(TObjString &code);

event::channel get_mode_genie(const TObjString &code);

template <typename T>
std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
  TFile root_file{file_path.c_str(), "READ"};
  auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
  assert(objptr);
  return std::unique_ptr<T>{objptr};
}

std::pair<double, double> get_xsec(TH1 *h_rate, TGraph *spline);

void BinLogX(TAxis* axis);