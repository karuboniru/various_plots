#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <THStack.h>
#include <TLegend.h>
#include <array>
#include <memory>
#include <plot.h>
using std::string_literals::operator""s;

const size_t nbins = 25;

const std::array<std::string, 3> variables{"Wtrue", "xBj", "Q2"},
    axis_titles{"W (GeV)", "x_{Bj}", "Q^{2} (GeV^{2})"};

const std::array<std::tuple<double, double>, 3> edges{
    std::tuple<double, double>{0.9, 3.5}, {0, 1.}, {0., 2.}};

const std::array<std::tuple<std::string, std::string, std::string>, 2>
    name_and_path{
        std::tuple<std::string, std::string, std::string>{
            "pi0", "CC#pi 0",
            "anaNuWro/outplot/"
            "outAna7_MINERvANuWro_test_GFSPIZEROa7nuCH_Filelist_NuWro_test_"
            "BeamNu_EnuMINERvA_TargetCarbon.root"},
        std::tuple<std::string, std::string, std::string>{
            "0pi", "CC0#pi",
            "anaNuWro/outplot/"
            "outAna9_MINERvANuWro_test_GFS0PIa9nuCH_Filelist_NuWro_test_"
            "BeamNu_EnuMINERvA_TargetCarbon.root"}};

void saveplot(auto &&h2) {
  auto c1 = getCanvas();
  ResetStyle(h2);
  h2->Draw("colz");
  c1->SaveAs((h2->GetName() + ".pdf"s).c_str());
  c1->SaveAs((h2->GetName() + ".eps"s).c_str());
  c1->SaveAs((h2->GetName() + ".png"s).c_str());
}

int main() {
  ROOT::EnableImplicitMT();
  for (auto &&[name, title, path] : name_and_path) {
    // ROOT::RDataFrame df("tree", path);
    ROOT::RDataFrame dfr("tree", path);
    auto count = dfr.Count();
    auto xsec = dfr.Mean("xsec");
    auto normalize_factor = xsec.GetValue() / count.GetValue();
    auto df = dfr.Filter([](double xBj) { return xBj < 1.; }, {"xBj"});
    for (size_t i{}; i < variables.size(); ++i) {
      auto variable_x = variables[i];
      auto title_x = axis_titles[i];
      auto [x_min, x_max] = edges[i];
      for (size_t j{}; j < variables.size(); ++j) {
        if (i <= j)
          continue;
        auto variable_y = variables[j];
        auto title_y = axis_titles[j];
        auto [y_min, y_max] = edges[j];
        auto h2 = df.Histo2D<double, double>(
            {(name + variable_x + variable_y).c_str(), "title", nbins, x_min,
             x_max, nbins, y_min, y_max},
            variable_x, variable_y);
        h2->GetXaxis()->SetTitle(title_x.c_str());
        h2->GetYaxis()->SetTitle(title_y.c_str());
        h2->Scale(normalize_factor, "width");
        saveplot(h2);
        saveplot(normalize_slice(h2, true));
        saveplot(normalize_slice(h2, false));
      }
    }
  }

  return 0;
}