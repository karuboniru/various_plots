#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <THStack.h>
#include <TLegend.h>
#include <memory>
#include <plot.h>
#include <vector>

int main(int argc, char **argv) {
  const Int_t cols[] = {kBlack,    kRed,         kGray,       kMagenta,
                        kOrange,   kBlue,        kGreen,      kRed - 9,
                        kGray + 3, kMagenta + 3, kOrange + 3, kBlue + 3,
                        kGreen + 3};
  const int *current_col = cols;

  // list of parameters
  // v: variable to draw
  std::string var;
  // f: file to draw from
  std::vector<std::string> file;
  // c: overall cut
  std::string cut;
  // s: condition to split the data
  std::vector<std::string> split;
  // o: output file name
  std::string output;
  for (char t; (t = getopt(argc, argv, "v:f:c:s:o:h")) != -1;) {
    switch (t) {
    case 'v':
      var = optarg;
      break;
    case 'f':
      file.push_back(optarg);
      break;
    case 'c':
      cut = optarg;
      break;
    case 's':
      split.push_back(optarg);
      break;
    case 'o':
      output = optarg;
      break;
    case 'h':
    default:
      std::cerr << "Usage: " << argv[0]
                << " -v var -f file ... -c cut -s split ... -o output"
                << std::endl;
      return 1;
    }
  }
  if (var.empty() || file.empty() || cut.empty() || split.empty() ||
      output.empty()) {
    std::cerr << "Missing arguments" << std::endl;
    return 1;
  }
  auto file_sample = std::make_unique<TFile>(
      "/media/storage/neutrino/nuwro_data/mhtest/minerva_ghent_alter_blending/"
      "TKI/anaNuWro/outStack/MINERvANuWroGFSPIZEROa7t4nuCarbon/"
      "MINERvANuWroGFSPIZEROa7t4nuCarbon.root",
      "READ");
  auto hist_sample = file_sample->Get<TH1D>(
      ("MINERvANuWroGFSPIZEROa7t4nuCarbon/" + var + "all").c_str());

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d("tree", file);
  auto d_cut = d.Filter(cut);
  std::vector<ROOT::RDF::RResultPtr<TH1>> plots_to_stack{};
  std::vector<std::unique_ptr<TH1>> fraction_plots{};
  auto h_all = hist_sample ? d_cut.Histo1D<double>(*hist_sample, var)
                           : d_cut.Histo1D<double>(var);

  auto stack = std::make_unique<THStack>("hs", "");
  auto legend = std::make_unique<TLegend>(0.6, 0.6, 0.9, 0.9);
  ResetStyle(legend);

  for (auto &s : split) {
    auto d_split = d_cut.Filter(s);
    auto h = d_split.Histo1D(*h_all, var);
    h->SetLineColor(*current_col);
    h->SetLineWidth(2);
    fraction_plots
        .emplace_back(static_cast<TH1 *>(
            h->Clone((h->GetName() + std::string("_fraction")).c_str())))
        ->Divide(h_all.GetPtr());
    // h->SetFillColor(*current_col);
    current_col++;
    plots_to_stack.push_back(h);
    stack->Add(h.GetPtr());
    legend->AddEntry(h.GetPtr(), s.c_str(), "l");
  }
  {
    auto c = getCanvas();
    stack->Draw("hist");
    stack->GetXaxis()->SetTitle(var.c_str());
    stack->GetYaxis()->SetTitle("Events");
    legend->Draw();
    c->SaveAs(output.c_str());
  }
  {
    auto c = getCanvas();
    stack->Draw("hist NOSTACK");
    stack->GetXaxis()->SetTitle(var.c_str());
    stack->GetYaxis()->SetTitle("Events");
    legend->Draw();
    auto output1 = output.substr(0, output.find_last_of('.')) + "_nostack.pdf";
    c->SaveAs(output1.c_str());
  }
  {
    auto c = getCanvas();
    bool first{true};
    for (auto &h : fraction_plots) {
      h->GetXaxis()->SetTitle(var.c_str());
      h->GetYaxis()->SetTitle("Fraction");
      h->SetMinimum(0.);
      h->SetMaximum(1.);
      h->Draw(first ? "HIST" : "HIST SAME");
      first = false;
    }
    legend->Draw();
    c->SaveAs(
        (output.substr(0, output.find_last_of('.')) + "_fraction.pdf").c_str());
  }
  {
    auto c = getCanvas();
    auto hs = std::make_unique<THStack>("stacked_fraction", "");
    // hs->GetXaxis()->SetTitle(var.c_str());
    // hs->GetYaxis()->SetTitle("Fraction");
    for (auto &h : fraction_plots) {
      h->GetXaxis()->SetTitle(var.c_str());
      h->GetYaxis()->SetTitle("Fraction");
      h->SetFillColor(h->GetLineColor());
      hs->Add((TH1 *)h.get()->Clone());
    }
    // hs->Draw();
    // hs->GetXaxis()->SetTitle(var.c_str());
    // hs->GetYaxis()->SetTitle("Fraction");
    hs->SetTitle((";" + var + ";Fraction").c_str());
    hs->SetMaximum(1.);
    hs->SetMinimum(0.);
    hs->Draw();
    legend->Draw();
    c->SaveAs(
        (output.substr(0, output.find_last_of('.')) + "_fraction_stacked.pdf")
            .c_str());
  }
  return 0;
}