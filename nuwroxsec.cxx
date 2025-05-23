#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <algorithm>
#include <memory>

int main(int argc, char **argv) {
  ROOT::EnableImplicitMT();
  TH1::AddDirectory(false);
  std::vector<std::string> files{};
  std::for_each(argv + 1, argv + argc,
                [&files](const char *arg) { files.push_back(arg); });
  auto df = ROOT::RDataFrame{"nRooTracker", files}.Define(
      "E", [](const ROOT::RVecD &StdHepP4) { return StdHepP4[3]; },
      {"StdHepP4"});
  ROOT::RDF::Experimental::AddProgressBar(df);
  auto xsec = df.Mean("EvtWght");
  auto count = df.Count();
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hists{};
  std::vector<double> bin_profile{};
  double logEmax = TMath::Log10(1.);
  double logEmin = TMath::Log10(0.1);
  int nstep = 5;
  double dlogE = (logEmax - logEmin) / (double)nstep;
  for (int i = 0; i <= nstep * 3; i++) {
    bin_profile.push_back(TMath::Power(10, logEmin + i * dlogE));
  }
  assert(bin_profile.last() == 100.);
  const int nbins = bin_profile.size() - 1;
  hists.emplace_back(
      df.Histo1D({"tot_cc", "E;E;xsec", nbins, bin_profile.data()}, "E"));
  hists.emplace_back(
      df.Filter(
            [](const TObjString &EvtCode) {
              return EvtCode.GetString() == "1";
            },
            {"EvtCode"})
          .Histo1D({"qel_cc", "E;E;xsec", nbins, bin_profile.data()}, "E"));
  hists.emplace_back(
      df.Filter(
            [](const TObjString &EvtCode) {
              return EvtCode.GetString() == "11";
            },
            {"EvtCode"})
          .Histo1D({"res_cc", "E;E;xsec", nbins, bin_profile.data()}, "E"));
  hists.emplace_back(
      df.Filter(
            [](const TObjString &EvtCode) {
              return EvtCode.GetString() == "26";
            },
            {"EvtCode"})
          .Histo1D({"dis_cc", "E;E;xsec", nbins, bin_profile.data()}, "E"));
  hists.emplace_back(
      df.Filter(
            [](const TObjString &EvtCode) {
              return EvtCode.GetString() == "2";
            },
            {"EvtCode"})
          .Histo1D({"mec_cc", "E;E;xsec", nbins, bin_profile.data()}, "E"));

  auto file = std::make_unique<TFile>("output.root", "RECREATE");

  std::ranges::for_each(hists, [&](auto &&hist) {
    hist->Scale(xsec.GetValue() / count.GetValue() * (100. - 0.1), "WIDTH");
    file->Add(hist.GetPtr());
  });
  // hist->Scale(xsec.GetValue() / count.GetValue());
  file->Write();
  file->Close();
  // hist->SaveAs("output.root");
  return 0;
}