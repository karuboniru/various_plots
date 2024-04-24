#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <THStack.h>
#include <TLegend.h>
#include <array>
#include <cmath>
#include <memory>
#include <plot.h>

auto plotW(const char *name, auto &&df_in) {
  return df_in.template Histo1D<double>({name, name, 360, 0, 3.4}, "Wtrue",
                                        "xsec");
}

auto plotpnratio(const char *name, auto &&df_in) {
  return df_in.template Histo1D<double>({name, name, 360, -1, 5.}, "ratio1",
                                        "xsec");
}

auto plotpnratio2(const char *name, auto &&df_in) {
  return df_in.template Histo1D<double>({name, name, 800, -1., 2.}, "ratio2",
                                        "xsec");
}

auto plotangle(const char *name, auto &&df_in) {
  return df_in.template Histo1D<double>({name, name, 800, 0, M_1_PI / 2.},
                                        "thetaL", "xsec");
}

auto plotdiff(const char *name, auto &&df_in) {
  return df_in.template Histo1D<double>({name, name, 360, -1., 5.}, "diff",
                                        "xsec");
}

void scale_plot(double EventRate, TH1D *ptr) {
  ptr->Scale(1e-38 / EventRate, "WIDTH");
}

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <input file> <output file> [snapshot_out]" << std::endl;
    exit(1);
  }
  TH1::AddDirectory(false);
  std::string input_file = argv[1];
  std::string output_file = argv[2];

  ROOT::EnableImplicitMT();
  auto fin = std::unique_ptr<TFile>(TFile::Open(input_file.c_str()));
  double EventRate = ((TH1F *)fin->Get("hCCrate"))->GetSumOfWeights();
  auto df = ROOT::RDataFrame{"tree", input_file}
                .Define("dpl", "TMath::Sqrt(IApN*IApN - dpt * dpt)")
                .Define("ratio1", "IApN/dpt - 1")
                .Define("ratio2", "dpl/IApN")
                .Define("thetaL", "TMath::ASin(ratio2)")
                .Define("diff", "IApN-dpt");
  // auto nevent = df.Count();
  auto hw = plotW("Wtrue_all", df);
  auto hw_hydrogen = plotW("Wtrue_hydrogen", df.Filter("targetZ==1"));
  auto hw_hydrogen_res = plotW(
      "Wtrue_hydrogen_res", df.Filter("targetZ==1 && (evtMode == 2 && "
                                      "(flag_delta == 1 || Wtrue <= 1.210))"));
  auto hw_nonres =
      plotW("Wtrue_nonres", df.Filter("evtMode!=2 && evtMode != 3"));
  auto hw_resdis =
      plotW("Wtrue_resdis", df.Filter("evtMode==2 || evtMode == 3"));
  auto hw_mpi = plotW("Wtrue_mpi", df.Filter("npi >1 "));
  auto hw_res =
      plotW("Wtrue_res",
            df.Filter(
                [](int evtMode, int flag_delta, double Wtrue) {
                  return evtMode == 2 && (flag_delta == 1 || Wtrue <= 1.210);
                },
                {"evtMode", "flag_delta", "Wtrue"}, "RES Delta Cut"));

  auto ratio1 = plotpnratio("ratio1", df);
  auto ratio2 = plotpnratio2("ratio2", df);
  auto ratio2_nonres =
      plotpnratio2("ratio2_nonres", df.Filter("evtMode!=2 && evtMode != 3"));
  auto diff = plotdiff("diff", df);
  auto theta = plotangle("thetaL", df);

  scale_plot(EventRate, hw.GetPtr());
  scale_plot(EventRate, hw_hydrogen.GetPtr());
  scale_plot(EventRate, hw_hydrogen_res.GetPtr());
  scale_plot(EventRate, hw_nonres.GetPtr());
  scale_plot(EventRate, hw_resdis.GetPtr());
  scale_plot(EventRate, hw_mpi.GetPtr());
  scale_plot(EventRate, hw_res.GetPtr());

  scale_plot(EventRate, ratio1.GetPtr());
  scale_plot(EventRate, ratio2.GetPtr());
  scale_plot(EventRate, ratio2_nonres.GetPtr());
  scale_plot(EventRate, diff.GetPtr());
  scale_plot(EventRate, theta.GetPtr());

  auto file = std::make_unique<TFile>(output_file.c_str(), "RECREATE");
  file->cd();
  hw->Write();
  hw_hydrogen->Write();
  hw_hydrogen_res->Write();
  hw_nonres->Write();
  hw_resdis->Write();
  hw_mpi->Write();
  hw_res->Write();
  ratio1->Write();
  ratio2->Write();
  ratio2_nonres->Write();
  diff->Write();
  theta->Write();
  file->Close();
  if (argc == 4) {
    std::string snapshot_out = argv[3];
    df.Snapshot("tree", snapshot_out, {"IApN", "dpt", "dpl"});
  }
  return 0;
}