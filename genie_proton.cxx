#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <event.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <tools.h>

// constexpr size_t Nbins = 300;
int main(int argc, const char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <intput.json>" << std::endl;
    return 1;
  }
  nlohmann::json config;
  {
    std::ifstream ifs(argv[1]);
    if (!ifs.is_open()) {
      std::cerr << "Could not open file " << argv[1] << std::endl;
      return 1;
    }
    ifs >> config;
  }
  std::filesystem::path output_path = config["output_path"];
  if (!output_path.parent_path().empty()) {
    std::cout << "Creating output directory " << output_path.parent_path()
              << std::endl;
    std::filesystem::create_directories(output_path.parent_path());
  }
  auto file = std::make_unique<TFile>(output_path.c_str(), "RECREATE");
  ROOT::EnableImplicitMT();
  std::vector<std::string> names{};
  for (auto name : config["input_file_list"]) {
    // std::cout << "Adding file " << name << std::endl;
    names.push_back(name);
  }
  ROOT::RDataFrame d("gRooTracker", names);
  // d.Define("event", []() { return event{}; });
  auto dataset =
      d.Filter(
           [](ROOT::RVec<int> &StdHepPdg, TObjString &EvtCode) {
             return EvtCode.GetString().Contains("Weak[CC]") &&
                    StdHepPdg[0] == 14;
           },
           {"StdHepPdg", "EvtCode"}) // nu_mu CC filter
          .Define(
              "event",
              [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                 ROOT::RVec<int> &StdHepStatus, ROOT::RVec<double> &StdHepP4_,
                 TObjString &EvtCode) {
                double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                event e{};
                e.set_mode(get_mode_genie(EvtCode));
                for (int i = 0; i < StdHepN; ++i) {
                  auto pdg = StdHepPdg[i];
                  if (StdHepPdg[i] == 1000000010) {
                    pdg = 2112;
                  }
                  switch (StdHepStatus[i]) {
                  case 0:
                    e.add_particle_in(
                        pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                  case 1:
                    e.add_particle_out(
                        pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                  default:
                    break;
                  }
                }
                return e;
              },
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"})
          .Define("enu", [](event &e) { return e.get_enu(); }, {"event"})
          .Define("kaonE",
                  [](event &e) {
                    ROOT::RVec<double> ret;
                    if (e.count_out(321) > 0) {
                      auto p = e.get_leading_out(321);
                      ret.push_back(p.E() - p.M());
                    }
                    return ret;
                  },
                  {"event"})
          .Define("ProtonP",
                  [](event &e) {
                    ROOT::RVec<double> ret;
                    if (e.count_out(2212) > 0) {
                      auto p = e.get_leading_out(2212);
                      ret.push_back(p.P());
                    }
                    return ret;
                  },
                  {"event"});

  std::vector<ROOT::RDF::RResultPtr<TH1>> objs_list{};
  auto spline_file =
      get_object<TGraph>(config["spline_file"].get<std::string>(),
                         config["spline_path"].get<std::string>());
  auto enu = dataset.Histo1D("enu");
  auto add_plot = [&](auto &&dataset_, const char *name, const char *varname) {
    auto [tot, xsec] = get_xsec(enu.GetPtr(), spline_file.get());
    xsec *= 1 / 12. * 1e-38;
    objs_list
        .emplace_back(dataset_.Histo1D({name, varname, 200, 0., 0.}, varname))
        ->Scale(xsec / tot, "WIDTH");
  };
  add_plot(dataset, "kaonE", "kaonE");
  auto lowke = dataset.Filter(
      [](ROOT::RVec<double> &kaonE) {
        if (kaonE.size() == 0) {
          return false;
        }
        return kaonE[0] < 0.3;
      },
      {"kaonE"});
  objs_list.emplace_back(
      lowke.Histo1D({"enu_lowke", "enu_lowke", 200, 0., 0.}, "enu"));
  add_plot(lowke, "ProtonP", "ProtonP");
  save(objs_list, file);
  return 0;
}