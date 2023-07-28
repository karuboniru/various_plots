#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <plot.h>
#include <string>
#include <tools.h>
#include <unordered_map>

int main(int argc, const char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <intput.json>" << "\n";
    return 1;
  }
  nlohmann::json config;
  {
    std::ifstream ifs(argv[1]);
    if (!ifs.is_open()) {
      std::cerr << "Could not open file " << argv[1] << "\n";
      return 1;
    }
    ifs >> config;
  }
  std::filesystem::path output_path = config["output_path"];
  if (!output_path.parent_path().empty()) {
    std::cout << "Creating output directory " << output_path.parent_path()
              << "\n";
    std::filesystem::create_directories(output_path.parent_path());
  }
  auto file = std::make_unique<TFile>(output_path.c_str(), "RECREATE");
  ROOT::EnableImplicitMT();
  std::vector<std::string> names{};
  for (auto name : config["input_file_list"]) {
    // std::cout << "Adding file " << name << "\n";
    if (std::filesystem::exists(name))
      names.push_back(name);
  }
  ROOT::RDataFrame d("outtree", names);
  // d.Define("event", []() { return event{}; });
  auto dataset =
      d.Define(
           "event",
           [](int nparticles, ROOT::RVec<double> &P_u, ROOT::RVec<int> &status,
              ROOT::RVec<int> &pdg) {
             double(*P)[4] = (double(*)[4]) & P_u[0];
             event e{};
             for (int i = 0; i < nparticles; ++i) {
               switch (status[i]) {
               case 0:
                 e.add_particle_in(pdg[i], TLorentzVector(P[i][1], P[i][2],
                                                          P[i][3], P[i][0]));
                 break;
               case 1:
                 e.add_particle_out(pdg[i], TLorentzVector(P[i][1], P[i][2],
                                                           P[i][3], P[i][0]));
                 break;
               case 2:
                 e.add_particle_nofsi(pdg[i], TLorentzVector(P[i][1], P[i][2],
                                                             P[i][3], P[i][0]));
                 break;
               default:
                 break;
               }
             }
             return e;
           },
           {"nparticles", "P", "status", "pdg"})
          .Define("pion_count",
                  [](event &e) -> int {
                    return e.count_out(211) + e.count_out(111) +
                           e.count_out(-211);
                  },
                  {"event"})
          .Define("pipcount", [](event &e) { return e.count_out(211); },
                  {"event"})
          .Define("pi0count", [](event &e) { return e.count_out(111); },
                  {"event"})
          .Define("pimcount", [](event &e) { return e.count_out(-211); },
                  {"event"})
          .Define("sum_p_pion_system",
                  [](event &e) {
                    TLorentzVector sum{};
                    for (const auto &particle : e.get_particle_out()) {
                      if ((std::abs(particle.first) == 211) ||
                          (particle.first == 111)) {
                        sum += particle.second;
                      }
                    }
                    return sum;
                  },
                  {"event"})
          .Define("W_pion_system", [](TLorentzVector &p) { return p.M(); },
                  {"sum_p_pion_system"})
          .Define("eta_P",
                  [](event &e) {
                    return e.get_particle_nofsi(221).begin()->second;
                  },
                  {"event"})
          .Define("peta", [](TLorentzVector &p) { return p.P(); }, {"eta_P"})
          .Define("ppion", [](TLorentzVector &p) { return p.P(); },
                  {"sum_p_pion_system"})
          .Define("meta", [](TLorentzVector &p) { return p.M(); }, {"eta_P"})
          .Define("channelname",
                  [](event &e) { return e.get_channelname_no_nucleon(); },
                  {"event"});

  auto dataset_3pi =
      dataset
          .Filter([](int pion_count) { return pion_count == 3; },
                  {"pion_count"})
          .Define("pip_count", [](event &e) { return e.count_out(211); },
                  {"event"})
          .Define("deltaP",
                  [](TLorentzVector &sum_p, TLorentzVector &eta_P) {
                    return (sum_p - eta_P).P();
                  },
                  {"sum_p_pion_system", "eta_P"})
          .Define("deltaM",
                  [](TLorentzVector &sum_p, TLorentzVector &eta_P) {
                    return sum_p.M() - eta_P.M();
                  },
                  {"sum_p_pion_system", "eta_P"});

  // hist->SetDirectory(file.get());
  // hist->SaveAs("test.root");
  std::vector<ROOT::RDF::RResultPtr<TH1D>> objs_list{};
  auto count = dataset.Count().GetValue();
  objs_list
      .emplace_back(dataset_3pi.Histo1D(
          {"W_3pi",
           "#it{W} of #pi system, 3#pi cut; #it{W} (GeV); Count (1/GeV)", 128u,
           0, 0},
          "W_pion_system"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(dataset.Filter("pion_count == 2")
                        .Histo1D({"W_2pi",
                                  "#it{W} of #pi system, 2#pi cut; #it{W} "
                                  "(GeV); Count (1/GeV)",
                                  128u, 0, 0},
                                 "W_pion_system"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(dataset.Histo1D(
          {"W_all", "#it{W} of #pi system; #it{W} (GeV); Count (1/GeV)", 128u,
           0, 0},
          "W_pion_system"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(dataset_3pi.Histo1D(
          {"deltaP", "; #Delta #it{p} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u,
           0, 0},
          "deltaP"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(dataset_3pi.Histo1D(
          {"deltaM",
           "; #Delta #it{m} (GeV/#it{c^{2}}); Count (1/GeV/#it{c^{2}})", 128u,
           0, 0},
          "deltaM"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(dataset.Histo1D(
          {"peta", "; #it{p_{#eta}} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u,
           0, 0},
          "peta"))
      ->Scale(1. / count, "WIDTH");
  objs_list.emplace_back(dataset_3pi.Histo1D("ppion"))
      ->Scale(1. / count, "WIDTH");
  objs_list.emplace_back(dataset_3pi.Histo1D("meta"))
      ->Scale(1. / count, "WIDTH");

  auto event_pi_nofsi = dataset_3pi.Filter("deltaM < 3e-5");
  objs_list
      .emplace_back(event_pi_nofsi.Histo1D(
          {"deltaP_select",
           "; #Delta #it{p} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u, 0, 0},
          "deltaP"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(
          event_pi_nofsi.Filter("deltaP < 1e-4")
              .Histo1D({"deltaP_select_1",
                        "; #Delta #it{p} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "deltaP"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(event_pi_nofsi.Histo1D(
          {"peta_pi_nofsi",
           "; #it{p}_{#eta} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u, 0, 0},
          "peta"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(
          event_pi_nofsi.Filter("deltaP < 3e-5")
              .Histo1D({"peta_pi_nofsi_eta_nofsi",
                        "; #it{p}_{#eta} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "peta"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(
          dataset.Filter("pion_count == 0")
              .Histo1D({"peta_all_absorb",
                        "; #it{p}_{#eta} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "peta"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(
          dataset_3pi.Filter("pip_count == 0")
              .Histo1D({"peta_3pi0",
                        "; #it{p}_{#eta} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "peta"))
      ->Scale(1. / count, "WIDTH");
  objs_list
      .emplace_back(
          dataset_3pi.Filter("pip_count == 1")
              .Histo1D({"peta_pipm0",
                        "; #it{p}_{#eta} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "peta"))
      ->Scale(1. / count, "WIDTH");
  save(objs_list, file);

  {
    auto array_channels = dataset.Take<std::string>("channelname");
    std::unordered_map<std::string, int> counts{};
    std::vector<std::pair<std::string, double>> counts_vec;
    for (auto &channel : array_channels) {
      counts[channel]++;
    }
    for (auto &channel : counts) {
      counts_vec.emplace_back(channel.first, channel.second / ((double)count));
    }
    std::sort(counts_vec.begin(), counts_vec.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    for (auto &channel : counts_vec) {
      std::cout << channel.first << ": " << channel.second << "\n";
    }
  }

  auto count_3pi = dataset_3pi.Count().GetValue();
  auto count_all = dataset.Count().GetValue();
  std::cout << "3pi events: " << count_3pi << "\n";
  std::cout << "all events: " << count_all << "\n";
  std::cout << "3pi/all: " << (double)count_3pi / count_all << "\n";
  return 0;
  // file->Write();
}