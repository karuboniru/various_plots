#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <event.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <tools.h>

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
          .Define("rho_P",
                  [](event &e) {
                    return e.get_particle_nofsi(113).begin()->second;
                  },
                  {"event"})
          .Define("prho", [](TLorentzVector &p) { return p.P(); }, {"rho_P"})
          .Define("ppion", [](TLorentzVector &p) { return p.P(); },
                  {"sum_p_pion_system"})
          .Define("mrho", [](TLorentzVector &p) { return p.M(); }, {"rho_P"})
          .Define("W", [](event &e) { return e.getW(); }, {"event"})
          .Define("channelname",
                  [](event &e) { return e.get_channelname_no_nucleon(); },
                  {"event"});

  auto dataset_2pi =
      dataset
          .Filter([](int pion_count) { return pion_count == 2; },
                  {"pion_count"})
          .Define("deltaP",
                  [](TLorentzVector &sum_p, TLorentzVector &rho_P) {
                    return (sum_p - rho_P).P();
                  },
                  {"sum_p_pion_system", "rho_P"})
          .Define("M2pi", [](TLorentzVector &sum_p) { return sum_p.M(); },
                  {"sum_p_pion_system"})
          .Define("deltaM",
                  [](double M2pi, TLorentzVector &rho_P) {
                    return M2pi - rho_P.M();
                  },
                  {"M2pi", "rho_P"})
          .Define("angle_2pi",
                  [](event &e) {
                    std::array<const TLorentzVector *, 2> pions{};
                    int i = 0;
                    for (const auto &[id, p4] : e.get_particle_out()) {
                      if (id == 211 || id == 111 || id == -211) {
                        pions[i++] = &p4;
                      }
                      if (i == 2) {
                        break;
                      }
                    }
                    return pions[0]->Angle(pions[1]->Vect());
                  },
                  {"event"});
  auto dataset_2pi_bc = dataset_2pi.Filter("pipcount - pimcount == 0");
  auto dataset_2pi_bc_nofsi = dataset_2pi_bc.Filter("deltaP < 1e-4");
  // hist->SetDirectory(file.get());
  // hist->SaveAs("test.root");
  std::vector<ROOT::RDF::RResultPtr<TH1D>> objs_list{};
  auto count_2pi = dataset_2pi.Count().GetValue();
  auto count_all = dataset.Count().GetValue();
  objs_list
      .emplace_back(dataset_2pi.Histo1D(
          {"W_2pi", "W(2#pi) ; #it{W} (GeV); Count (1/GeV)", 128u, 0.2, 0.9},
          "W_pion_system"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset_2pi_bc_nofsi.Histo1D(
          {"W_nofsi", "W(2#pi) ; #it{W} (GeV); Count (1/GeV)", 128u, 0.2, 0.9},
          "W_pion_system"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset_2pi_bc.Histo1D(
          {"W_2pi_noce",
           "W(#pi^{+}#pi^{-}/ 2#pi^{0}) ; #it{W} (GeV); Count (1/GeV)", 128u,
           0.2, 0.9},
          "W_pion_system"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(
          dataset.Histo1D({"W_all", "", 128u, 0.2, 0.9}, "W_pion_system"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset_2pi.Histo1D(
          {"deltaP", "; #Delta #it{p} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u,
           0, 0},
          "deltaP"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(
          dataset_2pi.Filter("deltaP < 1e-4")
              .Histo1D({"deltaP_filt",
                        "; #Delta #it{p} (GeV/#it{c}); Count (1/GeV/#it{c})",
                        128u, 0, 0},
                       "deltaP"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list.emplace_back(dataset_2pi.Histo1D("deltaM"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset.Histo1D(
          {"prho", "; #it{p_{#rho}} (GeV/#it{c}); Count (1/GeV/#it{c})", 128u,
           0, 0},
          "prho"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list.emplace_back(dataset.Histo1D("W"))->Scale(1. / count_all, "WIDTH");
  objs_list.emplace_back(dataset_2pi.Histo1D("ppion"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list.emplace_back(dataset_2pi.Histo1D("mrho"))
      ->Scale(1. / count_all, "WIDTH");

  objs_list.emplace_back(dataset_2pi.Histo1D("angle_2pi"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset_2pi_bc.Histo1D(
          {"angle_2pi_noce", "angle_2pi", 128, 0., 0.}, "angle_2pi"))
      ->Scale(1. / count_all, "WIDTH");
  objs_list
      .emplace_back(dataset_2pi_bc_nofsi.Histo1D(
          {"angle_2pi_nofsi", "angle_2pi", 128, 0., 0.}, "angle_2pi"))
      ->Scale(1. / count_all, "WIDTH");

  // objs_list.emplace_back(dataset.Histo2D("ppion", "prho"));
  save(objs_list, file);

  std::cout << "2pi events: " << count_2pi << std::endl;
  std::cout << "all events: " << count_all << std::endl;
  std::cout << "2pi/all: " << (double)count_2pi / count_all << std::endl;

  {
    auto array_channels = dataset.Take<std::string>("channelname");
    std::unordered_map<std::string, int> counts{};
    std::vector<std::pair<std::string, double>> counts_vec;
    for (auto &channel : array_channels) {
      counts[channel]++;
    }
    for (auto &channel : counts) {
      counts_vec.emplace_back(channel.first,
                              channel.second / ((double)count_all));
    }
    std::sort(counts_vec.begin(), counts_vec.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    for (auto &channel : counts_vec) {
      std::cout << channel.first << ": " << channel.second << "\n";
    }
  }
  return 0;
  // file->Write();
}