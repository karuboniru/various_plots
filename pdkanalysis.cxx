#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TPie.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <analysis.h>
#include <chain_helper.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
// #include <tools.h>
// #define MAX_COUNT 128
// using std::string_literals::operator""s;

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
}

int main(int argc, char *argv[]) {
  TH1::AddDirectory(kFALSE);
  ROOT::EnableThreadSafety();
  gStyle->SetOptStat(0);
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  // std::array<int, 10> col{kRed, kGreen, kBlue, kMagenta, kCyan, kOrange,
  // kViolet, kGray, kYellow, kBlack};
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
  // analysis::plots instance;
  std::initializer_list<analysis::plots::cut_group_t> cut_group{
      {{"0pi",
        [](event &e) {
          return e.count_out(211) == 0 && e.count_out(-211) == 0 &&
                 e.count_out(111) == 0;
        }},
       {"1pip",
        [](event &e) {
          return e.count_out(211) == 1 && e.count_out(111) == 0 &&
                 e.count_out(-211) == 0;
        }},
       {"1pim",
        [](event &e) {
          return e.count_out(-211) == 1 && e.count_out(111) == 0 &&
                 e.count_out(211) == 0;
        }},
       {"1pi0",
        [](event &e) {
          return e.count_out(211) == 0 && e.count_out(111) == 1 &&
                 e.count_out(-211) == 0;
        }},
       {"1pi",
        [](event &e) {
          return e.count_out(211) + e.count_out(-211) + e.count_out(111) == 1;
        }},
       {"mpi",
        [](event &e) {
          return e.count_out(211) + e.count_out(-211) + e.count_out(111) >= 2;
        }},
       {"epi",
        [](event &e) {
          size_t count_e{}, count_pi{};
          for (const auto &[id, p] : e.get_particle_out()) {
            if (id == -11) {
              count_e++;
            } else if (id == 111) {
              count_pi++;
            } else {
              return false;
            }
          }
          return (count_e == 1) && (count_pi == 1);
        }}},
      {
          {"nofsi",
           [](event &e) {
             return (e.get_particle_nofsi().find(111)->second -
                     e.get_leading_out(111))
                        .P() < 2e-3;
           }},
          {"fsi",
           [](event &e) {
             return (e.get_particle_nofsi().find(111)->second -
                     e.get_leading_out(111))
                        .P() >= 2e-3;
           }},
          {"nointeraction",
           [](event &e) { return e.is_no_pion_interaction(); }},
          {"interaction", [](event &e) { return !e.is_no_pion_interaction(); }},
      },
      // {{"1ce",
      //   [](event &e) {
      //     auto &&nod = e.get_nod();
      //     return nod[5] == 1 && (!nod[6]) && (!nod[7]) && (!nod[8]);
      //   }}}
      // {{"0nucleon",
      //   [](event &e) { return e.count_out(2112) + e.count_out(2212) == 0; }},
      //  {"1nucleon",
      //   [](event &e) { return e.count_out(2112) + e.count_out(2212) == 1; }},
      //  {"2nucleon",
      //   [](event &e) { return e.count_out(2112) + e.count_out(2212) == 2; }},
      //  {"morenucleon",
      //   [](event &e) { return e.count_out(2112) + e.count_out(2212) > 2; }}},
      {{"nosorb", [](event &e) { return e.get_nod()[8] == 0; }}} //
  };
  std::initializer_list<analysis::plots::cut_t> cuts{};
  std::unordered_map<std::string, size_t> counts{}, counts_nonucleon{};
  std::mutex count_mutex;
  std::initializer_list<analysis::plots::analysis_funcs_t> analysis_funcs{
      {"leadingPi0",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(111)) {
           //  h->Fill(e.get_leading_out(111).P(), e.get_weight());
           vars.push_back(e.get_leading_out(111).P());
         }
         return vars;
       },
       {0., 0.8, 200}},
      {"before_fsi_P_pi0",
       [](event &e) {
         std::vector<double> vars;
         // h->Fill(e.get_particle_nofsi().find(111)->second.P(),
         // e.get_weight());
         vars.push_back(e.get_particle_nofsi().find(111)->second.P());
         return vars;
       },
       {0., 1.2, 100}},
      {"fsi_deltaP",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(111)) {
           //  h->Fill((e.get_particle_nofsi().find(111)->second -
           //  e.get_leading_out(111)).P(), e.get_weight());
           vars.push_back((e.get_particle_nofsi().find(111)->second -
                           e.get_leading_out(111))
                              .P());
         }
         return vars;
       },
       {0, 2.0e-1, 1000}},
      {"delta_inv_m",
       [](event &e) {
         std::vector<double> vars{};
         if (e.count_particle_out(111)) {
           auto po = e.get_leading_out(-11);
           auto pi = e.get_leading_out(111);
           auto p = e.get_particle_in().begin()->second;
           vars.push_back(p.M() - (po + pi).M());
         }
         return vars;
       },
       {-2.0e-2, 2.0e-2, 1000}},
      {"leadingPip",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(211)) {
           //  h->Fill(e.get_leading_out(211).P(), e.get_weight());
           vars.push_back(e.get_leading_out(211).P());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"leadingPim",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(-211)) {
           //  h->Fill(e.get_leading_out(-211).P(), e.get_weight());
           vars.push_back(e.get_leading_out(-211).P());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"pe",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(-11)) {
           vars.push_back(e.get_particle_out().find(-11)->second.P());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"energy_positron",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(-11)) {
           auto &&particle = e.get_particle_out().find(-11)->second;
           vars.push_back(particle.E());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"energy_pi0",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(111)) {
           auto &&particle = e.get_leading_out(111);
           vars.push_back(particle.E());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"angle_e",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(-11) && e.count_particle_out(111)) {
           //  h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_out().find(111)->second.Vect()),
           //  e.get_weight());
           vars.push_back(e.get_particle_out().find(-11)->second.Angle(
                              e.get_particle_out().find(111)->second.Vect()) /
                          M_PI * 180.);
         }
         return vars;
       },
       {0., 180., 100}},
      {"angle_e_beforefsi",
       [](event &e) {
         std::vector<double> vars;
         //  h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_nofsi().find(111)->second.Vect()),
         //  e.get_weight());
         vars.push_back(e.get_particle_out().find(-11)->second.Angle(
                            e.get_particle_nofsi().find(111)->second.Vect()) /
                        M_PI * 180.);
         return vars;
       },
       {0., 180., 100}},
      {"ep_inv_mass",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(-11) == 1 && e.count_particle_out(111)) {
           //  h->Fill((e.get_particle_out().find(-11)->second +
           //  e.get_leading_out(111)).M(), e.get_weight());
           vars.push_back(
               (e.get_particle_out().find(-11)->second + e.get_leading_out(111))
                   .M());
         }
         return vars;
       },
       {0., 1.2, 100}},
      {"leading_nucleon_momentum",
       [](event &e) {
         double max{-1.};
         if (e.count_out(2212)) {
           max = std::max(max, e.get_leading_out(2212).P());
         }
         if (e.count_out(2112)) {
           max = std::max(max, e.get_leading_out(2112).P());
         }

         return std::vector<double>{max};
       },
       {0., 1.2, 100}},
      {"pi0_nucleon_inv_mass",
       [](event &e) {
         std::vector<double> vars;
         TLorentzVector nucleon{};
         double maxp{};
         auto pi0 = e.get_leading_out(111);
         bool is_nucleon{false};
         for (auto &&[id, p4] : e.get_particle_out()) {
           if (std::abs(id) == 2212 || std::abs(id) == 2112) {
             is_nucleon = true;
             if (p4.P() > maxp) {
               maxp = p4.P();
               nucleon = p4;
             }
           }
         }
         if (is_nucleon) {
           vars.push_back((nucleon + pi0).M());
         }
         return vars;
       },
       {0.93, 1.5, 100}},
      {"allpi0_nucleon_inv_mass",
       [](event &e) {
         std::vector<double> vars;
         TLorentzVector nucleon{}, pi0{};
         double maxp{};
         //  auto pi0 = e.get_leading_out(111);
         bool is_nucleon{false};
         for (auto &&[id, p4] : e.get_particle_out()) {
           if (std::abs(id) == 2212 || std::abs(id) == 2112) {
             is_nucleon = true;
             if (p4.P() > maxp) {
               maxp = p4.P();
               nucleon = p4;
             }
           }
           if (id == 111 || id == 211 || id == -211) {
             pi0 += p4;
           }
         }
         if (is_nucleon) {
           vars.push_back((nucleon + pi0).M());
         }
         return vars;
       },
       {0.93, 1.5, 100}},
      {"sum_of_nucleon_momentum",
       [](event &e) {
         double sum{0};
         for (auto &&p : e.get_particle_out(2212)) {
           sum += p.second.P();
         }
         for (auto &&p : e.get_particle_out(2112)) {
           sum += p.second.P();
         }

         return std::vector<double>{sum};
       },
       {0., 1.2, 100}},
      {"leading_proton_momentum",
       [](event &e) {
         if (e.count_out(2212)) {
           return std::vector<double>{e.get_leading_out(2212).P()};
         }
         return std::vector<double>{};
       },
       {0., 1.2, 100}},
      {"leading_neutron_momentum",
       [](event &e) {
         if (e.count_out(2112)) {
           return std::vector<double>{e.get_leading_out(2112).P()};
         }
         return std::vector<double>{};
       },
       {0., 1.2, 100}},
      {"deltaE",
       [](event &e) {
         auto &&p4 = e.get_particle_in().begin()->second;
         double p = p4.P();

         return std::vector<double>{sqrt(p * p + 0.938 * 0.938) - p4.E()};
       },
       {0., 0.1, 100}},
      {"init_proton_momentum",
       [](event &e) {
         auto &&p4 = e.get_particle_in().begin()->second;
         double p = p4.P();

         return std::vector<double>{p};
       },
       {0., 1.2, 500}},
      {"init_proton_energy",
       [](event &e) {
         auto &&p4 = e.get_particle_in().begin()->second;
         double energy = p4.E();

         return std::vector<double>{energy};
       },
       {0.6, 1.1, 500}},
      {"init_proton_invm",
       [](event &e) {
         auto &&p4 = e.get_particle_in().begin()->second;
         double invm = p4.M();

         return std::vector<double>{invm};
       },
       {0.6, 1.2, 500}},
      {"init_proton_kineticenergy",
       [](event &e) {
         auto &&p4 = e.get_particle_in().begin()->second;
         double energy = p4.E() - p4.M();

         return std::vector<double>{energy};
       },
       {0., 0.5, 500}},
      {"sumepi",
       [](event &e) {
         std::vector<double> vars{};
         //  auto &&p4 = e.get_particle_in().begin()->second;
         if (e.count_particle_out(-11) == 1 && e.count_particle_out(111)) {
           //  h->Fill((e.get_particle_out().find(-11)->second +
           //  e.get_leading_out(111)).M(), e.get_weight());
           vars.push_back(e.get_particle_out().find(-11)->second.E() +
                          e.get_leading_out(111).E());
         }

         return vars;
       },
       {0., 1.2, 100}},
      {"Q2",
       [](event &e) {
         if (e.count_particle_out(111)) {
           auto &&p_init = e.get_particle_nofsi(111).begin()->second;
           auto &&p_final = e.get_leading_out(111);
           return std::vector<double>{-(p_init - p_final).M2()};
         }
         return std::vector<double>{};
       },
       {0., 1.2, 100}},
      {"",
       [&](event &e) -> std::vector<double> {
         auto channel = e.get_channelname();
         auto channel_nonucleon = e.get_channelname_no_nucleon();
         {
           std::lock_guard<std::mutex> lk(count_mutex);
           counts[channel]++;
           counts_nonucleon[channel_nonucleon]++;
         }
         return {};
       },
       {0., 1.2, 100}}};
  std::initializer_list<analysis::plots::axis_pair_t> histsnames{
      {"ep_inv_mass", "angle_e_beforefsi"},
      {"ep_inv_mass", "leadingPi0"},
      {"leadingPi0", "leading_nucleon_momentum"},
      {"leadingPi0", "leading_proton_momentum"},
      {"leadingPi0", "leading_neutron_momentum"},
      {"leadingPi0", "sum_of_nucleon_momentum"},
      {"leadingPi0", "pi0_nucleon_inv_mass"},
      {"pi0_nucleon_inv_mass", "leadingPi0"},
      {"leadingPi0", "allpi0_nucleon_inv_mass"},
      {"leadingPi0", "angle_e"},
      {"energy_positron", "energy_pi0"},
      {"init_proton_momentum", "deltaE"}};
  size_t total{};
  if (config["type"] == "genie") {
    auto man = run_manager_genie::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    total = man.event_count;
    man.plot.save(output_path);
  } else if (config["type"] == "nuwro") {
    auto man = run_manager_nuwro::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    total = man.event_count;
    man.plot.save(output_path);
  } else if (config["type"] == "gibuu") {
    auto man = run_manager_gibuu::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    // total = man.event_count;
    man.plot.save(output_path);
  } else if (config["type"] == "gibuupdk") {
    auto man = run_manager_gibuu_pdk::run_analysis(config, analysis_funcs,
                                                   cut_group, cuts, histsnames);
    total = man.count;
    man.plot.save(output_path);
  } else {
    std::cerr << "Unknown type: " << config["type"] << std::endl;
    return 1;
  }
  // std::cout << "Writing to " << output_path << std::endl;
  // instance.save(output_path);
  {
    std::vector<std::pair<std::string, double>> counts_vec;
    for (const auto &[channel, count] : counts) {
      // std::cout << channel << ": " << count << std::endl;
      counts_vec.push_back({channel, (double)count / total});
    }
    std::sort(counts_vec.begin(), counts_vec.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    double others{};
    std::vector<std::pair<std::string, double>> counts_vec_pie;
    for (const auto &[channel, count] : counts_vec) {
      std::cout << channel << ": \t\t\t" << count << std::endl;
      if (count > 0.02) {
        counts_vec_pie.push_back({channel, count});
      } else {
        others += count;
      }
    }
    counts_vec_pie.push_back({"others", others});
    make_pie_plot(counts_vec_pie, output_path.string() + "counts.pdf");
  }
  {
    std::vector<std::pair<std::string, double>> counts_vec;
    for (const auto &[channel, count] : counts_nonucleon) {
      // std::cout << channel << ": " << count << std::endl;
      counts_vec.push_back({channel, (double)count / total});
    }
    std::sort(counts_vec.begin(), counts_vec.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    double others{};
    std::vector<std::pair<std::string, double>> counts_vec_pie;
    size_t c{};
    for (const auto &[channel, count] : counts_vec) {
      std::cout << channel << ": \t\t\t" << count << std::endl;
      if (c < 4) {
        counts_vec_pie.push_back({channel, count});
      } else {
        others += count;
      }
      c++;
    }
    counts_vec_pie.push_back({"others", others});
    make_pie_plot(counts_vec_pie,
                  output_path.string() + "counts_nonucleon.pdf");
  }
  return 0;
}