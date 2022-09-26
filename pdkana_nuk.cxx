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
      {{"fsi",
        [](event &e) {
          return (e.get_particle_nofsi().find(321)->second -
                  e.get_leading_out(321))
                     .P() > 2e-4;
        }},
       {"nofsi",
        [](event &e) {
          return (e.get_particle_nofsi().find(321)->second -
                  e.get_leading_out(321))
                     .P() < 2e-4;
        }}},
      {{"nonucleon",
        [](event &e) { return e.count_out(2212) + e.count_out(2112) == 0; }},
       {"nucleon",
        [](event &e) { return e.count_out(2212) + e.count_out(2112) > 0; }}},
      {{"k0", [](event &e) { return e.count_out(311) != 0; }},
       {"Kp", [](event &e) { return e.count_out(321) != 0; }},
       {"nok",
        [](event &e) { return (e.count_out(321) + e.count_out(311)) == 0; }}}};
  std::initializer_list<analysis::plots::cut_t> cuts{};
  std::unordered_map<std::string, size_t> counts{}, counts_nonucleon{};
  std::mutex count_mutex;
  std::initializer_list<analysis::plots::analysis_funcs_t> analysis_funcs{
      {"kaonp_momentum",
       [](event &e) {
         std::vector<double> vars;
         for (auto &p : e.get_particle_out(321)) {
           vars.push_back(p.second.P());
         }
         return vars;
       },
       {0, .7, 50}},
      {"kaon0_momentum",
       [](event &e) {
         std::vector<double> vars;
         for (auto &p : e.get_particle_out(311)) {
           vars.push_back(p.second.P());
         }
         return vars;
       },
       {0, .7, 50}},
      {"kaon_beforefsi_momentum",
       [](event &e) {
         std::vector<double> vars;
         for (auto &p : e.get_particle_nofsi(321)) {
           vars.push_back(p.second.P());
         }
         return vars;
       },
       {0, .7, 50}},
      {"FSI_delta_P",
       [](event &e) {
         return std::vector<double>{
             (e.get_particle_nofsi().find(321)->second - e.get_leading_out(321))
                 .P()};
       },
       {0, 1., 500e2}},
      {"WnuK",
       [](event &e) {
         std::vector<double> vars{};
         if (e.count_out(321) == 0 || e.count_out(14) == 0)
           return vars;
         auto &&nu = e.get_particle_nofsi().find(14)->second;
         auto &&k = e.get_particle_nofsi().find(321)->second;
         auto W = nu + k;
         vars.push_back(W.M());
         return vars;
       },
       {0.85, 0.95, 100}},
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
      {"angle",
       [](event &e) {
         std::vector<double> vars;
         if (e.count_particle_out(14) && e.count_particle_out(321)) {
           //  h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_out().find(111)->second.Vect()),
           //  e.get_weight());
           vars.push_back(e.get_particle_out().find(14)->second.Angle(
                              e.get_particle_out().find(321)->second.Vect()) /
                          M_PI * 180.);
         }
         return vars;
       },
       {0., 180., 100}},
      {"angle_beforefsi",
       [](event &e) {
         std::vector<double> vars;
         //  h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_nofsi().find(111)->second.Vect()),
         //  e.get_weight());
         vars.push_back(e.get_particle_out().find(14)->second.Angle(
                            e.get_particle_nofsi().find(321)->second.Vect()) /
                        M_PI * 180.);
         return vars;
       },
       {0., 180., 100}},
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
      {"Wknucl",
       [](event &e) -> std::vector<double> {
         std::vector<double> vars;
         int kid = e.count_out(321) ? 321 : 311;
         auto p = e.get_leading_out(2212);
         auto n = e.get_leading_out(2112);
         auto k = e.get_leading_out(kid);
         auto &leading_nucleon = p.P() > n.P() ? p : n;
         auto W = leading_nucleon + k;
         vars.push_back(W.M());
         return vars;
       },
       {1.35, 1.6, 300}},
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
      {"WnuK", "kaonp_momentum"},
      {"WnuK", "angle"},
      {"kaon0_momentum", "leading_nucleon_momentum"},
      {"kaonp_momentum", "leading_nucleon_momentum"},
  };
  size_t total{};
  if (config["type"] == "genie") {
    auto man = run_manager_genie::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    man.plot.save(output_path);
    total = man.event_count;
  } else if (config["type"] == "nuwro") {
    auto man = run_manager_nuwro::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    man.plot.save(output_path);
    total = man.event_count;
  } else if (config["type"] == "gibuu") {
    auto man = run_manager_gibuu::run_analysis(config, analysis_funcs,
                                               cut_group, cuts, histsnames);
    // total = man.event_count;
    man.plot.save(output_path);
  } else if (config["type"] == "gibuupdk") {
    auto man = run_manager_gibuu_pdk::run_analysis(config, analysis_funcs,
                                                   cut_group, cuts, histsnames);
    man.plot.save(output_path);
    total = man.count;
  } else {
    std::cerr << "Unknown type: " << config["type"] << std::endl;
    return 1;
  }
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