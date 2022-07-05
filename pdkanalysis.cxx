#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
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
// #include <tools.h>
// #define MAX_COUNT 128
using std::string_literals::operator""s;

void make_pie_plot(std::vector<std::pair<std::string, double>> data, std::string filename) {
    constexpr std::array<int, 10> col{kRed, kBlue, kViolet, kYellow, kOrange, kGreen, kGray, kTeal, kPink};
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
    // std::array<int, 10> col{kRed, kGreen, kBlue, kMagenta, kCyan, kOrange, kViolet, kGray, kYellow, kBlack};
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
        std::cout << "Creating output directory " << output_path.parent_path() << std::endl;
        std::filesystem::create_directories(output_path.parent_path());
    }
    // analysis::plots instance;
    std::initializer_list<analysis::plots::cut_group_t> cut_group{
        {{"0pi", [](event &e) { return e.count_out(211) == 0 && e.count_out(-211) == 0 && e.count_out(111) == 0; }},
         {"1pip", [](event &e) { return e.count_out(211) == 1 && e.count_out(111) == 0 && e.count_out(-211) == 0; }},
         {"1pim", [](event &e) { return e.count_out(-211) == 1 && e.count_out(111) == 0 && e.count_out(211) == 0; }},
         {"1pi0", [](event &e) { return e.count_out(211) == 0 && e.count_out(111) == 1 && e.count_out(-211) == 0; }},
         {"1pi", [](event &e) { return e.count_out(211) + e.count_out(-211) + e.count_out(111) == 1; }},
         {"mpi", [](event &e) { return e.count_out(211) > 0 && e.count_out(111) > 0; }},
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
        {{"nofsi", [](event &e) { return (e.get_particle_nofsi().find(111)->second - e.get_leading_out(111)).P() < 2e-3; }},
         {"fsi", [](event &e) { return (e.get_particle_nofsi().find(111)->second - e.get_leading_out(111)).P() >= 2e-3; }}}};
    std::initializer_list<analysis::plots::cut_t> cuts{};
    std::unordered_map<std::string, size_t> counts{}, counts_nonucleon{};
    std::mutex count_mutex;
    std::initializer_list<analysis::plots::analysis_funcs_t> analysis_funcs{
        {"leadingPi0",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(111)) {
                 h->Fill(e.get_leading_out(111).P(), e.get_weight());
             }
         }},
        {"before_fsi_P_pi0", [](event &e, TH1 *h) { h->Fill(e.get_particle_nofsi().find(111)->second.P(), e.get_weight()); }},
        {"fsi_deltaP",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(111)) {
                 h->Fill((e.get_particle_nofsi().find(111)->second - e.get_leading_out(111)).P(), e.get_weight());
             }
         }},
        {"leadingPip",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(211)) {
                 h->Fill(e.get_leading_out(211).P(), e.get_weight());
             }
         }},
        {"leadingPim",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-211)) {
                 h->Fill(e.get_leading_out(-211).P(), e.get_weight());
             }
         }},
        {"pe",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11)) {
                 h->Fill(e.get_particle_out().find(-11)->second.P(), e.get_weight());
             }
         }},
        {"angle_e",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11) && e.count_particle_out(111)) {
                 h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_out().find(111)->second.Vect()), e.get_weight());
             }
         }},
        {"angle_e_beforefsi",
         [](event &e, TH1 *h) { h->Fill(e.get_particle_out().find(-11)->second.Angle(e.get_particle_nofsi().find(111)->second.Vect()), e.get_weight()); }},
        {"ep_inv_mass",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11) == 1 && e.count_particle_out(111)) {
                 h->Fill((e.get_particle_out().find(-11)->second + e.get_leading_out(111)).M(), e.get_weight());
             }
         }},
        {"2dep_inv_mass_ppi",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11) == 1 && e.count_particle_out(111) == 1) {
                 dynamic_cast<TH2D *>(h)->Fill((e.get_particle_out().find(-11)->second + e.get_particle_out().find(111)->second).M(),
                                               e.get_particle_out().find(111)->second.P(), e.get_weight());
             }
         }},
        {"2dpep",
         [](event &e, TH1 *h) {
             dynamic_cast<TH2D *>(h)->Fill(e.get_particle_in().find(1000060120)->second.E(), e.get_particle_in().find(1000060120)->second.P(), e.get_weight());
         }},
        {"2d_pp0_mass_v_leading_nucl_momentum",
         [](event &e, TH1 *h) {
             double pp0_mass{-1.};
             if (e.count_particle_out(-11) == 1 && e.count_particle_out(111)) {
                 pp0_mass = (e.get_particle_out().find(-11)->second + e.get_leading_out(111)).M();
             }
             double nucl_momentum{};
             nucl_momentum = std::max(nucl_momentum, e.get_leading_out(2212).P());
             nucl_momentum = std::max(nucl_momentum, e.get_leading_out(2112).P());
             if (pp0_mass > 0)
                 dynamic_cast<TH2D *>(h)->Fill(pp0_mass, nucl_momentum, e.get_weight());
         }},
        {"2d_pp0_v_leading_nucl_momentum",
         [](event &e, TH1 *h) {
             //  double pp0_mass{-1.};
             //  if (e.count_particle_out(-11) == 1 && e.count_particle_out(111)) {
             //      pp0_mass = (e.get_particle_out().find(-11)->second + e.get_leading_out(111)).M();
             //  }
             double nucl_momentum{};
             nucl_momentum = std::max(nucl_momentum, e.get_leading_out(2212).P());
             nucl_momentum = std::max(nucl_momentum, e.get_leading_out(2112).P());
             if (nucl_momentum > 0 && e.count_particle_out(111))
                 dynamic_cast<TH2D *>(h)->Fill(e.get_particle_out().find(111)->second.P(), nucl_momentum, e.get_weight());
         }},
        {"2d_pp0_mass_v_epi_angle",
         [](event &e, TH1 *h) {
             double epi_angle{-1.};
             if (e.count_particle_out(-11) && e.count_particle_out(111)) {
                 epi_angle = e.get_particle_out().find(-11)->second.Angle(e.get_particle_out().find(111)->second.Vect());
             }
             if (epi_angle > 0)
                 dynamic_cast<TH2D *>(h)->Fill(e.get_particle_out().find(-11)->second.P(), epi_angle, e.get_weight());
         }},
        {"2dmass_v_mom",
         [](event &e, TH1 *h) {
             auto &&particle = e.get_particle_in().find(1000060120)->second;
             dynamic_cast<TH2D *>(h)->Fill(particle.M(), particle.P(), e.get_weight());
         }},
        {"2dE_v_mom",
         [](event &e, TH1 *h) {
             auto &&particle = e.get_particle_in().find(1000060120)->second;
             dynamic_cast<TH2D *>(h)->Fill(particle.E(), particle.P(), e.get_weight());
         }},
        {"2ddeltaE_v_mom",
         [](event &e, TH1 *h) {
             auto &&particle = e.get_particle_in().find(1000060120)->second;
             dynamic_cast<TH2D *>(h)->Fill(-particle.E() + sqrt(particle.P() * particle.P() + 0.953 * 0.953), particle.P(), e.get_weight());
         }},
        {"2dep_epi",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11) && e.count_particle_out(111)) {
                 auto &&p = e.get_leading_out(-11);
                 auto &&pi = e.get_leading_out(111);
                 dynamic_cast<TH2D *>(h)->Fill(p.E() - p.M(), pi.E() - pi.M(), e.get_weight());
             }
         }},
        {"ep_o_epi",
         [](event &e, TH1 *h) {
             if (e.count_particle_out(-11) && e.count_particle_out(111)) {
                 auto &&p = e.get_leading_out(-11);
                 auto &&pi = e.get_leading_out(111);
                 h->Fill((p.E() - p.M()) / (pi.E() - pi.M()), e.get_weight());
             }
         }},
        {"init_p_inv_mass", [](event &e, TH1 *h) { h->Fill(e.get_particle_in().find(1000060120)->second.M(), e.get_weight()); }},
        {"", [&](event &e, TH1 *) {
             auto channel = e.get_channelname();
             auto channel_nonucleon = e.get_channelname_no_nucleon();
             {
                 std::lock_guard<std::mutex> lk(count_mutex);
                 counts[channel]++;
                 counts_nonucleon[channel_nonucleon]++;
             }
         }}};
    size_t total{};
    if (config["type"] == "genie") {
        auto man = run_manager_genie::run_analysis(config, analysis_funcs, cut_group, cuts);
        total = man.event_count;
        man.plot.save(output_path);
    } else if (config["type"] == "nuwro") {
        auto man = run_manager_nuwro::run_analysis(config, analysis_funcs, cut_group, cuts);
        total = man.event_count;
        man.plot.save(output_path);
    } else if (config["type"] == "gibuu") {
        auto man = run_manager_gibuu::run_analysis(config, analysis_funcs, cut_group, cuts);
        // total = man.event_count;
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
        std::sort(counts_vec.begin(), counts_vec.end(), [](const auto &a, const auto &b) { return a.second > b.second; });
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
        std::sort(counts_vec.begin(), counts_vec.end(), [](const auto &a, const auto &b) { return a.second > b.second; });
        double others{};
        std::vector<std::pair<std::string, double>> counts_vec_pie;
        for (const auto &[channel, count] : counts_vec) {
            std::cout << channel << ": \t\t\t" << count << std::endl;
            if (count > 0.05) {
                counts_vec_pie.push_back({channel, count});
            } else {
                others += count;
            }
        }
        counts_vec_pie.push_back({"others", others});
        make_pie_plot(counts_vec_pie, output_path.string() + "counts_nonucleon.pdf");
    }
    return 0;
}