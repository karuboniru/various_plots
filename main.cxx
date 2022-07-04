#include <TFile.h>
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
        {{"qe", [](event &e) { return e.get_mode() == event::channel::QE; }},
         {"res", [](event &e) { return e.get_mode() == event::channel::RES; }},
         {"dis", [](event &e) { return e.get_mode() == event::channel::DIS; }},
         {"2p2h", [](event &e) { return e.get_mode() == event::channel::MEC; }}},
        {{"0pi", [](event &e) { return e.count_out(211) == 0 && e.count_out(-211) == 0 && e.count_out(111) == 0; }},
         {"1pip", [](event &e) { return e.count_out(211) == 1 && e.count_out(111) == 0 && e.count_out(-211) == 0; }},
         {"1pim", [](event &e) { return e.count_out(-211) == 1 && e.count_out(111) == 0 && e.count_out(211) == 0; }},
         {"1pi0", [](event &e) { return e.count_out(211) == 0 && e.count_out(111) == 1 && e.count_out(-211) == 0; }},
         {"1pi", [](event &e) { return e.count_out(211) + e.count_out(-211) + e.count_out(111) == 1; }},
         {"mpi", [](event &e) { return e.count_out(211) > 0 && e.count_out(111) > 0; }}},
        {{"with_spe", [](event &e) { return e.get_mode() == event::channel::QE && e.get_particle_nofsi().size() == 3; }},
         {"without_spe", [](event &e) { return e.get_mode() == event::channel::QE && e.get_particle_nofsi().size() != 3; }}}};
    // std::initializer_list<analysis::plots::cut_t> cuts{{"mucut", [](event &e) { return e.get_leading_out(14).P() > 1.0; }}};
    std::initializer_list<analysis::plots::cut_t> cuts{};
    std::initializer_list<analysis::plots::analysis_funcs_t> analysis_funcs{{"protonP",
                                                                             [](event &e, TH1 *h) {
                                                                                 for (const auto &[_, p] : e.get_particle_out(2212)) {
                                                                                     h->Fill(p.P(), e.get_weight());
                                                                                 }
                                                                             }}, // proton momentum (particle by particle)
                                                                            {"leadingP",
                                                                             [](event &e, TH1 *h) {
                                                                                 if (e.count_particle_out(2212)) {
                                                                                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                                                                                 }
                                                                             }}, // leading proton momentum (event by event), no cut
                                                                            {"muonP",
                                                                             [](event &e, TH1 *h) {
                                                                                 if (e.count_particle_out(13)) {
                                                                                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                                                                                 }
                                                                             }},
                                                                            {"sum_of_ke_P",
                                                                             [](event &e, TH1 *h) {
                                                                                 if (e.count_particle_out(2212)) {
                                                                                     double p_sum{};
                                                                                     for (const auto &[_, p] : e.get_particle_out(2212)) {
                                                                                         p_sum += p.E() - p.M();
                                                                                     }
                                                                                     h->Fill(p_sum, e.get_weight());
                                                                                 }
                                                                             }},
                                                                            {"protonP_nofsi",
                                                                             [](event &e, TH1 *h) {
                                                                                 for (const auto &[_, p] : e.get_particle_nofsi(2212)) {
                                                                                     if (p.P() != 0)
                                                                                         h->Fill(p.P(), e.get_weight());
                                                                                 }
                                                                             }},
                                                                            {"leadingP_nofsi",
                                                                             [](event &e, TH1 *h) {
                                                                                 if (e.count_particle_nofsi(2212)) {
                                                                                     h->Fill(e.get_leading_nofsi(2212).P(), e.get_weight());
                                                                                 }
                                                                             }},
                                                                            {"enu", [](event &e, TH1 *h) { h->Fill(e.get_enu(), e.get_weight()); }},
                                                                            {"kaonE", [](event &e, TH1 *h) {
                                                                                 for (const auto &[_, p] : e.get_particle_out(321)) {
                                                                                     h->Fill(p.E() - p.M(), e.get_weight());
                                                                                 }
                                                                             }}};
    if (config["type"] == "genie") {
        run_manager_genie::run_analysis(config, analysis_funcs, cut_group, cuts).plot.save(output_path);
    } else if (config["type"] == "nuwro") {
        run_manager_nuwro::run_analysis(config, analysis_funcs, cut_group, cuts).plot.save(output_path);
    } else if (config["type"] == "gibuu") {
        run_manager_gibuu::run_analysis(config, analysis_funcs, cut_group, cuts).plot.save(output_path);
    } else {
        std::cerr << "Unknown type: " << config["type"] << std::endl;
        return 1;
    }
    // TFile f(output_path.c_str(), "RECREATE");
    // for (auto &[name, hist] : instance.histos) {
    //     hist->SetDirectory(&f);
    // }
    // f.Write();
    // for (auto &[name, hist] : instance.histos) {
    //     hist->SetDirectory(0);
    // }
    // f.Close();
    // std::cout << "Writing to " << output_path << std::endl;
    // instance.save(output_path);

    return 0;
}