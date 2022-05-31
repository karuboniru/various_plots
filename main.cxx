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
    analysis::plots instance;
    if (config["type"] == "genie") {
        instance = run_manager_genie::run_analysis(config).plot;
    } else if (config["type"] == "nuwro") {
        instance = run_manager_nuwro::run_analysis(config).plot;
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
    std::cout << "Writing to " << output_path << std::endl;
    instance.save(output_path);

    return 0;
}