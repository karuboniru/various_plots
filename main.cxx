#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <analysis.h>
#include <chain_helper.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <TObjString.h>
#include <regex>
#include <string>
#include <tools.h>
#include <unordered_map>
#include <vector>
// #define MAX_COUNT 128
constexpr int MAX_COUNT = 128;
using std::string_literals::operator""s;
template <typename T> std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
    TFile root_file{file_path.c_str(), "READ"};
    auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
    assert(objptr);
    return std::unique_ptr<T>{objptr};
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
    const std::string output_prefix = config["output_prefix"];
    std::vector<std::string> input_files;
    {
        // std::fstream input_file_list_stream(config["input_file_list"].get<std::string>());

        // std::string input_file;
        for (auto & file:config["input_file_list"]) {
            input_files.push_back(file);
        }
    }
    auto spline_file = get_object<TGraph>(config["spline_file"], config["spline_path"]);
    // analysis instance(input_files, do_cut);
    chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], TObjString> chain(
        input_files, "gRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"}, std::thread::hardware_concurrency());
    auto instance = chain.run<run_manager>(chain.get_entries(), spline_file.get());
    auto &protonE = instance.plot.protonE;
    auto &protonP = instance.plot.protonP;
    auto &enu = instance.plot.enu;
    auto &protonE_nocut = instance.plot.protonE_nocut;
    auto &protonP_nocut = instance.plot.protonP_nocut;
    auto &leadingP = instance.plot.leadingP;
    auto &leadingP_nocut = instance.plot.leadingP_nocut;
    TFile f((output_prefix + ".root").c_str(), "RECREATE");
    // protonE.SaveAs((output_prefix + "_protonE.root").c_str());
    // protonP.SaveAs((output_prefix + "_protonP.root").c_str());
    // enu.SaveAs((output_prefix + "_enu.root").c_str());
    protonE.SetDirectory(&f);
    protonP.SetDirectory(&f);
    enu.SetDirectory(&f);
    protonE_nocut.SetDirectory(&f);
    protonP_nocut.SetDirectory(&f);
    leadingP.SetDirectory(&f);
    leadingP_nocut.SetDirectory(&f);
    f.Write();
    protonE.SetDirectory(0);
    protonP.SetDirectory(0);
    enu.SetDirectory(0);
    protonE_nocut.SetDirectory(0);
    protonP_nocut.SetDirectory(0);
    leadingP.SetDirectory(0);
    leadingP_nocut.SetDirectory(0);
    f.Close();

    return 0;
}