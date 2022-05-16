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
#include <regex>
#include <string>
#include <tools.h>
#include <unordered_map>
#include <vector>
// #define MAX_COUNT 128
constexpr int MAX_COUNT = 128;
using std::string_literals::operator""s;

int main(int argc, char *argv[]) {
    TH1::AddDirectory(kFALSE);
    ROOT::EnableThreadSafety();
    gStyle->SetOptStat(0);
    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigIllegalInstruction);
    // std::array<int, 10> col{kRed, kGreen, kBlue, kMagenta, kCyan, kOrange, kViolet, kGray, kYellow, kBlack};
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input_file_list> <output_prefix>" << std::endl;
        return 1;
    }
    const std::string output_prefix = argv[2];
    std::vector<std::string> input_files;
    {
        std::fstream input_file_list_stream(argv[1]);

        std::string input_file;
        while (input_file_list_stream >> input_file) {
            input_files.push_back(input_file);
        }
    }
    // analysis instance(input_files, do_cut);
    chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double> chain(
        input_files, "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"}, std::thread::hardware_concurrency());
    auto instance = chain.run<run_manager>(chain.get_entries());
    auto &protonE = instance.p_all.protonE;
    auto &protonP = instance.p_all.protonP;
    auto &enu = instance.p_all.enu;
    auto &protonE_nocut = instance.p_all.protonE_nocut;
    auto &protonP_nocut = instance.p_all.protonP_nocut;
    TFile f((output_prefix + ".root").c_str(), "RECREATE");
    // protonE.SaveAs((output_prefix + "_protonE.root").c_str());
    // protonP.SaveAs((output_prefix + "_protonP.root").c_str());
    // enu.SaveAs((output_prefix + "_enu.root").c_str());
    protonE.SetDirectory(&f);
    protonP.SetDirectory(&f);
    enu.SetDirectory(&f);
    protonE_nocut.SetDirectory(&f);
    protonP_nocut.SetDirectory(&f);
    f.Write();
    protonE.SetDirectory(0);
    protonP.SetDirectory(0);
    enu.SetDirectory(0);
    protonE_nocut.SetDirectory(0);
    protonP_nocut.SetDirectory(0);
    f.Close();

    return 0;
}