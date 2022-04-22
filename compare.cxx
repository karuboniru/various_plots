#include <TCanvas.h>
#include <TDatabasePDG.h>
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

void plot_compare(TH1D *res2, TH1D *res2_wcut, TH1D *hyb, std::string output_prefix, std::array<double, 3> xsec)
{
    // std::unique_ptr<TCanvas> c1(new TCanvas("c1", "c1", 800, 600));

    auto c1 = getCanvas();
    PadSetup(c1);
    std::unique_ptr<TLegend> leg(new TLegend);
    if (xsec[0])
    {
        leg->AddEntry(res2, ("resevent2, #sigma = "s + std::to_string(xsec[0]) + "x10^{-38} cm^{2}"s).c_str(), "f");
        leg->AddEntry(res2_wcut, ("resevent2, w. cut, #sigma = "s + std::to_string(xsec[1]) + "x10^{-38} cm^{2}"s).c_str(), "f");
        leg->AddEntry(hyb, ("hybrid, #sigma = "s + std::to_string(xsec[2]) + "x10^{-38} cm^{2}"s).c_str(), "f");
    }
    else
    {
        leg->AddEntry(res2, "resevent2", "f");
        leg->AddEntry(res2_wcut, "resevent2_cut", "f");
        leg->AddEntry(hyb, "hybrid", "f");
    }

    double max = std::max({res2->GetMaximum(), res2_wcut->GetMaximum(), hyb->GetMaximum()});
    res2->SetMaximum(max * 1.1);
    res2->SetMinimum(0);
    res2->SetLineColor(kRed);
    res2_wcut->SetLineColor(kBlack);
    hyb->SetLineColor(kBlue);
    ResetStyle(res2, c1->GetPad(0), true);
    res2->Draw("hist");
    res2_wcut->Draw("hist same");
    hyb->Draw("hist same");
    leg->Draw();
    auto name = output_prefix + "compare_"s + res2->GetName();
    c1->SaveAs((name + ".pdf").c_str());
    c1->SaveAs((name + ".png").c_str());
}


int main(int argc, char *argv[])
{
    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigIllegalInstruction);
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    // std::array<int, 10> col{kRed, kGreen, kBlue, kMagenta, kCyan,
    //                         kOrange, kViolet, kGray, kYellow, kBlack};
    if (argc < 4)
    {
        std::cout << "Usage: " << argv[0] << " <input_file_list_res2> <input_file_list_hyb> <output_prefix>" << std::endl;
        return 1;
    }
    const std::string output_prefix = argv[3];
    // const bool do_cut = argc > 4 ? std::stoi(argv[4]) : false;
    std::vector<std::string> input_files, input_files_hyb;
    {
        std::fstream input_file_list_stream(argv[1]);
        std::string input_file;
        while (input_file_list_stream >> input_file)
        {
            input_files.push_back(input_file);
        }
    }
    {
        std::fstream input_file_list_stream(argv[2]);
        std::string input_file;
        while (input_file_list_stream >> input_file)
        {
            input_files_hyb.push_back(input_file);
        }
    }
    analysis resevent2(input_files, true);
    analysis hybrid(input_files_hyb, false);
    for (const auto &[channelname, title] : resevent2.single_pion_channels)
    {
        std::array<double, 3> xsecp{resevent2.p_all.xsecs[channelname], resevent2.p_cut.xsecs[channelname], hybrid.p_all.xsecs[channelname]};
        plot_compare(resevent2.p_all.Qhistograms[channelname].get(), resevent2.p_cut.Qhistograms[channelname].get(), hybrid.p_all.Qhistograms[channelname].get(), output_prefix, xsecp);
        plot_compare(resevent2.p_all.Whistograms[channelname].get(), resevent2.p_cut.Whistograms[channelname].get(), hybrid.p_all.Whistograms[channelname].get(), output_prefix, xsecp);
        plot_compare(resevent2.p_all.Ehistograms[channelname].get(), resevent2.p_cut.Ehistograms[channelname].get(), hybrid.p_all.Ehistograms[channelname].get(), output_prefix, {0, 0, 0.});
        plot_compare(resevent2.p_all.p_mu[channelname].get(), resevent2.p_cut.p_mu[channelname].get(), hybrid.p_all.p_mu[channelname].get(), output_prefix, xsecp);
        plot_compare(resevent2.p_all.pl_mu[channelname].get(), resevent2.p_cut.pl_mu[channelname].get(), hybrid.p_all.pl_mu[channelname].get(), output_prefix, xsecp);
        plot_compare(resevent2.p_all.pt_mu[channelname].get(), resevent2.p_cut.pt_mu[channelname].get(), hybrid.p_all.pt_mu[channelname].get(), output_prefix, xsecp);
        plot_compare(resevent2.p_all.angle_mu[channelname].get(), resevent2.p_cut.angle_mu[channelname].get(), hybrid.p_all.angle_mu[channelname].get(), output_prefix, xsecp);
    }
    return 0;
}