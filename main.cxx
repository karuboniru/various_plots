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

int main(int argc, char *argv[])
{
    gStyle->SetOptStat(0);
    gSystem->ResetSignal(kSigBus);
    gSystem->ResetSignal(kSigSegmentationViolation);
    gSystem->ResetSignal(kSigIllegalInstruction);
    std::array<int, 10> col{kRed, kGreen, kBlue, kMagenta, kCyan,
                            kOrange, kViolet, kGray, kYellow, kBlack};
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <input_file_list> <output_prefix> [do_cut (default 0, 1 to enable)]" << std::endl;
        return 1;
    }
    bool do_cut = argc > 3 ? atoi(argv[3]) : 0;
    const std::string output_prefix = argv[2];
    std::vector<std::string> input_files;
    {
        std::fstream input_file_list_stream(argv[1]);

        std::string input_file;
        while (input_file_list_stream >> input_file)
        {
            // std::cout << "Adding " << input_file << " to input file list" << std::endl;
            input_files.push_back(input_file);
        }
    }
    analysis instance(input_files, do_cut);
    for (const auto &[channelname, histogram] : instance.QWhistograms)
    {
        plot(histogram, output_prefix, "colz");
        // plot(normalize_slice(histogram), output_prefix, "colz");
        // plot(normalize_slice(histogram, false), output_prefix, "colz");
        // plot_with_normalized(normalize_slice(histogram), std::unique_ptr<TH1D>(reinterpret_cast<TH1D *>(Qhistograms[channelname]->Clone())), output_prefix, "col", true);
        // plot_with_normalized(normalize_slice(histogram, false), std::unique_ptr<TH1D>(reinterpret_cast<TH1D *>(Whistograms[channelname]->Clone())), output_prefix, "col", false);
        normalize_plot(histogram, output_prefix, "col", true);
        normalize_plot(histogram, output_prefix, "col", false);
    }
    for (const auto &[channelname, histogram] : instance.p_all.Qhistograms)
    {
        // plot(histogram, output_prefix, "hist");
        plot(histogram, output_prefix, instance.p_all.xsecs[channelname], "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.Whistograms)
    {
        // plot(histogram, output_prefix, "hist");
        plot(histogram, output_prefix, instance.p_all.xsecs[channelname], "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.Ehistograms)
    {
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.p_mu)
    {
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.pl_mu)
    {
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.pt_mu)
    {
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : instance.p_all.angle_mu)
    {
        plot(histogram, output_prefix, "hist");
    }
    {
        std::size_t i{0};
        for (const auto &[channelname, title] : instance.single_pion_channels)
        {
            instance.p_all.Qhistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Qhistograms[channelname]->SetFillColor(col[i]);
            instance.p_all.Whistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Whistograms[channelname]->SetFillColor(col[i]);
            instance.p_all.Ehistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Ehistograms[channelname]->SetFillColor(col[i]);
            ++i;
        }
        for (const auto &[channelname, title] : instance.double_pion_channels)
        {
            instance.p_all.Qhistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Qhistograms[channelname]->SetFillColor(col[i]);
            instance.p_all.Whistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Whistograms[channelname]->SetFillColor(col[i]);
            instance.p_all.Ehistograms[channelname]->SetLineColor(col[i]);
            instance.p_all.Ehistograms[channelname]->SetFillColor(col[i]);
            ++i;
        }
    }
    {
        std::unique_ptr<THStack> stack = std::make_unique<THStack>("stackQ", "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2})");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        // std::size_t i{0};
        for (const auto &[channelname, histogram] : instance.p_all.Qhistograms)
        {
            // histogram->SetLineColor(col[i]);
            // histogram->SetFillColor(col[i++]);
            stack->Add(histogram.get());
            std::string title(histogram->GetTitle());
            auto x = std::find(title.begin(), title.end(), ';');
            title.erase(x, title.end());
            legend->AddEntry(histogram.get(), title.c_str(), "f");
        }
        plot(stack, output_prefix, "hist", legend.get());
    }
    {
        std::unique_ptr<THStack> stack2 = std::make_unique<THStack>("stackW", "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        // std::size_t i{0};
        for (const auto &[channelname, histogram] : instance.p_all.Whistograms)
        {
            // histogram->SetLineColor(col[i]);
            // histogram->SetFillColor(col[i++]);
            stack2->Add(histogram.get());
            std::string title(histogram->GetTitle());
            auto x = std::find(title.begin(), title.end(), ';');
            title.erase(x, title.end());
            legend->AddEntry(histogram.get(), title.c_str(), "f");
        }
        plot(stack2, output_prefix, "hist", legend.get());
    }
    {
        std::unique_ptr<THStack> stack3 = std::make_unique<THStack>("stackE", "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        for (const auto &[channelname, histogram] : instance.p_all.Ehistograms)
        {
            // histogram->SetLineColor(col[i]);
            // histogram->SetFillColor(col[i++]);
            stack3->Add(histogram.get());
            std::string title(histogram->GetTitle());
            auto x = std::find(title.begin(), title.end(), ';');
            title.erase(x, title.end());
            legend->AddEntry(histogram.get(), title.c_str(), "l");
        }
        plot(stack3, output_prefix, "hist", legend.get());
    }
    {
        std::unique_ptr<THStack> stackQ = std::make_unique<THStack>("stackQ1pi", "Stack Cross Section for 1#pi Channels; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2})");
        std::unique_ptr<THStack> stackW = std::make_unique<THStack>("stackW1pi", "Stack Cross Section for 1#pi Channels; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV)");
        std::unique_ptr<THStack> stackE = std::make_unique<THStack>("stackE1pi", "Stack Cross Section for 1#pi Channels; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legendQ = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        std::unique_ptr<TLegend> legendW = std::make_unique<TLegend>(0.1, 0.7, 0.3, 0.9);
        std::unique_ptr<TLegend> legendE = std::make_unique<TLegend>(0.1, 0.1, 0.3, 0.3);
        for (const auto &[channelname, title] : instance.single_pion_channels)
        {
            // std::cout << "adding channel " << channelname << std::endl;
            stackQ->Add(instance.p_all.Qhistograms[channelname].get());
            instance.p_all.Qhistograms[channelname]->SaveAs(("Q_" + channelname + ".root").c_str());
            stackW->Add(instance.p_all.Whistograms[channelname].get());
            instance.p_all.Whistograms[channelname]->SaveAs(("W_" + channelname + ".root").c_str());
            stackE->Add(instance.p_all.Ehistograms[channelname].get());
            instance.p_all.Ehistograms[channelname]->SaveAs(("E_" + channelname + ".root").c_str());
            legendQ->AddEntry(instance.p_all.Qhistograms[channelname].get(), title.c_str(), "f");
            legendW->AddEntry(instance.p_all.Whistograms[channelname].get(), title.c_str(), "f");
            legendE->AddEntry(instance.p_all.Ehistograms[channelname].get(), title.c_str(), "f");
        }
        plot(stackQ, output_prefix, "hist", legendQ.get());
        plot(stackW, output_prefix, "hist", legendW.get());
        plot(stackE, output_prefix, "hist", legendE.get());
    }
    {
        std::unique_ptr<THStack> stackQ = std::make_unique<THStack>("stackQ2pi", "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2})");
        std::unique_ptr<THStack> stackW = std::make_unique<THStack>("stackW2pi", "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV)");
        std::unique_ptr<THStack> stackE = std::make_unique<THStack>("stackE2pi", "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legendQ = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        std::unique_ptr<TLegend> legendW = std::make_unique<TLegend>(0.1, 0.7, 0.3, 0.9);
        std::unique_ptr<TLegend> legendE = std::make_unique<TLegend>(0.1, 0.1, 0.3, 0.3);
        for (const auto &[channelname, title] : instance.double_pion_channels)
        {
            stackQ->Add(instance.p_all.Qhistograms[channelname].get());
            stackW->Add(instance.p_all.Whistograms[channelname].get());
            stackE->Add(instance.p_all.Ehistograms[channelname].get());
            legendQ->AddEntry(instance.p_all.Qhistograms[channelname].get(), title.c_str(), "f");
            legendW->AddEntry(instance.p_all.Whistograms[channelname].get(), title.c_str(), "f");
            legendE->AddEntry(instance.p_all.Ehistograms[channelname].get(), title.c_str(), "f");
        }
        plot(stackQ, output_prefix, "hist", legendQ.get());
        plot(stackW, output_prefix, "hist", legendW.get());
        plot(stackE, output_prefix, "hist", legendE.get());
    }
    return 0;
}