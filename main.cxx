#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <algorithm>
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
#include <analysis.h>
#include <vector>
// #define MAX_COUNT 128
constexpr int MAX_COUNT = 128;
using std::string_literals::operator""s;

int main(int argc, char *argv[])
{
    gStyle->SetOptStat(0);
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
    root_chain<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double>
        chain(input_files, "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"});
    double xsec{};
    std::size_t event_count = chain.get_entries();
    std::cout << "event_count = " << event_count << std::endl;

    std::unordered_map<std::string, std::unique_ptr<TH2D>> QWhistograms{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms{}, Whistograms{}, Ehistograms{},
        p_mu{}, pl_mu{}, pt_mu{}, angle_mu{};
    std::unordered_map<std::string, double> xsecs{};
    std::vector<std::pair<std::string, std::string>> single_pion_channels{}, double_pion_channels{};
    int Q_bin_count{80},
        W_bin_count{70},
        E_bin_count{80},
        p_mu_bin_count{80},
        angle_bin{80};
    double Q_min{0}, Q_max{1.8},
        W_min{1.00}, W_max{1.70},
        E_min{1}, E_max{10},
        p_mu_min{0}, p_mu_max{1.5},
        angle_min{0}, angle_max{M_PI};
    for (auto [StdHepN, StdHepPdg, StdHepStatus, StdHepP4, EvtWght] : chain)
    {
        event e;
        for (int i = 0; i < StdHepN; ++i)
        {
            auto pdg = StdHepPdg[i];
            if (StdHepPdg[i] == 1000000010)
            {
                pdg = 2112;
            }

            if (StdHepStatus[i] == 0)
            {
                e.add_particle_in(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
            }
            else if (StdHepStatus[i] == 1)
            {
                e.add_particle_out(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
            }
        }
        xsec += EvtWght / event_count * 1e-38;
        if (e.is_good_event())
        {
            const auto [q2, w] = e.get_q2_w();
            const auto channelname = e.get_event_info();
            if (do_cut)
            {
                if (w > 1.5)
                {
                    continue;
                }
            }
            auto title = channelname;
            title = std::regex_replace(title, std::regex("pi-"), "#pi^{-}");
            title = std::regex_replace(title, std::regex("pi\\+"), "#pi^{+}");
            title = std::regex_replace(title, std::regex("pi0"), "#pi^{0}");
            title = std::regex_replace(title, std::regex("mu-"), "#mu^{-}");
            xsecs[channelname] += EvtWght / event_count;
            if (!QWhistograms[channelname])
            {
                auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); W (GeV)";
                auto thisname = "QW_"s + channelname;
                QWhistograms[channelname] = std::make_unique<TH2D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max, W_bin_count, W_min, W_max);
                if (e.get_pions() == 1) // first time seeing this channel
                {
                    single_pion_channels.push_back(std::make_pair(channelname, title));
                }
                else if (e.get_pions() == 2)
                {
                    double_pion_channels.push_back(std::make_pair(channelname, title));
                }
            }
            QWhistograms[channelname]->Fill(q2, w);
            if (!Qhistograms[channelname])
            {
                auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2}) "s;
                auto thisname = "Q_"s + channelname;
                Qhistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max);
            }
            Qhistograms[channelname]->Fill(q2);
            if (!Whistograms[channelname])
            {
                auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV) "s;
                auto thisname = "W_"s + channelname;
                Whistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), W_bin_count, W_min, W_max);
            }
            Whistograms[channelname]->Fill(w);
            if (!Ehistograms[channelname])
            {
                auto thistitle = title + "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV) "s;
                auto thisname = "E_"s + channelname;
                Ehistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), E_bin_count, E_min, E_max);
            }
            Ehistograms[channelname]->Fill(e.get_enu());
            if (!p_mu[channelname])
            {
                auto thistitle = title + "; #it{p}_{#mu} (GeV); d#it{#sigma}/d#it{p}_{#mu} (cm^{2}/GeV) "s;
                auto thisname = "p_mu_"s + channelname;
                p_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
            }
            p_mu[channelname]->Fill(e.get_p_mu());
            if (!pl_mu[channelname])
            {
                auto thistitle = title + "; #it{pt}_{#mu} (GeV); d#it{#sigma}/d#it{p}_{l#mu} (cm^{2}/GeV) "s;
                auto thisname = "pl_mu_"s + channelname;
                pl_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
            }
            pl_mu[channelname]->Fill(e.get_pl_mu());
            if (!pt_mu[channelname])
            {
                auto thistitle = title + "; #it{pt}_{t#mu} (GeV); d#it{#sigma}/d#it{p}_{t#mu} (cm^{2}/GeV) "s;
                auto thisname = "pt_mu_"s + channelname;
                pt_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
            }
            pt_mu[channelname]->Fill(e.get_pt_mu());
            if (!angle_mu[channelname])
            {
                auto thistitle = title + "; #it{#theta} (rad); d#it{#sigma}/d#it{#theta} (cm^{2}/rad) "s;
                auto thisname = "angle_"s + channelname;
                angle_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), angle_bin, angle_min, angle_max);
            }
            angle_mu[channelname]->Fill(e.get_angle_mu());
        }
    }
    {
        std::fstream outfile(output_prefix + "xsecs.txt", std::ios::out | std::ios::trunc);
        for (const auto &[channelname, xsec] : xsecs)
        {
            outfile << channelname << " " << xsec << std::endl;
        }
    }
    std::cout << "xsec = " << xsec << std::endl;
    std::sort(single_pion_channels.begin(), single_pion_channels.end(),
              [](const auto &lhs, const auto &rhs)
              { return std::hash<std::string>{}(lhs.first) > std::hash<std::string>{}(rhs.first); });
    std::sort(double_pion_channels.begin(), double_pion_channels.end(),
              [](const auto &lhs, const auto &rhs)
              { return std::hash<std::string>{}(lhs.first) > std::hash<std::string>{}(rhs.first); });
    auto Q_bin_width = (Q_max - Q_min) / Q_bin_count;
    auto W_bin_width = (W_max - W_min) / W_bin_count;
    auto pmu_bin_width = (p_mu_max - p_mu_min) / p_mu_bin_count;
    auto angle_bin_width = (angle_max - angle_min) / angle_bin;
    auto E_bin_width = 1. / E_bin_count;

    for (const auto &[channelname, histogram] : QWhistograms)
    {
        histogram->Scale(xsec / event_count / Q_bin_width / W_bin_width);
        plot(histogram, output_prefix, "colz");
        // plot(normalize_slice(histogram), output_prefix, "colz");
        // plot(normalize_slice(histogram, false), output_prefix, "colz");
        // plot_with_normalized(normalize_slice(histogram), std::unique_ptr<TH1D>(reinterpret_cast<TH1D *>(Qhistograms[channelname]->Clone())), output_prefix, "col", true);
        // plot_with_normalized(normalize_slice(histogram, false), std::unique_ptr<TH1D>(reinterpret_cast<TH1D *>(Whistograms[channelname]->Clone())), output_prefix, "col", false);
        normalize_plot(histogram, output_prefix, "col", true);
        normalize_plot(histogram, output_prefix, "col", false);
    }
    for (const auto &[channelname, histogram] : Qhistograms)
    {
        histogram->Scale(xsec / event_count / Q_bin_width);
        // plot(histogram, output_prefix, "hist");
        plot(histogram, output_prefix, xsecs[channelname], "hist");
    }
    for (const auto &[channelname, histogram] : Whistograms)
    {
        histogram->Scale(xsec / event_count / W_bin_width);
        // plot(histogram, output_prefix, "hist");
        plot(histogram, output_prefix, xsecs[channelname], "hist");
    }
    for (const auto &[channelname, histogram] : Ehistograms)
    {
        histogram->Scale(xsec / event_count / E_bin_width);
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : p_mu)
    {
        histogram->Scale(xsec / event_count / pmu_bin_width);
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : pl_mu)
    {
        histogram->Scale(xsec / event_count / pmu_bin_width);
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : pt_mu)
    {
        histogram->Scale(xsec / event_count / pmu_bin_width);
        plot(histogram, output_prefix, "hist");
    }
    for (const auto &[channelname, histogram] : angle_mu)
    {
        histogram->Scale(xsec / event_count / angle_bin_width);
        plot(histogram, output_prefix, "hist");
    }
    {
        std::size_t i{0};
        for (const auto &[channelname, title] : single_pion_channels)
        {
            Qhistograms[channelname]->SetLineColor(col[i]);
            Qhistograms[channelname]->SetFillColor(col[i]);
            Whistograms[channelname]->SetLineColor(col[i]);
            Whistograms[channelname]->SetFillColor(col[i]);
            Ehistograms[channelname]->SetLineColor(col[i]);
            Ehistograms[channelname]->SetFillColor(col[i]);
            ++i;
        }
        for (const auto &[channelname, title] : double_pion_channels)
        {
            Qhistograms[channelname]->SetLineColor(col[i]);
            Qhistograms[channelname]->SetFillColor(col[i]);
            Whistograms[channelname]->SetLineColor(col[i]);
            Whistograms[channelname]->SetFillColor(col[i]);
            Ehistograms[channelname]->SetLineColor(col[i]);
            Ehistograms[channelname]->SetFillColor(col[i]);
            ++i;
        }
    }
    {
        std::unique_ptr<THStack> stack = std::make_unique<THStack>("stackQ", "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2})");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        // std::size_t i{0};
        for (const auto &[channelname, histogram] : Qhistograms)
        {
            // histogram->SetLineColor(col[i]);
            // histogram->SetFillColor(col[i++]);
            stack->Add(histogram.get());
            std::string title(histogram->GetTitle());
            auto x = std::find(title.begin(), title.end(), ';');
            title.erase(x, title.end());
            legend->AddEntry(histogram.get(), title.c_str(), "l");
        }
        plot(stack, output_prefix, "hist", legend.get());
    }
    {
        std::unique_ptr<THStack> stack2 = std::make_unique<THStack>("stackW", "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        // std::size_t i{0};
        for (const auto &[channelname, histogram] : Whistograms)
        {
            // histogram->SetLineColor(col[i]);
            // histogram->SetFillColor(col[i++]);
            stack2->Add(histogram.get());
            std::string title(histogram->GetTitle());
            auto x = std::find(title.begin(), title.end(), ';');
            title.erase(x, title.end());
            legend->AddEntry(histogram.get(), title.c_str(), "l");
        }
        plot(stack2, output_prefix, "hist", legend.get());
    }
    {
        std::unique_ptr<THStack> stack3 = std::make_unique<THStack>("stackE", "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV)");
        std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);
        for (const auto &[channelname, histogram] : Ehistograms)
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
        for (const auto &[channelname, title] : single_pion_channels)
        {
            // std::cout << "adding channel " << channelname << std::endl;
            stackQ->Add(Qhistograms[channelname].get());
            Qhistograms[channelname]->SaveAs(("Q_" + channelname + ".root").c_str());
            stackW->Add(Whistograms[channelname].get());
            Whistograms[channelname]->SaveAs(("W_" + channelname + ".root").c_str());
            stackE->Add(Ehistograms[channelname].get());
            Ehistograms[channelname]->SaveAs(("E_" + channelname + ".root").c_str());
            legendQ->AddEntry(Qhistograms[channelname].get(), title.c_str(), "l");
            legendW->AddEntry(Whistograms[channelname].get(), title.c_str(), "l");
            legendE->AddEntry(Ehistograms[channelname].get(), title.c_str(), "l");
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
        for (const auto &[channelname, title] : double_pion_channels)
        {
            stackQ->Add(Qhistograms[channelname].get());
            stackW->Add(Whistograms[channelname].get());
            stackE->Add(Ehistograms[channelname].get());
            legendQ->AddEntry(Qhistograms[channelname].get(), title.c_str(), "l");
            legendW->AddEntry(Whistograms[channelname].get(), title.c_str(), "l");
            legendE->AddEntry(Ehistograms[channelname].get(), title.c_str(), "l");
        }
        plot(stackQ, output_prefix, "hist", legendQ.get());
        plot(stackW, output_prefix, "hist", legendW.get());
        plot(stackE, output_prefix, "hist", legendE.get());
    }
    return 0;
}