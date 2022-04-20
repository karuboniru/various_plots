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
#include <vector>
// #define MAX_COUNT 128
constexpr int MAX_COUNT = 128;
using std::string_literals::operator""s;

void plot_compare(TH1D *res2, TH1D *hyb, std::string output_prefix, std::pair<double, double> xsec)
{
    // std::unique_ptr<TCanvas> c1(new TCanvas("c1", "c1", 800, 600));
    auto c1 = getCanvas();
    PadSetup(c1);
    std::unique_ptr<TLegend> leg(new TLegend);
    if (xsec.first)
    {
        leg->AddEntry(res2, ("resevent2, #sigma = "s + std::to_string(xsec.first) + "x10^{-38} cm^{2}"s).c_str());
        leg->AddEntry(hyb, ("hybrid, #sigma = "s + std::to_string(xsec.second) + "x10^{-38} cm^{2}"s).c_str());
    }
    else
    {
        leg->AddEntry(res2, "resevent2");
        leg->AddEntry(hyb, "hybrid");
    }

    double max = std::max(res2->GetMaximum(), hyb->GetMaximum());
    res2->SetMaximum(max * 1.1);
    res2->SetMinimum(0);
    res2->SetLineColor(kRed);
    hyb->SetLineColor(kBlue);
    ResetStyle(res2, c1->GetPad(0), true);
    res2->Draw("hist");
    hyb->Draw("hist same");
    leg->Draw();
    auto name = output_prefix + "compare_"s + res2->GetName();
    c1->SaveAs((name + ".pdf").c_str());
    c1->SaveAs((name + ".png").c_str());
}

void plot_compare(TH1D *res2, TH1D *res2_cut, TH1D *hyb, std::string output_prefix, std::array<double, 3> xsec)
{
    // std::unique_ptr<TCanvas> c1(new TCanvas("c1", "c1", 800, 600));
    auto c1 = getCanvas();
    PadSetup(c1);
    std::unique_ptr<TLegend> leg(new TLegend);
    if (xsec[0])
    {
        leg->AddEntry(res2, ("resevent2, #sigma = "s + std::to_string(xsec[0]) + "x10^{-38} cm^{2}"s).c_str());
        leg->AddEntry(res2_cut, ("resevent2, w/o. cut, #sigma = "s + std::to_string(xsec[1]) + "x10^{-38} cm^{2}"s).c_str());
        leg->AddEntry(hyb, ("hybrid, #sigma = "s + std::to_string(xsec[2]) + "x10^{-38} cm^{2}"s).c_str());
    }
    else
    {
        leg->AddEntry(res2, "resevent2");
        leg->AddEntry(res2, "resevent2_cut");
        leg->AddEntry(hyb, "hybrid");
    }

    double max = std::max(res2->GetMaximum(), res2_cut->GetMaximum(), hyb->GetMaximum());
    res2->SetMaximum(max * 1.1);
    res2->SetMinimum(0);
    res2->SetLineColor(kRed);
    hyb->SetLineColor(kBlue);
    ResetStyle(res2, c1->GetPad(0), true);
    res2->Draw("hist");
    res2_cut->Draw("hist same");
    hyb->Draw("hist same");
    leg->Draw();
    auto name = output_prefix + "compare_"s + res2->GetName();
    c1->SaveAs((name + ".pdf").c_str());
    c1->SaveAs((name + ".png").c_str());
}

int main(int argc, char *argv[])
{
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

    double xsec{}, xsec_hyb{}, xsec_cut{};
    std::size_t event_count{}, event_count_hyb{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms{}, Whistograms{}, Ehistograms{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms_cut{}, Whistograms_cut{}, Ehistograms_cut{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms_hyb{}, Whistograms_hyb{}, Ehistograms_hyb{};
    std::vector<std::pair<std::string, std::string>> single_pion_channels{};
    std::unordered_map<std::string, double> xsecs{}, xsecs_hyb{}, xsecs_cut{};
    int Q_bin_count{80}, W_bin_count{70}, E_bin_count{80};
    double Q_min{0}, Q_max{1.8}, W_min{1.0}, W_max{1.7}, E_min{0}, E_max{10};
    {
        root_chain<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double>
            chain(input_files, "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"});
        event_count = chain.get_entries();
        std::cout << "event_count = " << event_count << std::endl;
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
            if (e.get_pions() == 1 && e.is_good_event())
            {
                const auto [q2, w] = e.get_q2_w();
                const auto channelname = e.get_event_info();
                auto title = channelname;
                title = std::regex_replace(title, std::regex("pi-"), "#pi^{-}");
                title = std::regex_replace(title, std::regex("pi\\+"), "#pi^{+}");
                title = std::regex_replace(title, std::regex("pi0"), "#pi^{0}");
                title = std::regex_replace(title, std::regex("mu-"), "#mu^{-}");
                xsecs[channelname] += EvtWght / event_count;
                if (!Qhistograms[channelname])
                {
                    single_pion_channels.push_back(std::make_pair(channelname, title));
                    auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2}) "s;
                    auto thisname = "Q_"s + channelname;
                    Qhistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max);
                }
                Qhistograms[channelname]->Fill(q2);
                if (!Whistograms[channelname])
                {
                    auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/dW (cm^{2}/GeV) "s;
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
                if (w < 1.5) // The W < 1.5 GeV Part, for comparing with hybrid data
                {
                    xsecs_cut[channelname] += EvtWght / event_count;
                    if (!Qhistograms_cut[channelname])
                    {
                        auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2}) "s;
                        auto thisname = "Q_"s + channelname;
                        Qhistograms_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max);
                    }
                    Qhistograms_cut[channelname]->Fill(q2);
                    if (!Whistograms_cut[channelname])
                    {
                        auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/dW (cm^{2}/GeV) "s;
                        auto thisname = "W_"s + channelname;
                        Whistograms_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), W_bin_count, W_min, W_max);
                    }
                    Whistograms_cut[channelname]->Fill(w);
                    if (!Ehistograms_cut[channelname])
                    {
                        auto thistitle = title + "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV) "s;
                        auto thisname = "E_"s + channelname;
                        Ehistograms_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), E_bin_count, E_min, E_max);
                    }
                    Ehistograms_cut[channelname]->Fill(e.get_enu());
                }
            }
        }
    }
    {
        root_chain<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double>
            chain_hyb(input_files_hyb, "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"});
        event_count_hyb = chain_hyb.get_entries();
        std::cout << "event_count_hyb = " << event_count_hyb << std::endl;
        for (auto [StdHepN, StdHepPdg, StdHepStatus, StdHepP4, EvtWght] : chain_hyb)
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
            xsec_hyb += EvtWght / event_count_hyb * 1e-38;
            if (e.get_pions() == 1 && e.is_good_event())
            {

                const auto [q2, w] = e.get_q2_w();
                const auto channelname = e.get_event_info();
                auto title = channelname;
                title = std::regex_replace(title, std::regex("pi-"), "#pi^{-}");
                title = std::regex_replace(title, std::regex("pi\\+"), "#pi^{+}");
                title = std::regex_replace(title, std::regex("pi0"), "#pi^{0}");
                title = std::regex_replace(title, std::regex("mu-"), "#mu^{-}");
                xsecs_hyb[channelname] += EvtWght / event_count;

                if (!Qhistograms_hyb[channelname])
                {
                    auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2}) "s;
                    auto thisname = "Q_"s + channelname;
                    Qhistograms_hyb[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max);
                }
                Qhistograms_hyb[channelname]->Fill(q2);
                if (!Whistograms_hyb[channelname])
                {
                    auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/dW (cm^{2}/GeV) "s;
                    auto thisname = "W_"s + channelname;
                    Whistograms_hyb[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), W_bin_count, W_min, W_max);
                }
                Whistograms_hyb[channelname]->Fill(w);
                if (!Ehistograms_hyb[channelname])
                {
                    auto thistitle = title + "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}/GeV) "s;
                    auto thisname = "E_"s + channelname;
                    Ehistograms_hyb[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), E_bin_count, E_min, E_max);
                }
                Ehistograms_hyb[channelname]->Fill(e.get_enu());
            }
        }
    }
    std::cout << "xsec = " << xsec << std::endl;
    std::cout << "xsec_hyb = " << xsec_hyb << std::endl;
    std::sort(single_pion_channels.begin(), single_pion_channels.end(),
              [](const auto &lhs, const auto &rhs)
              { return std::hash<std::string>{}(lhs.first) > std::hash<std::string>{}(rhs.first); });

    auto Q_bin_width = (Q_max - Q_min) / Q_bin_count;
    auto W_bin_width = (W_max - W_min) / W_bin_count;
    auto E_bin_width = 1. / E_bin_count;

    for (const auto &[channelname, title] : single_pion_channels)
    {
        Qhistograms[channelname]->Scale(xsec / event_count / Q_bin_width);
        Whistograms[channelname]->Scale(xsec / event_count / W_bin_width);
        Ehistograms[channelname]->Scale(xsec / event_count / E_bin_width);
        Qhistograms_hyb[channelname]->Scale(xsec_hyb / event_count_hyb / Q_bin_width);
        Whistograms_hyb[channelname]->Scale(xsec_hyb / event_count_hyb / W_bin_width);
        Ehistograms_hyb[channelname]->Scale(xsec_hyb / event_count_hyb / E_bin_width);
        auto xsecp = std::array<double, 3>{xsecs[channelname], xsecs_cut[channelname], xsecs_hyb[channelname]};
        plot_compare(Qhistograms[channelname].get(), Qhistograms_cut[channelname].get(), Qhistograms_hyb[channelname].get(), output_prefix, xsecp);
        plot_compare(Whistograms[channelname].get(), Qhistograms_cut[channelname].get(), Whistograms_hyb[channelname].get(), output_prefix, xsecp);
        plot_compare(Ehistograms[channelname].get(), Qhistograms_cut[channelname].get(), Ehistograms_hyb[channelname].get(), output_prefix, {0, 0, 0});
    }

    return 0;
}