#pragma once

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

class analysis
{
    // private:
public:
    std::unordered_map<std::string, std::unique_ptr<TH2D>> QWhistograms{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms{}, Whistograms{}, Ehistograms{},
        p_mu{}, pl_mu{}, pt_mu{}, angle_mu{};
    std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms_cut{}, Whistograms_cut{}, Ehistograms_cut{},
        p_mu_cut{}, pl_mu_cut{}, pt_mu_cut{}, angle_mu_cut{};
    std::unordered_map<std::string, double> xsecs{}, xsecs_cut{};
    double xsec{};
    std::vector<std::pair<std::string, std::string>> single_pion_channels{}, double_pion_channels{};
    root_chain<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double> chain;
    const int Q_bin_count{80},
        W_bin_count{70},
        E_bin_count{80},
        p_mu_bin_count{80},
        angle_bin{80};
    const double Q_min{0}, Q_max{1.8},
        W_min{1.00}, W_max{1.70},
        E_min{1}, E_max{10},
        p_mu_min{0}, p_mu_max{1.5},
        angle_min{0}, angle_max{M_PI};

public:
    analysis(const std::vector<std::string> files, bool do_cut = false) : chain(files, "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"})
    {
        auto event_count = chain.get_entries();
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
                if (do_cut)
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
                        auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV) "s;
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
                    if (!p_mu_cut[channelname])
                    {
                        auto thistitle = title + "; #it{p}_{#mu} (GeV); d#it{#sigma}/d#it{p}_{#mu} (cm^{2}/GeV) "s;
                        auto thisname = "p_mu_"s + channelname;
                        p_mu_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
                    }
                    p_mu[channelname]->Fill(e.get_p_mu());
                    if (!pl_mu_cut[channelname])
                    {
                        auto thistitle = title + "; #it{pt}_{#mu} (GeV); d#it{#sigma}/d#it{p}_{l#mu} (cm^{2}/GeV) "s;
                        auto thisname = "pl_mu_"s + channelname;
                        pl_mu_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
                    }
                    pl_mu_cut[channelname]->Fill(e.get_pl_mu());
                    if (!pt_mu_cut[channelname])
                    {
                        auto thistitle = title + "; #it{pt}_{t#mu} (GeV); d#it{#sigma}/d#it{p}_{t#mu} (cm^{2}/GeV) "s;
                        auto thisname = "pt_mu_"s + channelname;
                        pt_mu_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
                    }
                    pt_mu_cut[channelname]->Fill(e.get_pt_mu());
                    if (!angle_mu_cut[channelname])
                    {
                        auto thistitle = title + "; #it{#theta} (rad); d#it{#sigma}/d#it{#theta} (cm^{2}/rad) "s;
                        auto thisname = "angle_"s + channelname;
                        angle_mu_cut[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), angle_bin, angle_min, angle_max);
                    }
                    angle_mu_cut[channelname]->Fill(e.get_angle_mu());
                }
            }
        }
        auto Q_bin_width = (Q_max - Q_min) / Q_bin_count;
        auto W_bin_width = (W_max - W_min) / W_bin_count;
        auto pmu_bin_width = (p_mu_max - p_mu_min) / p_mu_bin_count;
        auto angle_bin_width = (angle_max - angle_min) / angle_bin;
        auto E_bin_width = 1. / E_bin_count;
        for (const auto &[channelname, histogram] : QWhistograms)
        {
            histogram->Scale(xsec / event_count / Q_bin_width / W_bin_width);
        }
        for (const auto &[channelname, histogram] : Qhistograms)
        {
            histogram->Scale(xsec / event_count / Q_bin_width);
        }
        for (const auto &[channelname, histogram] : Whistograms)
        {
            histogram->Scale(xsec / event_count / W_bin_width);
        }
        for (const auto &[channelname, histogram] : Ehistograms)
        {
            histogram->Scale(xsec / event_count / E_bin_width);
        }
        for (const auto &[channelname, histogram] : p_mu)
        {
            histogram->Scale(xsec / event_count / pmu_bin_width);
        }
        for (const auto &[channelname, histogram] : pl_mu)
        {
            histogram->Scale(xsec / event_count / pmu_bin_width);
        }
        for (const auto &[channelname, histogram] : pt_mu)
        {
            histogram->Scale(xsec / event_count / pmu_bin_width);
        }
        for (const auto &[channelname, histogram] : angle_mu)
        {
            histogram->Scale(xsec / event_count / angle_bin_width);
        }
        if (do_cut)
        {
            for (const auto &[channelname, histogram] : Qhistograms_cut)
            {
                histogram->Scale(xsec / event_count / Q_bin_width);
            }
            for (const auto &[channelname, histogram] : Whistograms_cut)
            {
                histogram->Scale(xsec / event_count / W_bin_width);
            }
            for (const auto &[channelname, histogram] : Ehistograms_cut)
            {
                histogram->Scale(xsec / event_count / E_bin_width);
            }
            for (const auto &[channelname, histogram] : p_mu_cut)
            {
                histogram->Scale(xsec / event_count / pmu_bin_width);
            }
            for (const auto &[channelname, histogram] : pl_mu_cut)
            {
                histogram->Scale(xsec / event_count / pmu_bin_width);
            }
            for (const auto &[channelname, histogram] : pt_mu_cut)
            {
                histogram->Scale(xsec / event_count / pmu_bin_width);
            }
            for (const auto &[channelname, histogram] : angle_mu_cut)
            {
                histogram->Scale(xsec / event_count / angle_bin_width);
            }
        }
    }
};
