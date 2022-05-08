#pragma once

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
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
#include <deque>

using std::string_literals::operator""s;

class analysis {
public:
    double xsec{};
    std::deque<std::pair<std::string, std::string>> single_pion_channels{}, double_pion_channels{};
    std::unordered_map<std::string, std::unique_ptr<TH2D>> QWhistograms{};
    constexpr static const int Q_bin_count{80}, W_bin_count{70}, E_bin_count{80}, p_mu_bin_count{80}, angle_bin{80};
    constexpr static double Q_min{0}, Q_max{1.8}, W_min{1.00}, W_max{1.70}, E_min{0}, E_max{10}, p_mu_min{0}, p_mu_max{1.5}, p_mu_t_min{0}, p_mu_t_max{0.8},
        angle_min{0}, angle_max{M_PI};
    class plots {
    public:
        std::unordered_map<std::string, std::unique_ptr<TH1D>> Qhistograms{}, Whistograms{}, Ehistograms{}, p_mu{}, pl_mu{}, pt_mu{}, angle_mu{};
        std::unordered_map<std::string, double> xsecs{};
        void add_event(const double q2, const double w, const event &e, const double dxsec, const std::string channelname, const std::string title);
        void finalize_plot(const double xsec, const std::size_t event_count);
    } p_all, p_cut;

public:
    analysis(const std::vector<std::string> files, bool do_cut = false);
};
