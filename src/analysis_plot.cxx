#include <analysis.h>

void analysis::plots::finalize_plot(double xsec, std::size_t event_count) {
    auto Q_bin_width = (Q_max - Q_min) / Q_bin_count;
    auto W_bin_width = (W_max - W_min) / W_bin_count;
    auto pmu_bin_width = (p_mu_max - p_mu_min) / p_mu_bin_count;
    auto pmu_t_bin_width = (p_mu_t_max - p_mu_t_min) / p_mu_bin_count;
    auto angle_bin_width = (angle_max - angle_min) / angle_bin;
    auto E_bin_width = 1. / E_bin_count;
    for (const auto &[channelname, histogram] : Qhistograms) {
        histogram->Scale(xsec / event_count / Q_bin_width);
    }
    for (const auto &[channelname, histogram] : Whistograms) {
        histogram->Scale(xsec / event_count / W_bin_width);
    }
    for (const auto &[channelname, histogram] : Ehistograms) {
        histogram->Scale(xsec / event_count / E_bin_width);
    }
    for (const auto &[channelname, histogram] : p_mu) {
        histogram->Scale(xsec / event_count / pmu_bin_width);
    }
    for (const auto &[channelname, histogram] : pl_mu) {
        histogram->Scale(xsec / event_count / pmu_bin_width);
    }
    for (const auto &[channelname, histogram] : pt_mu) {
        histogram->Scale(xsec / event_count / pmu_t_bin_width);
    }
    for (const auto &[channelname, histogram] : angle_mu) {
        histogram->Scale(xsec / event_count / angle_bin_width);
    }
}

void analysis::plots::add_event(const double q2, const double w, const event &e, double dxsec, std::string channelname, std::string title) {
    xsecs[channelname] += dxsec;
    if (!Qhistograms[channelname]) {
        auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); d#it{#sigma}/d#it{Q}^{2} (cm^{2}/GeV^{2}) "s;
        auto thisname = "Q_"s + channelname;
        Qhistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max);
    }
    Qhistograms[channelname]->Fill(q2);
    if (!Whistograms[channelname]) {
        auto thistitle = title + "; #it{W} (GeV); d#it{#sigma}/d#it{W} (cm^{2}/GeV) "s;
        auto thisname = "W_"s + channelname;
        Whistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), W_bin_count, W_min, W_max);
    }
    Whistograms[channelname]->Fill(w);
    if (!Ehistograms[channelname]) {
        auto thistitle = title + "; #it{E} (GeV); #it{#sigma}(#it{E}) (cm^{2}) "s;
        auto thisname = "E_"s + channelname;
        Ehistograms[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), E_bin_count, E_min, E_max);
    }
    Ehistograms[channelname]->Fill(e.get_enu());
    if (!p_mu[channelname]) {
        auto thistitle = title + "; #it{p}_{#mu} (GeV/#it{c}); d#it{#sigma}/d#it{p}_{#mu} (cm^{2}#it{c}/GeV) "s;
        auto thisname = "p_mu_"s + channelname;
        p_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
    }
    p_mu[channelname]->Fill(e.get_p_mu());
    if (!pl_mu[channelname]) {
        auto thistitle = title + "; #it{p}_{l,#mu} (GeV/#it{c}); d#it{#sigma}/d#it{p}_{l,#mu} (cm^{2}#it{c}/GeV) "s;
        auto thisname = "pl_mu_"s + channelname;
        pl_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_min, p_mu_max);
    }
    pl_mu[channelname]->Fill(e.get_pl_mu());
    if (!pt_mu[channelname]) {
        auto thistitle = title + "; #it{p}_{t,#mu} (GeV/#it{c}); d#it{#sigma}/d#it{p}_{t,#mu} (cm^{2}#it{c}/GeV) "s;
        auto thisname = "pt_mu_"s + channelname;
        pt_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), p_mu_bin_count, p_mu_t_min, p_mu_t_max);
    }
    pt_mu[channelname]->Fill(e.get_pt_mu());
    if (!angle_mu[channelname]) {
        auto thistitle = title + "; #it{#theta} (rad); d#it{#sigma}/d#it{#theta} (cm^{2}/rad) "s;
        auto thisname = "angle_"s + channelname;
        angle_mu[channelname] = std::make_unique<TH1D>(thisname.c_str(), thistitle.c_str(), angle_bin, angle_min, angle_max);
    }
    angle_mu[channelname]->Fill(e.get_angle_mu());
}