#include <TLorentzVector.h>
#include <analysis.h>

analysis::analysis(std::vector<std::string> files, bool do_cut) {
    static constexpr int MAX_COUNT = 128;

    root_chain<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double> chain(files, "nRooTracker",
                                                                                        {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght"});
    auto event_count = chain.get_entries();
    for (auto [StdHepN, StdHepPdg, StdHepStatus, StdHepP4, EvtWght] : chain) {
        event e;
        for (int i = 0; i < StdHepN; ++i) {
            auto pdg = StdHepPdg[i];
            if (StdHepPdg[i] == 1000000010) {
                pdg = 2112;
            }

            if (StdHepStatus[i] == 0) {
                e.add_particle_in(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
            } else if (StdHepStatus[i] == 2) {
                e.add_particle_out(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
            }
        }
        xsec += EvtWght / event_count * 1e-38;
        if (e.is_good_event()) {
            const auto [q2, w] = e.get_q2_w();
            const auto channelname = e.get_event_info();
            // if (w - 1.6 > 0.01)
            // {
            //     std::cerr << "w = " << w << "@" << channelname << std::endl;
            // }
            auto title = channelname;
            title = std::regex_replace(title, std::regex("pi-"), "#pi^{-}");
            title = std::regex_replace(title, std::regex("pi\\+"), "#pi^{+}");
            title = std::regex_replace(title, std::regex("pi0"), "#pi^{0}");
            title = std::regex_replace(title, std::regex("mu-"), "#mu^{-}");
            if (!QWhistograms[channelname]) {
                auto thistitle = title + "; #it{Q}^{2} (GeV^{2}); W (GeV)";
                auto thisname = "QW_"s + channelname;
                QWhistograms[channelname] = std::make_unique<TH2D>(thisname.c_str(), thistitle.c_str(), Q_bin_count, Q_min, Q_max, W_bin_count, W_min, W_max);
                if (e.get_pions() == 1) { // first time seeing this channel
                    single_pion_channels.push_back(std::make_pair(channelname, title));
                } else if (e.get_pions() == 2) {
                    double_pion_channels.push_back(std::make_pair(channelname, title));
                }
            }
            QWhistograms[channelname]->Fill(q2, w);
            p_all.add_event(q2, w, e, EvtWght / event_count, channelname, title);
            if (do_cut && w < 1.5) { // W cut on 1.5 GeV
                p_cut.add_event(q2, w, e, EvtWght / event_count, channelname, title);
            }
        }
    }
    p_all.finalize_plot(xsec, event_count);
    p_cut.finalize_plot(xsec, event_count);
    std::sort(single_pion_channels.begin(), single_pion_channels.end(),
              [](const auto &lhs, const auto &rhs) { return std::hash<std::string>{}(lhs.first) > std::hash<std::string>{}(rhs.first); });
    std::sort(double_pion_channels.begin(), double_pion_channels.end(),
              [](const auto &lhs, const auto &rhs) { return std::hash<std::string>{}(lhs.first) > std::hash<std::string>{}(rhs.first); });
    auto Q_bin_width = (Q_max - Q_min) / Q_bin_count;
    auto W_bin_width = (W_max - W_min) / W_bin_count;
    // auto pmu_bin_width = (p_mu_max - p_mu_min) / p_mu_bin_count;
    // auto pmu_t_bin_width = (p_mu_t_max - p_mu_t_min) / p_mu_bin_count;
    // auto angle_bin_width = (angle_max - angle_min) / angle_bin;
    // auto E_bin_width = 1. / E_bin_count;
    for (const auto &[channelname, histogram] : QWhistograms) {
        histogram->Scale(xsec / event_count / Q_bin_width / W_bin_width);
    }
}