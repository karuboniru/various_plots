#pragma once

#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TStyle.h>
#include <chain_helper.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <unordered_set>
#include <string>
#include <unordered_map>
#include <vector>

using std::string_literals::operator""s;

double get_xsec(TH1 *h_rate, TGraph *spline) {
    double fluxint{};
    TF1 func(
        "spline", [&spline](double *x, double *) { return spline->Eval(*x); }, 0, h_rate->GetXaxis()->GetXmax(), 0);
    for (int ii = 1; ii <= h_rate->GetNbinsX(); ii++) {
        double bin_c = h_rate->GetBinContent(ii);
        double bin_up = h_rate->GetXaxis()->GetBinUpEdge(ii);
        double bin_low = h_rate->GetXaxis()->GetBinLowEdge(ii);
        double bin_width = bin_up - bin_low;
        if (bin_c < 1 || func.Integral(bin_low, bin_up) == 0) {
            continue;
        }

        fluxint += bin_c / func.Integral(bin_low, bin_up) * bin_width;
    }
    return h_rate->Integral() / fluxint;
}

template <typename T> std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
    TFile root_file{file_path.c_str(), "READ"};
    auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
    assert(objptr);
    return std::unique_ptr<T>{objptr};
}

class analysis {
public:
    constexpr static double pmin = 0, pmax = 20., bin_wid = 4e-3;
    constexpr static size_t bin_count = (pmax - pmin) / bin_wid;
    class plots {
    public:
        // TH1D enu;
        std::unordered_set<std::string> plot_names{};
        std::unordered_map<std::string, std::unique_ptr<TH1D>> histos;
        std::unordered_map<std::string, std::function<void(event &, TH1D *)>> plot_funcs{
            {"protonP_nocut",
             [](event &e, TH1D *h) {
                 for (const auto &[_, p] : e.get_particle_out(2212)) {
                     h->Fill(p.P(), e.get_weight());
                 }
             }}, // proton momentum (particle by particle)
            {"protonP",
             [](event &e, TH1D *h) {
                 if (e.TKI_mu_p_cut()) {
                     for (const auto &[_, p] : e.get_particle_out(2212)) {
                         h->Fill(p.P(), e.get_weight());
                     }
                 }
             }}, // proton momentum, with muon/p angle and momentum cut (particle by particle)
            {"leadingP_nocut",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(2212)) {
                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                 }
             }}, // leading proton momentum (event by event), no cut
            {"leadingP",
             [](event &e, TH1D *h) {
                 if (e.TKI_mu_p_cut() && e.count_particle_out(2212)) {
                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                 }
             }}, // leading proton momentum (event by event), with muon/p angle and momentum cut
            {"leadingP_mucut",
             [](event &e, TH1D *h) {
                 if (e.TKI_mu_cut() && e.count_particle_out(2212)) {
                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                 }
             }}, // leading proton momentum (event by event), with muon angle and momentum cut
            {"protonP_mucut",
             [](event &e, TH1D *h) {
                 if (e.TKI_mu_cut()) {
                     for (const auto &[_, p] : e.get_particle_out(2212)) {
                         h->Fill(p.P(), e.get_weight());
                     }
                 }
             }}, // proton momentum, with muon angle and momentum cut (particle by particle)
            {"leadingP_0pi",
             [](event &e, TH1D *h) {
                 if (e.TKI_mu_cut() && e.count_particle_out(111) == 0 && e.count_particle_out(2212)) {
                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                 }
             }}, // leading proton momentum (event by event), with muon angle and momentum cut, no pi0
            {"muon",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13)) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_qe",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::QE) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_qe_forward",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::QE && e.get_leading_out(13).Theta() < 20. / 180 * M_PI) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_res",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::RES) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_res_forward",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::RES && e.get_leading_out(13).Theta() < 20. / 180 * M_PI) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_dis",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::DIS) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_dis_forward",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::DIS && e.get_leading_out(13).Theta() < 20. / 180 * M_PI) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_2p2h",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::MEC) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"muon_2p2h_forward",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13) && e.get_mode() == event::channel::MEC && e.get_leading_out(13).Theta() < 20. / 180 * M_PI) {
                     h->Fill(e.get_leading_out(13).P(), e.get_weight());
                 }
             }},
            {"sum_of_ke_P",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(2212)) {
                     double p_sum{};
                     for (const auto &[_, p] : e.get_particle_out(2212)) {
                         p_sum += p.E() - p.M();
                     }
                     h->Fill(p_sum, e.get_weight());
                 }
             }},
            {"protonP_nofsi",
             [](event &e, TH1D *h) {
                 for (const auto &[_, p] : e.get_particle_nofsi(2212)) {
                     h->Fill(p.P(), e.get_weight());
                 }
             }},
            {"leadingP_nofsi",
             [](event &e, TH1D *h) {
                 if (e.count_particle_nofsi(2212)) {
                     h->Fill(e.get_leading_nofsi(2212).P(), e.get_weight());
                 }
             }},
            {"protonP_nofsi_qe",
             [](event &e, TH1D *h) {
                 if(e.get_mode() == event::channel::QE)
                 for (const auto &[_, p] : e.get_particle_nofsi(2212)) {
                     h->Fill(p.P(), e.get_weight());
                 }
             }},
            {"leadingP_nofsi_qe",
             [](event &e, TH1D *h) {
                 if (e.count_particle_nofsi(2212) && e.get_mode() == event::channel::QE) {
                     h->Fill(e.get_leading_nofsi(2212).P(), e.get_weight());
                 }
             }},
            {"protonP_nocut_qe",
             [](event &e, TH1D *h) {
                 if (e.get_mode() == event::channel::QE)
                 for (const auto &[_, p] : e.get_particle_out(2212)) {
                     h->Fill(p.P(), e.get_weight());
                 }
             }}
        };
        void add_event(event &e) {
            for (const auto &[name, func] : plot_funcs) {
                func(e, histos[name].get());
            }
        }
        plots() {
            for (const auto &[name, func] : plot_funcs) {
                histos.emplace(name, std::make_unique<TH1D>(name.c_str(), name.c_str(), bin_count, pmin, pmax));
                plot_names.insert(name);
            }
        }
        void save(std::filesystem::path file_path) {
            std::filesystem::create_directory(file_path.parent_path());
            TFile root_file{file_path.c_str(), "RECREATE"};
            for (const auto &[name, histo] : histos) {
                histo->SetDirectory(&root_file);
            }
            root_file.Write();
            for (const auto &[name, histo] : histos) {
                histo->SetDirectory(0);
            }
            root_file.Close();
        }
    } plot;

public:
    analysis(){};
    analysis(const analysis &) = delete;
    analysis(analysis &&) = default;
};
event::channel get_mode_genie(const TObjString &code) {
    if (code.GetString().Contains("QES")) {
        return event::channel::QE;
    } else if (code.GetString().Contains("RES")) {
        return event::channel::RES;
    } else if (code.GetString().Contains("DIS")) {
        return event::channel::DIS;
    } else if (code.GetString().Contains("MEC")) {
        return event::channel::MEC;
    }
    return event::channel::Other;
}
class run_manager_genie : public analysis {
public:
    const size_t event_count{};
    TGraph *sp;
    TH1D enu;
    run_manager_genie(size_t count, TGraph *spline) : event_count(count), sp(spline), enu("enu", "enu", bin_count, pmin, pmax) {}
    run_manager_genie(const run_manager_genie &) = delete;
    run_manager_genie(run_manager_genie &&) { plot = std::move(plot); }
    // std::mutex run_lock;
    class thread_object : public analysis {
    public:
        thread_object(const thread_object &) = default;
        thread_object(thread_object &&) = default;
        run_manager_genie &thisrun;
        TH1D enu;
        thread_object(run_manager_genie &run) : thisrun(run), enu("enu", "enu", bin_count, pmin, pmax) {}
        void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus, auto &&StdHepP4, auto &&EvtCode) {
            if (!EvtCode.GetString().Contains("Weak[CC]")) {
                return;
            } // focus on CC interactions
            event e;
            e.set_mode(get_mode_genie(EvtCode));
            for (int i = 0; i < StdHepN; ++i) {
                auto pdg = StdHepPdg[i];
                if (StdHepPdg[i] == 1000000010) {
                    pdg = 2112;
                }
                switch (StdHepStatus[i]) {
                case 0:
                    e.add_particle_in(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                case 1:
                    e.add_particle_out(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                default:
                    break;
                }
            }
            enu.Fill(e.get_enu());
            plot.add_event(e);
        }
        void finalize() {
            thisrun.enu.Add(&enu);
            for (const auto &[name, hist] : plot.histos) {
                thisrun.plot.histos[name]->Add(hist.get());
            }
        }
    };
    thread_object get_thread_object() { return thread_object(*this); }
    void finalize() {
        double xsec = get_xsec(&enu, sp) / 12. * 1e-38;
        std::cout << event_count << " events processed" << std::endl;
        std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;
        enu.Scale(xsec / event_count / bin_wid);
        for (const auto &[name, hist] : plot.histos) {
            hist->Scale(xsec / event_count / bin_wid);
        }
    }
    static run_manager_genie run_analysis(nlohmann::json &config) {
        constexpr size_t MAX_COUNT = 128;
        chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], TObjString> chain(
            config["input_file_list"], "gRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"}, std::thread::hardware_concurrency());
        auto spline_file = get_object<TGraph>(config["spline_file"], config["spline_path"]);

        return chain.run<run_manager_genie>(chain.get_entries(), spline_file.get());
    }
};
event::channel getmode_nuwro(TObjString &code) {
    const int mode = code.GetString().Atoi();
    // cout << "mode = " << mode << endl;
    switch (mode) {
    case 1:
        return event::channel::QE;
        break;
    case 11:
        return event::channel::RES;
        break;
    case 26:
        return event::channel::DIS;
        break;
    case 2:
        return event::channel::MEC;
        break;
    case 16:
    default:
        return event::channel::Other;
        break;
    }
}
class run_manager_nuwro : public analysis {
public:
    const size_t event_count{};
    double xsec{};
    run_manager_nuwro(size_t count) : event_count(count) {}
    run_manager_nuwro(const run_manager_nuwro &) = delete;
    run_manager_nuwro(run_manager_nuwro &&) { plot = std::move(plot); }
    // std::mutex run_lock;
    class thread_object : public analysis {
    public:
        thread_object(const thread_object &) = default;
        thread_object(thread_object &&) = default;
        run_manager_nuwro &thisrun;
        double xsec{};
        thread_object(run_manager_nuwro &run) : thisrun(run) {}
        void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus, auto &&StdHepP4, auto &&EvtWght, auto &&EvtCode) {
            event e;
            e.set_mode(getmode_nuwro(EvtCode));
            for (int i = 0; i < StdHepN; ++i) {
                auto pdg = StdHepPdg[i];
                if (StdHepPdg[i] == 1000000010) {
                    pdg = 2112;
                }
                switch (StdHepStatus[i]) {
                case 0:
                    e.add_particle_in(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                case 1:
                    e.add_particle_out(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                case 2:
                    e.add_particle_nofsi(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                    break;
                default:
                    break;
                }
            }
            xsec += EvtWght / thisrun.event_count * 1e-38;
            plot.add_event(e);
        }
        void finalize() {
            // thisrun.enu.Add(&enu);
            for (const auto &[name, hist] : plot.histos) {
                thisrun.plot.histos[name]->Add(hist.get());
            }
            thisrun.xsec += xsec;
        }
    };
    thread_object get_thread_object() { return thread_object(*this); }
    void finalize() {
        std::cout << event_count << " events processed" << std::endl;
        std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;

        // enu.Scale(xsec / event_count / bin_wid);
        for (const auto &[name, hist] : plot.histos) {
            hist->Scale(xsec / event_count / bin_wid);
        }
    }
    static run_manager_nuwro run_analysis(nlohmann::json &config) {
        constexpr size_t MAX_COUNT = 128;
        chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4], double, TObjString> chain(
            config["input_file_list"], "nRooTracker", {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght", "EvtCode"},
            std::thread::hardware_concurrency());

        return chain.run<run_manager_nuwro>(chain.get_entries());
    }
};