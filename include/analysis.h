#pragma once

#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TSpline.h>
#include <TStyle.h>
#include <chain_helper.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::string_literals::operator""s;

std::pair<double, double> get_xsec(TH1 *h_rate, TGraph *spline) {
    double fluxint{};
    // spline->SaveAs("wrong.root");
    TSpline3 sp("sp", spline);
    TF1 func(
        "spline", [&](double *x, double *) { return sp.Eval(*x); }, 0, h_rate->GetXaxis()->GetXmax(), 0);
    for (int ii = 1; ii <= h_rate->GetNbinsX(); ii++) {
        double bin_c = h_rate->GetBinContent(ii);
        double bin_up = h_rate->GetXaxis()->GetBinUpEdge(ii);
        double bin_low = h_rate->GetXaxis()->GetBinLowEdge(ii);
        double bin_width = bin_up - bin_low;
        if (bin_c < 1 || func.Integral(bin_low, bin_up) == 0) {
            continue;
        }
        // std::cout << "func.Integral(" << bin_low << ", " << bin_up << ")\t" << func.Integral(bin_low, bin_up) << std::endl;
        fluxint += bin_c / func.Integral(bin_low, bin_up) * bin_width;
    }
    // std::cout << "fluxint\t" << fluxint << std::endl;
    // std::cout << "h_rate->Integral()\t" << h_rate->Integral() << std::endl;
    return {h_rate->Integral(), h_rate->Integral() / fluxint};
}

template <typename T> std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
    TFile root_file{file_path.c_str(), "READ"};
    auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
    assert(objptr);
    return std::unique_ptr<T>{objptr};
}

class analysis {
public:
    constexpr static double pmin = 0, pmax = 40., bin_wid = 4e-3;
    constexpr static size_t bin_count = (pmax - pmin) / bin_wid;
    class plots {
    public:
        // TH1D enu;
        std::unordered_set<std::string> plot_names{};
        std::unordered_map<std::string, std::unique_ptr<TH1D>> histos;
        std::unordered_map<std::string, std::function<void(event &, TH1D *)>> plot_funcs{
            {"protonP",
             [](event &e, TH1D *h) {
                 for (const auto &[_, p] : e.get_particle_out(2212)) {
                     h->Fill(p.P(), e.get_weight());
                 }
             }}, // proton momentum (particle by particle)
            {"leadingP",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(2212)) {
                     h->Fill(e.get_leading_proton().P(), e.get_weight());
                 }
             }}, // leading proton momentum (event by event), no cut
            {"muonP",
             [](event &e, TH1D *h) {
                 if (e.count_particle_out(13)) {
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
                     if (p.P() != 0)
                         h->Fill(p.P(), e.get_weight());
                 }
             }},
            {"leadingP_nofsi",
             [](event &e, TH1D *h) {
                 if (e.count_particle_nofsi(2212)) {
                     h->Fill(e.get_leading_nofsi(2212).P(), e.get_weight());
                 }
             }},
            {"enu", [](event &e, TH1D *h) { h->Fill(e.get_enu(), e.get_weight()); }},
            {"kaonE",
             [](event &e, TH1D *h) {
                 for (const auto &[_, p] : e.get_particle_out(321)) {
                     h->Fill(p.E() - p.M(), e.get_weight());
                 }
             }},
            {"leadingPi0", [](event &e, TH1D *h) {
                 if (e.count_particle_out(111)) {
                     h->Fill(e.get_leading_out(111).P(), e.get_weight());
                 }
             }}};
        void add_event(event &e) {
            for (const auto &[name, func] : plot_funcs) {
                func(e, histos[name].get());
            }
        }
        void add_cut(std::string cutname, std::function<bool(event &)> func) {
            decltype(plot_funcs) cut_plot_funcs{};
            for (const auto &[name, h_func] : plot_funcs) {
                auto new_name = name + "_" + cutname;
                histos.emplace(new_name, std::make_unique<TH1D>(new_name.c_str(), new_name.c_str(), bin_count, pmin, pmax));
                plot_names.insert(new_name);
                cut_plot_funcs.emplace(new_name, [&, func](event &e, TH1D *h) {
                    if (func(e)) {
                        h_func(e, h);
                    }
                });
            }
            plot_funcs.merge(cut_plot_funcs);
        }
        void add_cut_group(std::initializer_list<std::pair<std::string, std::function<bool(event &)>>> cut_group) {
            decltype(plot_funcs) cut_plot_funcs{};
            for (const auto &[name, h_func] : plot_funcs) {
                for (auto &cut : cut_group) {
                    auto [cutname, cut_func] = cut;
                    auto new_name = name + "_" + cutname;
                    histos.emplace(new_name, std::make_unique<TH1D>(new_name.c_str(), new_name.c_str(), bin_count, pmin, pmax));
                    plot_names.insert(new_name);
                    cut_plot_funcs.emplace(new_name, [&, cut_func](event &e, TH1D *h) {
                        if (cut_func(e)) {
                            h_func(e, h);
                        }
                    });
                }
            }
            plot_funcs.merge(cut_plot_funcs);
        }
        plots() {
            for (const auto &[name, func] : plot_funcs) {
                histos.emplace(name, std::make_unique<TH1D>(name.c_str(), name.c_str(), bin_count, pmin, pmax));
                plot_names.insert(name);
            }
            // add_cut("TKI", [](event &e) { return e.TKI_mu_cut(); });
            // add_cut("TKIp", [](event &e) { return e.TKI_mu_cut(); });
            add_cut("mucut", [](event &e) { return e.get_leading_out(14).P() > 1.0; });
            add_cut("epi", [](event &e) {
                size_t count_e{}, count_pi{};
                for (const auto &[id, p] : e.get_particle_out()) {
                    if (id == -11) {
                        count_e++;
                    } else if (id == 111) {
                        count_pi++;
                    } else {
                        return false;
                    }
                }
                return (count_e == 1) && (count_pi == 1);
            }); // for PDK final state
            add_cut_group({{"qe", [](event &e) { return e.get_mode() == event::channel::QE; }},
                           {"res", [](event &e) { return e.get_mode() == event::channel::RES; }},
                           {"dis", [](event &e) { return e.get_mode() == event::channel::DIS; }},
                           {"2p2h", [](event &e) { return e.get_mode() == event::channel::MEC; }}});
            add_cut_group({{"1pi", [](event &e) { return e.count_particle_out(211) + e.count_particle_out(111) + e.count_particle_out(-211) == 1; }},
                           {"2pi", [](event &e) { return e.count_particle_out(211) + e.count_particle_out(111) + e.count_particle_out(-211) == 2; }},
                           {"0pi", [](event &e) { return e.count_particle_out(211) + e.count_particle_out(111) + e.count_particle_out(-211) == 0; }},
                           {"Mpi", [](event &e) { return e.count_particle_out(211) + e.count_particle_out(111) + e.count_particle_out(-211) > 2; }},
                           {"1pi0", [](event &e) { return e.count_particle_out(111) == 1; }},
                           {"2pi0", [](event &e) { return e.count_particle_out(111) == 2; }},
                           {"0pi0", [](event &e) { return e.count_particle_out(111) == 0; }},
                           {"Mpi0", [](event &e) { return e.count_particle_out(111) > 2; }},
                           {"1pip", [](event &e) { return e.count_particle_out(211) == 1; }},
                           {"2pip", [](event &e) { return e.count_particle_out(211) == 2; }},
                           {"0pip", [](event &e) { return e.count_particle_out(211) == 0; }},
                           {"Mpip", [](event &e) { return e.count_particle_out(211) > 2; }},
                           {"1pim", [](event &e) { return e.count_particle_out(-211) == 1; }},
                           {"2pim", [](event &e) { return e.count_particle_out(-211) == 2; }},
                           {"0pim", [](event &e) { return e.count_particle_out(-211) == 0; }},
                           {"Mpim", [](event &e) { return e.count_particle_out(-211) > 2; }}});
        }
        void save(std::filesystem::path file_path) {
            // std::filesystem::create_directory(file_path.parent_path());
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
    run_manager_genie(run_manager_genie &&lhs) = default;
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
            if (StdHepPdg[0] != 14)
                return;
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
        auto [tot, xsec] = get_xsec(&enu, sp);
        xsec *= 1 / 12. * 1e-38;
        std::cout << event_count << " events processed" << std::endl;
        std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;
        enu.Scale(xsec / tot / bin_wid);
        for (const auto &[name, hist] : plot.histos) {
            hist->Scale(xsec / tot / bin_wid);
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
    run_manager_nuwro(run_manager_nuwro &&lhs) = default;
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
            if (getmode_nuwro(EvtCode) == event::channel::Other) {
                return;
            }
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

event::channel getmode_gibuu(int modeid) {
    if (modeid == 1) {
        return event::channel::QE;
    } else if (modeid > 1 && modeid <= 33) {
        return event::channel::RES;
    } else if (modeid == 34) {
        return event::channel::DIS;
    } else if (modeid == 35) {
        return event::channel::MEC;
    }
    return event::channel::Other;
}

class run_manager_gibuu : public analysis {
public:
    const size_t nrun{};
    const int nupdg, nucpdg;
    double xsec{};
    run_manager_gibuu(size_t nrun, int nupdg, int nucpdg) : nrun(nrun), nupdg(nupdg), nucpdg(nucpdg) {}
    run_manager_gibuu(const run_manager_gibuu &) = delete;
    run_manager_gibuu(run_manager_gibuu &&lhs) = default;
    // std::mutex run_lock;
    class thread_object : public analysis {
    public:
        thread_object(const thread_object &) = default;
        thread_object(thread_object &&) = default;
        run_manager_gibuu &thisrun;
        thread_object(run_manager_gibuu &run) : thisrun(run) {}
        void run(auto &&weight, auto &&barcode, auto &&E, auto &&Px, auto &&Py, auto &&Pz, auto &&evType, auto &&lepIn_E, auto &&lepIn_Px, auto &&lepIn_Py,
                 auto &&lepIn_Pz, auto &&lepOut_E, auto &&lepOut_Px, auto &&lepOut_Py, auto &&lepOut_Pz, auto &&nuc_E, auto &&nuc_Px, auto &&nuc_Py,
                 auto &&nuc_Pz) {
            event e;
            e.set_weight(weight);
            e.set_mode(getmode_gibuu(evType));
            const auto size_o = E.size();
            for (size_t i = 0; i < size_o; ++i) {
                e.add_particle_out(barcode[i], TLorentzVector(Px[i], Py[i], Pz[i], E[i]));
            }
            e.add_particle_out(thisrun.nupdg > 0 ? thisrun.nupdg - 1 : thisrun.nupdg + 1, TLorentzVector(lepOut_Px, lepOut_Py, lepOut_Pz, lepOut_E));
            e.add_particle_in(thisrun.nupdg, TLorentzVector(lepIn_Px, lepIn_Py, lepIn_Pz, lepIn_E));
            e.add_particle_in(thisrun.nucpdg, TLorentzVector(nuc_Px, nuc_Py, nuc_Pz, nuc_E));
            plot.add_event(e);
        }
        void finalize() {
            // thisrun.enu.Add(&enu);
            for (const auto &[name, hist] : plot.histos) {
                thisrun.plot.histos[name]->Add(hist.get());
            }
        }
    };
    thread_object get_thread_object() { return thread_object(*this); }
    void finalize() {
        // std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;

        // enu.Scale(xsec / event_count / bin_wid);
        for (const auto &[name, hist] : plot.histos) {
            hist->Scale(1. / nrun * 1e-38 / bin_wid);
        }
    }
    static run_manager_gibuu run_analysis(nlohmann::json &config) {
        // constexpr size_t MAX_COUNT = 128;
        chain_runner<double,                                                                             // weight
                     std::vector<int>,                                                                   // barcode. seems to be pdgcode
                     std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, // E, Px, Py, Pz
                     int,                                                                                // evType
                     double, double, double, double,                                                     // lepIn_E, lepIn_px, lepIn_py, lepIn_pz
                     double, double, double, double,                                                     // lepOut_E, lepOut_px, lepOut_py, lepOut_pz
                     double, double, double, double                                                      // nuc_E, nuc_px, nuc_py, nuc_pz
                     >
            chain(config["input_file_list"], "RootTuple",
                  {"weight", "barcode", "E", "Px", "Py", "Pz", "evType", "lepIn_E", "lepIn_Px", "lepIn_Py", "lepIn_Pz", "lepOut_E", "lepOut_Px", "lepOut_Py",
                   "lepOut_Pz", "nuc_E", "nuc_Px", "nuc_Py", "nuc_Pz"},
                  std::thread::hardware_concurrency());

        return chain.run<run_manager_gibuu>(config.value("nrun", config["input_file_list"].size()), config.value("nupdg", 14),
                                            config.value("nucpdg", 1000060120));
    }
};