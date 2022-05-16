#pragma once

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <algorithm>
#include <chain_helper.h>
#include <deque>
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

using std::string_literals::operator""s;

class analysis {
public:
    double xsec{};
    constexpr static double pmin = 0, pmax = 20., bin_wid = 4e-3;
    constexpr static size_t bin_count = (pmax - pmin) / bin_wid;
    class plots {
    public:
        TH1D protonE, protonP, enu;
        TH1D protonE_nocut, protonP_nocut;
        void add_event(event &e) {
            enu.Fill(e.get_enu());
            auto phase = e.TKI_phase_cut();
            for (auto &p : e.get_particle(2212)) {
                protonE_nocut.Fill(p.E() - p.M());
                protonP_nocut.Fill(p.P());
                if (phase) {
                    protonE.Fill(p.E() - p.M());
                    protonP.Fill(p.P());
                }
            }
        }
        plots()
            : protonE("protonE", "protonE", bin_count, pmin, pmax),                   //
              protonP("protonP", "protonP", bin_count, pmin, pmax),                   //
              enu("enu", "enu", bin_count, pmin, pmax),                               //
              protonE_nocut("protonE_nocut", "protonE_nocut", bin_count, pmin, pmax), //
              protonP_nocut("protonP_nocut", "protonP_nocut", bin_count, pmin, pmax)  //
        {}
    } p_all;

public:
    analysis(){};
    analysis(const analysis &) = default;
    analysis(analysis &&) = default;
};

class run_manager : public analysis {
public:
    const size_t event_count{};
    run_manager(size_t count) : event_count(count) {}
    run_manager(const run_manager &) = delete;
    run_manager(run_manager &&) { p_all = std::move(p_all); }
    // std::mutex run_lock;
    class thread_object : public analysis {
    public:
        thread_object(const thread_object &) = default;
        thread_object(thread_object &&) = default;
        run_manager &thisrun;
        double xsec{};
        thread_object(run_manager &run) : thisrun(run) {}
        void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus, auto &&StdHepP4, auto &&EvtWght) {
            event e;
            for (int i = 0; i < StdHepN; ++i) {
                auto pdg = StdHepPdg[i];
                if (StdHepPdg[i] == 1000000010) {
                    pdg = 2112;
                }
                if (StdHepStatus[i] == 0) {
                    e.add_particle_in(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                } else if (StdHepStatus[i] == 1) {
                    e.add_particle_out(pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2], StdHepP4[i][3]));
                }
            }
            xsec += EvtWght / thisrun.event_count * 1e-38;
            p_all.add_event(e);
        }
        void finalize() {
            // p_all.finalize_plot(xsec, thisrun.event_count);
            {
                std::cout << "xsec per thread" << xsec << std::endl;
                // std::lock_guard<std::mutex> l(thisrun.run_lock);
                thisrun.p_all.protonE.Add(&p_all.protonE);
                thisrun.p_all.protonP.Add(&p_all.protonP);
                thisrun.p_all.enu.Add(&p_all.enu);
                thisrun.p_all.protonE_nocut.Add(&p_all.protonE_nocut);
                thisrun.p_all.protonP_nocut.Add(&p_all.protonP_nocut);
                thisrun.xsec += xsec;
            }
        }
    };
    thread_object get_thread_object() {
        // std::lock_guard lock(run_lock);
        return thread_object(*this);
    }
    void finalize() {
        std::cout << "xsec: " << xsec << std::endl;
        p_all.protonE.Scale(xsec / event_count / bin_wid);
        p_all.protonP.Scale(xsec / event_count / bin_wid);
        p_all.enu.Scale(xsec / event_count / bin_wid);
        p_all.protonE_nocut.Scale(xsec / event_count / bin_wid);
        p_all.protonP_nocut.Scale(xsec / event_count / bin_wid);
    }
};
