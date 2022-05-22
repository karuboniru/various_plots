#pragma once

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
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

class analysis {
public:
    constexpr static double pmin = 0, pmax = 20., bin_wid = 4e-3;
    constexpr static size_t bin_count = (pmax - pmin) / bin_wid;
    class plots {
    public:
        TH1D protonE, protonP, enu;
        TH1D protonE_nocut, protonP_nocut;
        TH1D leadingP, leadingP_nocut;
        void add_event(event &e) {
            enu.Fill(e.get_enu());
            auto phase = e.TKI_phase_cut();
            auto p_leading_p = e.get_leading_proton();
            if (e.count_particle_out(2212)) {
                leadingP_nocut.Fill(p_leading_p.P());
                if (phase) {
                    leadingP.Fill(p_leading_p.P());
                }
                for (const auto &[_, p] : e.get_particle_out(2212)) {
                    protonE_nocut.Fill(p.E() - p.M());
                    protonP_nocut.Fill(p.P());
                    if (phase) {
                        protonE.Fill(p.E() - p.M());
                        protonP.Fill(p.P());
                    }
                }
            }
        }
        plots()
            : protonE("protonE", "protonE", bin_count, pmin, pmax),                   //
              protonP("protonP", "protonP", bin_count, pmin, pmax),                   //
              enu("enu", "enu", bin_count, pmin, pmax),                               // aka. hccrate
              protonE_nocut("protonE_nocut", "protonE_nocut", bin_count, pmin, pmax), //
              protonP_nocut("protonP_nocut", "protonP_nocut", bin_count, pmin, pmax), //
              leadingP("leadingP", "leadingP", bin_count, pmin, pmax),                //
              leadingP_nocut("leadingP_nocut", "leadingP_nocut", bin_count, pmin, pmax) {}
    } plot;

public:
    analysis(){};
    analysis(const analysis &) = default;
    analysis(analysis &&) = default;
};

class run_manager : public analysis {
public:
    const size_t event_count{};
    TGraph *sp;
    run_manager(size_t count, TGraph *spline) : event_count(count), sp(spline) {}
    run_manager(const run_manager &) = delete;
    run_manager(run_manager &&) { plot = std::move(plot); }
    // std::mutex run_lock;
    class thread_object : public analysis {
    public:
        thread_object(const thread_object &) = default;
        thread_object(thread_object &&) = default;
        run_manager &thisrun;
        thread_object(run_manager &run) : thisrun(run) {}
        void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus, auto &&StdHepP4, auto &&EvtCode) {
            if (!EvtCode.GetString().Contains("Weak[CC]")) {
                return;
            }
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
            plot.add_event(e);
        }
        void finalize() {
            thisrun.plot.protonE.Add(&plot.protonE);
            thisrun.plot.protonP.Add(&plot.protonP);
            thisrun.plot.enu.Add(&plot.enu);
            thisrun.plot.protonE_nocut.Add(&plot.protonE_nocut);
            thisrun.plot.protonP_nocut.Add(&plot.protonP_nocut);
            thisrun.plot.leadingP.Add(&plot.leadingP);
            thisrun.plot.leadingP_nocut.Add(&plot.leadingP_nocut);
        }
    };
    thread_object get_thread_object() { return thread_object(*this); }
    void finalize() {
        std::cout << event_count << " events processed" << std::endl;
        double xsec = get_xsec(&plot.enu, sp);
        std::cout << "xsec = " << xsec << std::endl;
        plot.protonE.Scale(xsec / event_count / bin_wid);
        plot.protonP.Scale(xsec / event_count / bin_wid);
        plot.enu.Scale(xsec / event_count / bin_wid);
        plot.protonE_nocut.Scale(xsec / event_count / bin_wid);
        plot.protonP_nocut.Scale(xsec / event_count / bin_wid);
        plot.leadingP.Scale(xsec / event_count / bin_wid);
        plot.leadingP_nocut.Scale(xsec / event_count / bin_wid);
    }
};
