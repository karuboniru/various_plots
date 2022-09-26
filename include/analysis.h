#pragma once

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TObjString.h>
#include <TSpline.h>
#include <TStyle.h>
#include <chain_helper.h>
#include <event.h>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

using std::string_literals::operator""s;

std::pair<double, double> get_xsec(TH1 *h_rate, TGraph *spline) {
  double fluxint{};
  // spline->SaveAs("wrong.root");
  TSpline3 sp("sp", spline);
  TF1 func(
      "spline", [&](double *x, double *) { return sp.Eval(*x); }, 0,
      h_rate->GetXaxis()->GetXmax(), 0);
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
  double event_rate = h_rate->Integral();
  return {event_rate, event_rate / fluxint};
}

template <typename T>
std::unique_ptr<T> get_object(std::string file_path, std::string obj_path) {
  TFile root_file{file_path.c_str(), "READ"};
  auto objptr = static_cast<T *>(root_file.Get(obj_path.c_str())->Clone());
  assert(objptr);
  return std::unique_ptr<T>{objptr};
}

class analysis {
public:
  constexpr static double pmin = 0, pmax = 25, bin_wid = 1e-3;
  constexpr static int m_2d_factor = 1;
  constexpr static size_t bin_count = (pmax - pmin) / bin_wid;
  class plots {
  public:
    // TH1D enu;
    // using analysis_funcs_t = std::pair<std::string, std::function<void(event
    // &, TH1 *)>>;
    using variables_t = std::vector<double>;
    using cut_t = std::pair<std::string, std::function<bool(event &)>>;
    using axix_t = std::tuple<double, double, int>;

    using analysis_funcs_t =
        std::tuple<std::string, std::function<variables_t(event &)>, axix_t>;

    using axis_pair_t = std::pair<std::string, std::string>;

    using cut_group_t = std::list<cut_t>;
    std::list<analysis_funcs_t> m_analysis_funcs;
    std::list<cut_group_t> m_cut_groups;
    std::list<cut_t> m_cuts;
    std::list<axis_pair_t> hist_2d_pairs{};
    std::unordered_set<std::string> plot_names{};
    std::unordered_map<std::string, std::unique_ptr<TH1>> histos;
    std::unordered_map<std::string, std::function<variables_t(event &)>>
        plot_funcs{};
    std::unordered_map<std::string, axix_t> plot_axes{};
    std::vector<cut_t> all_cuts{};
    // std::vector<std::pair<std::string, std::string>> hist_2d_pairs{};
    template <class... Ts> struct overloaded : Ts... {
      using Ts::operator()...;
    };
    template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;
    void add_to_hist(std::string name, double var, event &e) {
      histos[name]->Fill(var);
      for (auto &[cutname, cutfunc] : all_cuts) {
        if (cutfunc(e)) {
          histos[name + "_" + cutname]->Fill(var, e.get_weight());
        }
      }
    }
    void add_event(event &e) {
      // for (const auto &[name, func] : plot_funcs) {
      //     func(e, histos[name].get());
      // }
      std::unordered_map<std::string, variables_t> vars{};
      for (const auto &[name, func] : plot_funcs) {
        auto &&var = func(e);
        for (auto &thisvar : var) {
          histos[name]->Fill(thisvar);
        }
        vars[name] = std::move(var);
      }
      vars.erase("");
      for (const auto &[name1, name2] : hist_2d_pairs) {
        auto name = "2d_" + name1 + "_" + name2;
        for (const auto &var1 : vars[name1]) {
          for (const auto &var2 : vars[name2]) {
            dynamic_cast<TH2 *>(histos[name].get())
                ->Fill(var1, var2, e.get_weight());
          }
        }
      }
      for (auto &[cutname, cutfunc] : all_cuts) {
        if (cutfunc(e)) {
          for (const auto &[name, vars] : vars) {
            for (const auto &var : vars) {
              histos[name + "_" + cutname]->Fill(var, e.get_weight());
            }
          }
          for (const auto &[name1, name2] : hist_2d_pairs) {
            auto name = "2d_" + name1 + "_" + name2 + "_" + cutname;
            for (const auto &var1 : vars[name1]) {
              for (const auto &var2 : vars[name2]) {
                dynamic_cast<TH2 *>(histos[name].get())
                    ->Fill(var1, var2, e.get_weight());
              }
            }
          }
        }
      }
    }
    void add_plot2d(std::string name, std::tuple<double, double, int> x,
                    std::tuple<double, double, int> y) {
      auto &&[downx, upx, binx] = x;
      auto &&[downy, upy, biny] = y;
      // std::cout << "Adding plot " << name << std::endl;
      plot_names.insert(name);
      histos.emplace(name,
                     std::make_unique<TH2D>(name.c_str(), name.c_str(), binx,
                                            downx, upx, biny, downy, upy));
    }
    void add_plot(std::string name, std::tuple<double, double, int> x) {
      auto &&[down, up, bin] = x;
      plot_names.insert(name);
      // std::cout << "Adding plot " << name << std::endl;
      histos.emplace(name, std::make_unique<TH1D>(name.c_str(), name.c_str(),
                                                  bin, down, up));
    }
    void add_plot(std::string name, std::string basename) {
      plot_names.insert(name);
      histos.emplace(
          name, dynamic_cast<TH1 *>(histos[basename]->Clone(name.c_str())));
      // if (basename.empty())
      //   basename = name;
      // if (name.substr(0, 2) == "2d") {
      //   // std::cerr << "Adding 2d cut plot " << name << std::endl;
      //   histos.emplace(name,
      //                  std::make_unique<TH2D>(
      //                      name.c_str(), basename.c_str(),
      //                      bin_count / m_2d_factor, pmin, pmax / m_2d_factor,
      //                      bin_count / m_2d_factor, pmin, pmax /
      //                      m_2d_factor));
      // } else {
      //   // std::cerr << "Adding cut plot " << name << std::endl;
      //   histos.emplace(name,
      //                  std::make_unique<TH1D>(name.c_str(), basename.c_str(),
      //                                         bin_count, pmin, pmax));
      // }
    }
    void add_cut(std::string cutname, std::function<bool(event &)> func) {
      decltype(all_cuts) cut_plot_funcs{};
      std::vector<std::string> cut_plot_names{};
      for (const auto &[name, h_func] : all_cuts) {
        auto new_name = name + "_" + cutname;
        cut_plot_names.push_back(new_name);
        cut_plot_funcs.push_back(
            cut_t{new_name, [=](event &e) { return h_func(e) && func(e); }});
      }
      cut_plot_funcs.push_back(std::make_pair(cutname, func));
      cut_plot_names.push_back(cutname);
      //   add_cut_plots(cut_plot_names);
      all_cuts.insert(all_cuts.end(), cut_plot_funcs.begin(),
                      cut_plot_funcs.end());
    }
    // void print_list(auto &&list) {
    //   for (auto &&entry : list) {
    //     std::cout << entry.first << std::endl;
    //   }
    // }
    void add_cut_group(
        const std::list<std::pair<std::string, std::function<bool(event &)>>>
            cut_group) {
      decltype(all_cuts) cut_plot_funcs{};
      std::vector<std::string> cut_plot_names{};
      // std::cout << "cut_group: ";
      // print_list(cut_group);
      // std::cout << "cuts" << std::endl;
      // print_list(all_cuts);
      for (auto &[cutname, cut_func] : cut_group) {
        for (const auto &[name, h_func] : all_cuts) {
          auto new_name = name + "_" + cutname;
          cut_plot_names.push_back(new_name);
          cut_plot_funcs.push_back(cut_t{
              new_name, [=](event &e) { return h_func(e) && cut_func(e); }});
        }
        cut_plot_funcs.push_back(cut_t{cutname, cut_func});
        cut_plot_names.push_back(cutname);
      }
      //   add_cut_plots(cut_plot_names);
      all_cuts.insert(all_cuts.end(), cut_plot_funcs.begin(),
                      cut_plot_funcs.end());
      // std::cout << "finish one" << std::endl;
    }

    plots(std::list<analysis_funcs_t> analysis_funcs,
          std::list<cut_group_t> cut_groups, std::list<cut_t> cuts,
          std::list<axis_pair_t> hist_2d_pairs_i) {
      m_analysis_funcs = analysis_funcs;
      m_cut_groups = cut_groups;
      m_cuts = cuts;
      hist_2d_pairs = hist_2d_pairs_i;

      for (const auto &[name, func, axis] : m_analysis_funcs) {
        // std::cerr << "registered" << name << std::endl;
        plot_funcs[name] = func;
        plot_axes[name] = axis;
        add_plot(name, axis);
      }
      for (const auto &[name1, name2] : hist_2d_pairs) {
        auto &&axis1 = plot_axes[name1];
        auto &&axis2 = plot_axes[name2];
        add_plot2d("2d_" + name1 + "_" + name2, axis1, axis2);
      }
      // std::cerr << "Cut groups: " << std::endl;
      for (const auto &cut_group : m_cut_groups) {
        add_cut_group(cut_group);
      }
      // std::cerr << "Cuts: " << std::endl;
      for (const auto &[name, func] : m_cuts) {
        add_cut(name, func);
      }
      auto baseset = plot_names;
      for (auto &&name : baseset) {
        if (name.empty())
          continue;
        for (auto &&[cut, func] : all_cuts) {
          // std::terminate();
          // std::cout << "adding " << name << " " << cut << std::endl;
          plot_names.insert(name + "_" + cut);
          add_plot(name + "_" + cut, name);
        }
      }
      // std::cerr << "finished " << std::endl;
    }
    void save(std::filesystem::path file_path) {
      // std::filesystem::create_directory(file_path.parent_path());
      TFile root_file{file_path.c_str(), "RECREATE"};
      for (const auto &[name, histo] : histos) {
        if (!name.empty())
          histo->SetDirectory(&root_file);
      }
      root_file.Write();
      for (const auto &[name, histo] : histos) {
        if (!name.empty())
          histo->SetDirectory(0);
      }
      root_file.Close();
    }
  } plot;

public:
  // analysis(){};
  template <typename... Args>
  analysis(Args &&...args) : plot(std::forward<decltype(args)>(args)...) {}
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
  // run_manager_genie(size_t count, TGraph *spline) : event_count(count),
  // sp(spline), enu("enu", "enu", bin_count, pmin, pmax) {}
  template <typename... Args>
  run_manager_genie(size_t count, TGraph *spline, Args &&...args)
      : analysis(std::forward<Args>(args)...), event_count(count), sp(spline),
        enu("enu", "enu", bin_count, pmin, pmax) {}
  run_manager_genie(const run_manager_genie &) = delete;
  run_manager_genie(run_manager_genie &&lhs) = default;
  // std::mutex run_lock;
  class thread_object : public analysis {
  public:
    thread_object(const thread_object &) = delete;
    thread_object(thread_object &&) = default;
    run_manager_genie &thisrun;
    TH1D enu;
    template <typename... Args>
    thread_object(run_manager_genie &run, Args &&...args)
        : analysis(std::forward<Args>(args)...), thisrun(run),
          enu("enu", "enu", bin_count, pmin, pmax) {}
    void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus,
             auto &&StdHepP4, auto &&EvtCode) {
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
          e.add_particle_in(pdg,
                            TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                           StdHepP4[i][2], StdHepP4[i][3]));
          break;
        case 1:
          e.add_particle_out(pdg,
                             TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
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
  thread_object get_thread_object() {
    return thread_object(*this, plot.m_analysis_funcs, plot.m_cut_groups,
                         plot.m_cuts, plot.hist_2d_pairs);
  }
  void finalize() {
    auto [tot, xsec] = get_xsec(&enu, sp);
    xsec *= 1 / 12. * 1e-38;
    // auto cc_event_count = enu.Integral("");
    std::cout << event_count << " events processed" << std::endl;
    std::cout << tot / event_count << " CC ratio" << std::endl;
    std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;
    enu.Scale(xsec / tot / bin_wid);
    for (const auto &[name, hist] : plot.histos) {
      if (name.substr(0, 2) == "2d") {
        auto bin_wid_x =
            dynamic_cast<TH2 *>(hist.get())->GetXaxis()->GetBinWidth(1);
        auto bin_wid_y =
            dynamic_cast<TH2 *>(hist.get())->GetYaxis()->GetBinWidth(1);
        hist->Scale(xsec / tot / bin_wid_x / bin_wid_y);
      } else {
        auto bin_widx = hist->GetXaxis()->GetBinWidth(1);
        hist->Scale(xsec / tot / bin_widx);
      }
    }
  }
  template <typename... Args>
  static run_manager_genie run_analysis(nlohmann::json &config,
                                        Args &&...args) {
    constexpr size_t MAX_COUNT = 128;
    chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4],
                 TObjString>
        chain(config["input_file_list"], "gRooTracker",
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"},
              std::thread::hardware_concurrency());
    auto spline_file =
        get_object<TGraph>(config["spline_file"], config["spline_path"]);

    return chain.run<run_manager_genie>(chain.get_entries(), spline_file.get(),
                                        std::forward<Args>(args)...);
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
  const size_t event_count;
  double xsec{};
  // run_manager_nuwro(size_t count) : event_count(count) {}
  template <typename... Args>
  run_manager_nuwro(size_t count, Args &&...args)
      : analysis(std::forward<Args>(args)...), event_count(count) {}
  run_manager_nuwro(const run_manager_nuwro &) = delete;
  run_manager_nuwro(run_manager_nuwro &&lhs) = default;
  // std::mutex run_lock;
  class thread_object : public analysis {
  public:
    thread_object(const thread_object &) = delete;
    thread_object(thread_object &&) = default;
    run_manager_nuwro &thisrun;
    double xsec{};
    template <typename... Args>
    thread_object(run_manager_nuwro &run, Args &&...args)
        : analysis(std::forward<Args>(args)...), thisrun(run) {}
    void run(auto &&StdHepN, auto &&StdHepPdg, auto &&StdHepStatus,
             auto &&StdHepP4, auto &&EvtWght, auto &&EvtCode, auto &&nod) {
      event e;
      e.set_nod(nod);
      if (getmode_nuwro(EvtCode) == event::channel::Other) {
        return;
      }
      e.set_mode(getmode_nuwro(EvtCode));
      size_t proton_count{};
      for (int i = 0; i < StdHepN; ++i) {
        auto pdg = StdHepPdg[i];
        if (StdHepPdg[i] == 1000000010) {
          pdg = 2112;
        }
        switch (StdHepStatus[i]) {
        case 0:
          e.add_particle_in(pdg,
                            TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                           StdHepP4[i][2], StdHepP4[i][3]));
          break;
        case 1:
          e.add_particle_out(pdg,
                             TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                            StdHepP4[i][2], StdHepP4[i][3]));
          break;
        case 2: {
          TLorentzVector p4(StdHepP4[i][0], StdHepP4[i][1], StdHepP4[i][2],
                            StdHepP4[i][3]);
          e.add_particle_nofsi(pdg, p4);
          if (pdg == 2212 || pdg == 2112) {
            if (proton_count == 0) {
              e.setprimaryP(p4);
            }
            if (proton_count == 1) {
              e.setspectatorP(p4);
            }
            proton_count++;
          }
        } break;
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
  thread_object get_thread_object() {
    return thread_object(*this, plot.m_analysis_funcs, plot.m_cut_groups,
                         plot.m_cuts, plot.hist_2d_pairs);
  }
  void finalize() {
    std::cout << event_count << " events processed" << std::endl;
    std::cout << "averaged overall xsec = " << xsec << " [cm^2]" << std::endl;

    // enu.Scale(xsec / event_count / bin_wid);
    for (const auto &[name, hist] : plot.histos) {
      if (name.substr(0, 2) == "2d") {
        auto bin_wid_x =
            dynamic_cast<TH2 *>(hist.get())->GetXaxis()->GetBinWidth(1);
        auto bin_wid_y =
            dynamic_cast<TH2 *>(hist.get())->GetYaxis()->GetBinWidth(1);
        hist->Scale(xsec / event_count / bin_wid_x / bin_wid_y);
      } else {
        auto bin_widx = hist->GetXaxis()->GetBinWidth(1);
        hist->Scale(xsec / event_count / bin_widx);
      }
    }
  }
  template <typename... Args>
  static run_manager_nuwro run_analysis(nlohmann::json &config,
                                        Args &&...args) {
    constexpr size_t MAX_COUNT = 128;
    chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4],
                 double, TObjString, int[18]>
        chain(config["input_file_list"], "nRooTracker",
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtWght",
               "EvtCode", "nod"},
              std::thread::hardware_concurrency());

    return chain.run<run_manager_nuwro>(chain.get_entries(),
                                        std::forward<Args>(args)...);
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
  // run_manager_gibuu(size_t nrun, int nupdg, int nucpdg) : nrun(nrun),
  // nupdg(nupdg), nucpdg(nucpdg) {}
  template <typename... Args>
  run_manager_gibuu(size_t nrun, int nupdg, int nucpdg, Args &&...args)
      : analysis(std::forward<Args>(args)...), nrun(nrun), nupdg(nupdg),
        nucpdg(nucpdg) {}
  run_manager_gibuu(const run_manager_gibuu &) = delete;
  run_manager_gibuu(run_manager_gibuu &&lhs) = default;
  // std::mutex run_lock;
  class thread_object : public analysis {
  public:
    thread_object(const thread_object &) = delete;
    thread_object(thread_object &&) = default;
    run_manager_gibuu &thisrun;
    template <typename... Args>
    thread_object(run_manager_gibuu &run, Args &&...args)
        : analysis(std::forward<Args>(args)...), thisrun(run) {}
    void run(auto &&weight, auto &&barcode, auto &&E, auto &&Px, auto &&Py,
             auto &&Pz, auto &&evType, auto &&lepIn_E, auto &&lepIn_Px,
             auto &&lepIn_Py, auto &&lepIn_Pz, auto &&lepOut_E,
             auto &&lepOut_Px, auto &&lepOut_Py, auto &&lepOut_Pz, auto &&nuc_E,
             auto &&nuc_Px, auto &&nuc_Py, auto &&nuc_Pz) {
      event e;
      e.set_weight(weight);
      e.set_mode(getmode_gibuu(evType));
      const auto size_o = E.size();
      for (size_t i = 0; i < size_o; ++i) {
        e.add_particle_out(barcode[i],
                           TLorentzVector(Px[i], Py[i], Pz[i], E[i]));
      }
      e.add_particle_out(
          thisrun.nupdg > 0 ? thisrun.nupdg - 1 : thisrun.nupdg + 1,
          TLorentzVector(lepOut_Px, lepOut_Py, lepOut_Pz, lepOut_E));
      e.add_particle_in(thisrun.nupdg,
                        TLorentzVector(lepIn_Px, lepIn_Py, lepIn_Pz, lepIn_E));
      e.add_particle_in(thisrun.nucpdg,
                        TLorentzVector(nuc_Px, nuc_Py, nuc_Pz, nuc_E));
      plot.add_event(e);
    }
    void finalize() {
      // thisrun.enu.Add(&enu);
      for (const auto &[name, hist] : plot.histos) {
        thisrun.plot.histos[name]->Add(hist.get());
      }
    }
  };
  thread_object get_thread_object() {
    return thread_object(*this, plot.m_analysis_funcs, plot.m_cut_groups,
                         plot.m_cuts, plot.hist_2d_pairs);
  }
  void finalize() {
    // std::cout << "averaged overall xsec = " << xsec << " [cm^2]" <<
    // std::endl;

    // enu.Scale(xsec / event_count / bin_wid);
    for (const auto &[name, hist] : plot.histos) {
      if (name.substr(0, 2) == "2d")
        hist->Scale(1. / nrun * 1e-38 / bin_wid / bin_wid);
      else
        hist->Scale(1. / nrun * 1e-38 / bin_wid);
    }
  }
  template <typename... Args>
  static run_manager_gibuu run_analysis(nlohmann::json &config,
                                        Args &&...args) {
    // constexpr size_t MAX_COUNT = 128;
    chain_runner<double,           // weight
                 std::vector<int>, // barcode. seems to be pdgcode
                 std::vector<double>, std::vector<double>, std::vector<double>,
                 std::vector<double>, // E, Px, Py, Pz
                 int,                 // evType
                 double, double, double,
                 double, // lepIn_E, lepIn_px, lepIn_py, lepIn_pz
                 double, double, double,
                 double, // lepOut_E, lepOut_px, lepOut_py, lepOut_pz
                 double, double, double, double // nuc_E, nuc_px, nuc_py, nuc_pz
                 >
        chain(config["input_file_list"], "RootTuple",
              {"weight", "barcode", "E", "Px", "Py", "Pz", "evType", "lepIn_E",
               "lepIn_Px", "lepIn_Py", "lepIn_Pz", "lepOut_E", "lepOut_Px",
               "lepOut_Py", "lepOut_Pz", "nuc_E", "nuc_Px", "nuc_Py", "nuc_Pz"},
              std::thread::hardware_concurrency());

    return chain.run<run_manager_gibuu>(
        config.value("nrun", config["input_file_list"].size()),
        config.value("nupdg", 14), config.value("nucpdg", 1000060120),
        std::forward<Args>(args)...);
  }
};

class run_manager_gibuu_pdk : public analysis {
public:
  // run_manager_gibuu_pdk(size_t nrun, int nupdg, int nucpdg) : nrun(nrun),
  // nupdg(nupdg), nucpdg(nucpdg) {}
  size_t count{};
  template <typename... Args>
  run_manager_gibuu_pdk(Args &&...args)
      : analysis(std::forward<Args>(args)...) {}
  run_manager_gibuu_pdk(const run_manager_gibuu_pdk &) = delete;
  run_manager_gibuu_pdk(run_manager_gibuu_pdk &&lhs) = default;
  // std::mutex run_lock;
  class thread_object : public analysis {
    size_t count{};

  public:
    thread_object(const thread_object &) = delete;
    thread_object(thread_object &&) = default;
    run_manager_gibuu_pdk &thisrun;
    template <typename... Args>
    thread_object(run_manager_gibuu_pdk &run, Args &&...args)
        : analysis(std::forward<Args>(args)...), thisrun(run) {}
    void run(auto &&nparticles, auto &&pdg, auto &&status, auto P) {
      event e;
      for (int i = 0; i < nparticles; ++i) {
        switch (status[i]) {
        case 0:
          e.add_particle_in(pdg[i],
                            TLorentzVector(P[i][1], P[i][2], P[i][3], P[i][0]));
          break;
        case 1:
          e.add_particle_out(
              pdg[i], TLorentzVector(P[i][1], P[i][2], P[i][3], P[i][0]));
          break;
        case 2:
          e.add_particle_nofsi(
              pdg[i], TLorentzVector(P[i][1], P[i][2], P[i][3], P[i][0]));
          break;
        default:
          break;
        }
      }
      count++;
      plot.add_event(e);
    }
    void finalize() {
      // thisrun.enu.Add(&enu);
      thisrun.count += count;
      for (const auto &[name, hist] : plot.histos) {
        thisrun.plot.histos[name]->Add(hist.get());
      }
    }
  };
  thread_object get_thread_object() {
    return thread_object(*this, plot.m_analysis_funcs, plot.m_cut_groups,
                         plot.m_cuts, plot.hist_2d_pairs);
  }
  void finalize() {
    // std::cout << "averaged overall xsec = " << xsec << " [cm^2]" <<
    // std::endl;

    // enu.Scale(xsec / event_count / bin_wid);

    for (const auto &[name, hist] : plot.histos) {
      // auto all = hist->Integral();

      // if (name.substr(0, 2) == "2d")
      //   hist->Scale(1. / all / bin_wid / bin_wid);
      // else
      //   hist->Scale(1. / all );
      hist->Scale(1. / count, "WIDTH");
    }
  }
  template <typename... Args>
  static run_manager_gibuu_pdk run_analysis(nlohmann::json &config,
                                            Args &&...args) {
    constexpr size_t MAX_COUNT = 128;
    chain_runner<int, int[MAX_COUNT], int[MAX_COUNT], double[MAX_COUNT][4]>
        chain(config["input_file_list"], "outtree",
              {"nparticles", "pdg", "status", "P"},
              std::thread::hardware_concurrency());

    return chain.run<run_manager_gibuu_pdk>(std::forward<Args>(args)...);
  }
};