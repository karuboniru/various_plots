#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TVector3.h>
#include <event.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <tools.h>
#include <type_traits>

auto get_kin_variables(ROOT::RDF::RNode input_node) {
  auto df = input_node.Define(
      "particle_selected_p",
      [&](ROOT::RVec<double> &StdHepP4_, std::set<int> &pdg_status_to_index) {
        double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
        ROOT::RVec<double> particle_selected_p{};
        for (auto &index : pdg_status_to_index) {
          TVector3 p3{StdHepP4[index][0], StdHepP4[index][1],
                      StdHepP4[index][2]};
          particle_selected_p.emplace_back(p3.Mag());
        }
        return particle_selected_p;
      },
      {"StdHepP4", "current_id"});
  auto hist = df.Histo1D({"", "", 100, 0, 1.}, "particle_selected_p");
  return (TH1D *)hist->Clone();
}

auto get_kin_variables_final(ROOT::RDF::RNode input_node) {
  auto df = input_node.Define(
      "particle_selected_p",
      [&](ROOT::RVec<double> &StdHepP4_, std::set<int> &pdg_status_to_index,
          std::map<int, int> &id_history) {
        double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
        ROOT::RVec<double> particle_selected_p{};
        for (auto &thisindex : pdg_status_to_index) {
          auto index = id_history[thisindex];
          if (index == -1 || index == 0)
            continue;
          TVector3 p3{StdHepP4[index][0], StdHepP4[index][1],
                      StdHepP4[index][2]};
          particle_selected_p.emplace_back(p3.Mag());
        }
        return particle_selected_p;
      },
      {"StdHepP4", "current_id", "id_history"});
  auto hist = df.Histo1D({"", "", 100, 0, 1.}, "particle_selected_p");
  return (TH1D *)hist->Clone();
}

ROOT::RDF::RNode re_define_current_id(ROOT::RDF::RNode df, int pdg_to_get,
                                      int status_to_get) {
  return df
      .Redefine(
          "current_id",
          [&](std::map<std::pair<int, int>, std::set<int>> &next_level_ids) {
            std::set<int> current_id_new{};
            for (auto &index : next_level_ids[{pdg_to_get, status_to_get}]) {
              current_id_new.insert(index);
            }
            return current_id_new;
          },
          {"next_level_ids"})
      .Redefine("next_level_keys",
                [](std::set<int> &current_id_old, ROOT::RVec<int> &StdHepPdg,
                   ROOT::RVec<int> &StdHepStatus, ROOT::RVec<int> &StdHepFm) {
                  std::set<std::pair<int, int>> next_level_keys;
                  for (auto &index : current_id_old) {
                    auto m_index = StdHepFm[index];
                    if (m_index != -1)
                      next_level_keys.insert(
                          {StdHepPdg[m_index], StdHepStatus[m_index]});
                  }
                  return next_level_keys;
                },
                {"current_id", "StdHepPdg", "StdHepStatus", "StdHepFm"})
      .Redefine(
          "next_level_ids",
          [](std::set<int> &current_id_old, ROOT::RVec<int> &StdHepPdg,
             ROOT::RVec<int> &StdHepStatus, ROOT::RVec<int> &StdHepFm) {
            std::map<std::pair<int, int>, std::set<int>> next_level_ids;
            for (auto &index : current_id_old) {
              auto m_index = StdHepFm[index];
              if (m_index != -1)
                next_level_ids[{StdHepPdg[m_index], StdHepStatus[m_index]}]
                    .insert(m_index);
            }
            return next_level_ids;
          },
          {"current_id", "StdHepPdg", "StdHepStatus", "StdHepFm"})
      .Redefine(
          "id_history",
          [](std::map<int, int> &id_history_old, ROOT::RVec<int> &StdHepFm) {
            std::map<int, int> id_history;
            for (auto &[hist, final] : id_history_old) {
              auto m_index = StdHepFm[hist];
              if (m_index != -1)
                id_history[m_index] = final;
            }
            return id_history;
          },
          {"id_history", "StdHepFm"});
}

std::set<std::pair<int, int>> get_all_key_list(ROOT::RDF::RNode df) {
  return df
      .Reduce(
          [](std::set<std::pair<int, int>> a, std::set<std::pair<int, int>> b) {
            std::set<std::pair<int, int>> ret{};
            // merge two sets
            std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                           std::inserter(ret, ret.begin()));
            return ret;
          },
          "next_level_keys")
      .GetValue();
}

void do_plot(ROOT::RDF::RNode df, TDirectory *dir, double normalize = 1.,
             int pdg = 0, int status = 0) {
  auto hist = get_kin_variables(df);
  auto hist1 = get_kin_variables_final(df);
  hist->SetName(Form("d%d_%d", pdg, status));
  hist->SetTitle(Form("d%d_%d", pdg, status));
  hist1->SetNameTitle("proton", "proton");
  hist->Scale(normalize, "width");
  hist1->Scale(normalize, "width");
  std::cout << "get hist of xsec: " << hist->Integral("width") << std::endl;
  if (hist->Integral("width") < 0.1) {
    std::cout << "skip this hist" << std::endl;
    return;
  }
  // hist->SetDirectory(dir);
  // dir->Add(hist);
  // dir->cd();
  hist->SetDirectory(dir);
  hist1->SetDirectory(dir);
  // hist->Write();
  auto list_of_next_key = get_all_key_list(df);
  for (auto &key : list_of_next_key) {
    std::cout << "pdg: " << key.first << " status: " << key.second << std::endl;
    auto next_df = re_define_current_id(df, key.first, key.second);
    auto next_dir = dir->mkdir(Form("d%d_%d", key.first, key.second));
    do_plot(next_df, next_dir, normalize, key.first, key.second);
  }
  std::cout << "end of do_plot" << std::endl;
}

int main(int argc, char **argv) {
  TH1::AddDirectory(kFALSE);
  if (argc == 1) {
    std::cout << "Usage: " << argv[0] << " <normalize.json> <ROOT files>"
              << std::endl;
    return 1;
  }
  ROOT::EnableImplicitMT();
  std::vector<std::string> names{};
  for (int i = 2; i < argc; ++i) {
    names.push_back(argv[i]);
  }
  nlohmann::json config;
  {
    std::string filename = argv[1];
    std::fstream file(filename);
    file >> config;
  }

  ROOT::RDataFrame d("gRooTracker", names);

  auto df =
      d.Define("pdglist_out",
               [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                  ROOT::RVec<int> &StdHepStatus) {
                 std::unordered_set<int> pdglist;
                 for (int i = 0; i < StdHepN; ++i) {
                   if (StdHepStatus[i] == 1) {
                     pdglist.insert(StdHepPdg[i]);
                   }
                 }
                 return pdglist;
               },
               {"StdHepN", "StdHepPdg", "StdHepStatus"})
          .Define("p3_nu",
                  [](ROOT::RVec<double> &StdHepP4_) {
                    return TVector3{StdHepP4_[0], StdHepP4_[1], StdHepP4_[2]};
                  },
                  {"StdHepP4"})
          .Define("enu",
                  [](ROOT::RVec<double> &StdHepP4_) { return StdHepP4_[3]; },
                  {"StdHepP4"})
          .Define("current_id",
                  [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                     ROOT::RVec<int> &StdHepStatus) {
                    std::set<int> current_id;
                    for (int i = 0; i < StdHepN; ++i) {
                      if (StdHepStatus[i] == 1) {
                        if (StdHepPdg[i] == 2212) {
                          current_id.insert(i);
                        }
                      }
                    }
                    return current_id;
                  },
                  {"StdHepN", "StdHepPdg", "StdHepStatus"})
          // .Define("current_final_id",
          //         [](std::set<int> current_id) { return current_id; },
          //         {"current_id"})
          .Define("next_level_keys",
                  [](std::set<int> &current_id_old, ROOT::RVec<int> &StdHepPdg,
                     ROOT::RVec<int> &StdHepStatus, ROOT::RVec<int> &StdHepFm) {
                    std::set<std::pair<int, int>> next_level_keys;
                    for (auto &index : current_id_old) {
                      auto m_index = StdHepFm[index];
                      if (m_index != -1)
                        next_level_keys.insert(
                            {StdHepPdg[m_index], StdHepStatus[m_index]});
                    }
                    return next_level_keys;
                  },
                  {"current_id", "StdHepPdg", "StdHepStatus", "StdHepFm"})
          .Define(
              "next_level_ids",
              [](std::set<int> &current_id_old, ROOT::RVec<int> &StdHepPdg,
                 ROOT::RVec<int> &StdHepStatus, ROOT::RVec<int> &StdHepFm) {
                std::map<std::pair<int, int>, std::set<int>> next_level_ids;
                for (auto &index : current_id_old) {
                  auto m_index = StdHepFm[index];
                  if (m_index != -1)
                    next_level_ids[{StdHepPdg[m_index], StdHepStatus[m_index]}]
                        .insert(m_index);
                }
                return next_level_ids;
              },
              {"current_id", "StdHepPdg", "StdHepStatus", "StdHepFm"})
          .Define("id_history",
                  [](std::set<int> &current_id_old, ROOT::RVec<int> &StdHepPdg,
                     ROOT::RVec<int> &StdHepStatus, ROOT::RVec<int> &StdHepFm) {
                    std::map<int, int> id_history;
                    for (auto &index : current_id_old) {
                      // auto m_index = StdHepFm[index];
                      // if (m_index != -1)
                      id_history[index] = index;
                    }
                    return id_history;
                  },
                  {"current_id", "StdHepPdg", "StdHepStatus", "StdHepFm"});

  auto spline_file =
      get_object<TGraph>(config["spline_file"].get<std::string>(),
                         config["spline_path"].get<std::string>());
  auto enu = df.Histo1D("enu");

  auto [tot, xsec] = get_xsec(enu.GetPtr(), spline_file.get());
  auto normalize_factor = xsec / tot;
  // auto list_states = get_all_key_list(df);
  auto file = TFile::Open("test.root", "RECREATE");
  auto dir = file->mkdir("plots");
  do_plot(df.Filter(
              [](TObjString &EvtCode) {
                return EvtCode.GetString().Contains(
                    "Weak[CC],QES"); // Just because I want to make things quick
              },
              {"EvtCode"}, "QE-CUT"),
          dir, normalize_factor, 2212, 1);
  file->Write();
  file->Close();
  return 0;
}