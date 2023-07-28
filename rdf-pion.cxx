#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <event.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <rdftools.h>
#include <tools.h>
constexpr size_t Nbins = 300;
int main(int argc, const char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <intput.json>" << std::endl;
    return 1;
  }
  nlohmann::json config;
  {
    std::ifstream ifs(argv[1]);
    if (!ifs.is_open()) {
      std::cerr << "Could not open file " << argv[1] << std::endl;
      return 1;
    }
    ifs >> config;
  }
  std::filesystem::path output_path = config["output_path"];
  if (!output_path.parent_path().empty()) {
    std::cout << "Creating output directory " << output_path.parent_path()
              << std::endl;
    std::filesystem::create_directories(output_path.parent_path());
  }
  auto file = std::make_unique<TFile>(output_path.c_str(), "RECREATE");
  ROOT::EnableImplicitMT();
  std::vector<std::string> names{};
  for (auto name : config["input_file_list"]) {
    // std::cout << "Adding file " << name << std::endl;
    names.push_back(name);
  }
  ROOT::RDataFrame d("nRooTracker", names);
  // d.Define("event", []() { return event{}; });
  auto dataset =
      d.Define("event",
               [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                  ROOT::RVec<int> &StdHepStatus, ROOT::RVec<double> &StdHepP4_,
                  TObjString &EvtCode) {
                 double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                 event e{};
                 if (getmode_nuwro(EvtCode) == event::channel::Other) {
                   return e;
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
                     e.add_particle_in(
                         pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                             StdHepP4[i][2], StdHepP4[i][3]));
                     break;
                   case 1:
                     e.add_particle_out(
                         pdg, TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                             StdHepP4[i][2], StdHepP4[i][3]));
                     break;
                   case 2: {
                     TLorentzVector p4(StdHepP4[i][0], StdHepP4[i][1],
                                       StdHepP4[i][2], StdHepP4[i][3]);
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
                 return e;
               },
               {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"})
          .Define("pion_count",
                  [](event &e) {
                    return e.count_out(211) + e.count_out(111) +
                           e.count_out(-211);
                  },
                  {"event"})
          .Define("pip_count", [](event &e) { return e.count_out(211); },
                  {"event"})
          .Define("pi0_count", [](event &e) { return e.count_out(111); },
                  {"event"})
          .Define("pim_count", [](event &e) { return e.count_out(-211); },
                  {"event"})
          .Define("pcount", [](event &e) { return e.count_out(2212); },
                  {"event"})
          .Define("ncount", [](event &e) { return e.count_out(2112); },
                  {"event"})
          .Define("out_count",
                  [](event &e) { return e.get_particle_out().size(); },
                  {"event"})
          // .Redefine("W", [](event &e) { return e.getW(); }, {"event"})
          // .Redefine("Q2", [](event &e) { return e.getQ2(); }, {"event"});
          .Redefine("W", [](double W) { return W / 1000.; }, {"W"})
          .Redefine("Q2", [](double Q2) { return Q2 / 1e6; }, {"Q2"});

  // spp_1: pi+ + p
  auto spp_1 =
      dataset.Filter("pip_count == 1 && pcount == 1 && out_count == 3");
  // spp_2: pi0 + p
  auto spp_2 =
      dataset.Filter("pi0_count == 1 && pcount == 1 && out_count == 3");
  // spp_3: pi+ + n
  auto spp_3 =
      dataset.Filter("pip_count == 1 && ncount == 1 && out_count == 3");
  auto xsec = d.Mean<double>("EvtWght").GetValue();
  xsec *= 1e-38;
  auto total = d.Count().GetValue();
  std::vector<ROOT::RDF::RResultPtr<TH1D>> objs_list{};
  objs_list.push_back(dataset.Histo1D("W"));
  objs_list.push_back(dataset.Histo1D("Q2"));
  objs_list.push_back(spp_1.Histo1D({"W_spp_1", "W_spp_1", Nbins, 0, 0}, "W"));
  objs_list.push_back(
      spp_1.Filter("flag_delta==1")
          .Histo1D({"W_spp_1_delta", "W_spp_1_delta", Nbins, 0, 0}, "W"));
  objs_list.push_back(spp_2.Histo1D({"W_spp_2", "W_spp_2", Nbins, 0, 0}, "W"));
  objs_list.push_back(
      spp_2.Filter("flag_delta==1")
          .Histo1D({"W_spp_2_delta", "W_spp_2_delta", Nbins, 0, 0}, "W"));
  objs_list.push_back(spp_3.Histo1D({"W_spp_3", "W_spp_3", Nbins, 0, 0}, "W"));
  objs_list.push_back(
      spp_3.Filter("flag_delta==1")
          .Histo1D({"W_spp_3_delta", "W_spp_3_delta", Nbins, 0, 0}, "W"));
  objs_list.push_back(
      spp_1.Histo1D({"Q2_spp_1", "Q2_spp_1", Nbins, 0, 0.}, "Q2"));
  objs_list.push_back(
      spp_2.Histo1D({"Q2_spp_2", "Q2_spp_2", Nbins, 0, 0.}, "Q2"));
  objs_list.push_back(
      spp_3.Histo1D({"Q2_spp_3", "Q2_spp_3", Nbins, 0, 0.}, "Q2"));

  objs_list.push_back(
      spp_1.Histo1D({"W_spp_1_detail", "W_spp_1", Nbins, 2.35, 2.45}, "W"));
  objs_list.push_back(
      spp_1.Filter("flag_delta==1")
          .Histo1D({"W_spp_1_delta_detail", "W_spp_1_delta", Nbins, 2.35, 2.45},
                   "W"));
  objs_list.push_back(
      spp_2.Histo1D({"W_spp_2_detail", "W_spp_2", Nbins, 2.35, 2.45}, "W"));
  objs_list.push_back(
      spp_2.Filter("flag_delta==1")
          .Histo1D({"W_spp_2_delta_detail", "W_spp_2_delta", Nbins, 2.35, 2.45},
                   "W"));
  objs_list.push_back(
      spp_3.Histo1D({"W_spp_3_detail", "W_spp_3", Nbins, 2.35, 2.45}, "W"));
  objs_list.push_back(
      spp_3.Filter("flag_delta==1")
          .Histo1D({"W_spp_3_delta_detail", "W_spp_3_delta", Nbins, 2.35, 2.45},
                   "W"));

  auto QW_spp_1 =
      spp_1.Histo2D({"QW_spp_1", "QW_spp_1; Q^{2} (GeV^{2}); W (GeV)", Nbins, 0,
                     10, Nbins, 1.08, 3.5},
                    "Q2", "W");
  auto model_log_x = (TH2D *)QW_spp_1->Clone();
  BinLogX(model_log_x->GetXaxis());
  auto QW_spp_1_logx = spp_1.Histo2D(*model_log_x, "Q2", "W");
  QW_spp_1_logx->SetName("QW_spp_1_logx");
  QW_spp_1_logx->SetTitle("QW_spp_1; Q^{2} (GeV^{2}); W (GeV)");
  auto QW_spp_2 =
      spp_2.Histo2D({"QW_spp_2", "QW_spp_2; Q^{2} (GeV^{2}); W (GeV)", Nbins, 0,
                     10, Nbins, 1.08, 3.5},
                    "Q2", "W");
  auto QW_spp_2_logx = spp_2.Histo2D(*model_log_x, "Q2", "W");
  QW_spp_2_logx->SetName("QW_spp_2_logx");
  QW_spp_2_logx->SetTitle("QW_spp_2; Q^{2} (GeV^{2}); W (GeV)");
  auto QW_spp_3 =
      spp_3.Histo2D({"QW_spp_3", "QW_spp_3; Q^{2} (GeV^{2}); W (GeV)", Nbins, 0,
                     10, Nbins, 1.08, 3.5},
                    "Q2", "W");
  auto QW_spp_3_logx = spp_3.Histo2D(*model_log_x, "Q2", "W");
  QW_spp_3_logx->SetName("QW_spp_3_logx");
  QW_spp_3_logx->SetTitle("QW_spp_3; Q^{2} (GeV^{2}); W (GeV)");

  auto plot_hist = [&](auto &&tt) {
    ResetStyle(tt);
    tt->Scale(xsec / total, "WIDTH");
    auto c1 = getCanvas();
    tt->Draw("COLZ");
    tt->SetLineWidth(2);
    c1->SaveAs((tt->GetName() + std::string(".pdf")).c_str());
    c1->SaveAs((tt->GetName() + std::string(".png")).c_str());
    c1->SaveAs((tt->GetName() + std::string(".eps")).c_str());
    c1->SetLogz();
    c1->SaveAs((tt->GetName() + std::string("_log.pdf")).c_str());
    c1->SaveAs((tt->GetName() + std::string("_log.png")).c_str());
    c1->SaveAs((tt->GetName() + std::string("_log.eps")).c_str());
  };

  auto plot_2d_with_normalized = [&](TH2D &hist) {
    auto nor_x = normalize_slice(&hist, true);
    auto nor_y = normalize_slice(&hist, false);
    plot_hist(nor_x);
    plot_hist(nor_y);
    plot_hist(&hist);
  };

  plot_2d_with_normalized(*QW_spp_1);
  plot_2d_with_normalized(*QW_spp_2);
  plot_2d_with_normalized(*QW_spp_3);
  plot_2d_with_normalized(*QW_spp_1_logx);
  plot_2d_with_normalized(*QW_spp_2_logx);
  plot_2d_with_normalized(*QW_spp_3_logx);

  for (auto &obj : objs_list) {
    // auto width = obj->GetBinWidth(1);
    obj->Scale(xsec / total, "WIDTH");
  }
  save(objs_list, file);
  return 0;
}