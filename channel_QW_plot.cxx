#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <event.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <tools.h>

// constexpr size_t Nbins = 300;
int main(int argc, const char **argv) {
  if (argc == 1)
    return 1;
  auto file = std::make_unique<TFile>("output.root", "RECREATE");
  ROOT::EnableImplicitMT();
  std::vector<std::string> names{};
  for (int i = 1; i < argc; ++i) {
    auto &name = argv[i];
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
                 bool lepton_set{false};
                 for (int i = 0; i < StdHepN; ++i) {
                   auto pdg = StdHepPdg[i];
                   if (StdHepPdg[i] == 1000000010) {
                     pdg = 2112;
                   }
                   auto thispartice = TLorentzVector(StdHepP4[i][0], StdHepP4[i][1],
                                             StdHepP4[i][2], StdHepP4[i][3]);
                   switch (StdHepStatus[i]) {
                   case 0:
                     e.add_particle_in(
                         pdg, thispartice);
                     break;
                   case 1:
                     e.add_particle_out(
                         pdg, thispartice);
                     if (!lepton_set){
                      lepton_set = true;
                      e.setPrimaryLepton(thispartice);
                     }
                     break;
                   case 2: {
                     e.add_particle_nofsi(pdg, thispartice);
                     if (pdg == 2212 || pdg == 2112) {
                       if (proton_count == 0) {
                         e.setprimaryP(thispartice);
                       }
                       if (proton_count == 1) {
                         e.setspectatorP(thispartice);
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
          .Define("W", [](event &e) { return e.getW_nofsi(); }, {"event"})
          .Define("Q2", [](event &e) { return e.getQ2(); }, {"event"});
  auto Wmax = dataset.Max<double>("W");
  // auto Qmax = dataset.Max("Q2");
  // auto W_plot_max = ((int)(Wmax.GetValue() * 2) + 1) / 2.;
  // auto Q_plot_max = ((int)(Qmax.GetValue() * 2) + 1) / 2.;
  auto dataset_qel = dataset.Filter(
      [](event &e) { return e.get_mode() == event::channel::QE; }, {"event"});
  auto dataset_res = dataset.Filter(
      [](event &e) { return e.get_mode() == event::channel::RES; }, {"event"});
  auto dataset_dis = dataset.Filter(
      [](event &e) { return e.get_mode() == event::channel::DIS; }, {"event"});
  auto dataset_2p2h = dataset.Filter(
      [](event &e) { return e.get_mode() == event::channel::MEC; }, {"event"});
  auto count = dataset.Count().GetValue();
  std::cout << "Number of events: " << count << std::endl;
  auto xsec = d.Mean("EvtWght").GetValue();
  std::cout << "Cross section: " << xsec << " +- " << xsec / sqrt(count)
            << std::endl;
  std::vector<ROOT::RDF::RResultPtr<TH1>> objs_list{};
  // auto enu = dataset.Mean<double>()
  auto add_plot = [&](auto &&dataset_, const char *name, const char *varname) {
    auto max = dataset_.Max(varname);
    double p_max{};
    std::string title{};
    switch (*varname) {
    case 'W':
      title = "d#sigma/d #it{W}; #it{W} (GeV); d#sigma/d #it{W} (#times "
              "10^{-38} cm^{2}/GeV)";
      p_max = 3.0;
      if (Wmax.GetValue() > 3.1) {
        p_max = 4.0;
      }
      break;
    case 'Q':
      p_max = 6.;
      if (Wmax.GetValue() > 3.1) {
        p_max = 10.0;
      }
      title = "d#sigma/d #it{Q^{2}}; #it{Q^{2}} (GeV^{2}); d#sigma/d "
              "#it{Q^{2}} (#times 10^{-38} cm^{2}/GeV^{2})";
      break;
    }
    objs_list
        .emplace_back(
            dataset_.Histo1D({name, title.c_str(), 200, 0., p_max}, varname))
        ->Scale(xsec / count, "WIDTH");
    auto &h1 = objs_list.emplace_back(dataset_.Histo1D(
        {(name + std::string{"_shape"}).c_str(), title.c_str(), 200, 0., p_max},
        varname));
    h1->Scale(1. / h1->Integral(), "WIDTH");
  };
  add_plot(dataset, "W", "W");
  add_plot(dataset_qel, "W_qel", "W");
  add_plot(dataset_res, "W_res", "W");
  add_plot(dataset_dis, "W_dis", "W");
  add_plot(dataset_2p2h, "W_2p2h", "W");
  add_plot(dataset, "Q2", "Q2");
  add_plot(dataset_qel, "Q2_qel", "Q2");
  add_plot(dataset_res, "Q2_res", "Q2");
  add_plot(dataset_dis, "Q2_dis", "Q2");
  add_plot(dataset_2p2h, "Q2_2p2h", "Q2");
  // add_plot(dataset.Filter("enu>1 && enu < 10"), "kaonE_cut", "kaonE");
  save(objs_list, file);
  return 0;
}