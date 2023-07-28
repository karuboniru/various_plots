#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <event.h>
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
                   case 1: {
                     TLorentzVector p4(StdHepP4[i][0], StdHepP4[i][1],
                                       StdHepP4[i][2], StdHepP4[i][3]);
                     e.add_particle_out(pdg, p4);
                     if (pdg == 2212 || pdg == 2112) {
                       if (proton_count == 0) {
                         e.setprimaryP(p4);
                       }
                       if (proton_count == 1) {
                         e.setspectatorP(p4);
                       }
                       proton_count++;
                     }
                     break;
                   }
                   case 2: {
                     TLorentzVector p4(StdHepP4[i][0], StdHepP4[i][1],
                                       StdHepP4[i][2], StdHepP4[i][3]);
                     e.add_particle_nofsi(pdg, p4);
                     //  if (pdg == 2212 || pdg == 2112) {
                     //    if (proton_count == 0) {
                     //      e.setprimaryP(p4);
                     //    }
                     //    if (proton_count == 1) {
                     //      e.setspectatorP(p4);
                     //    }
                     //    proton_count++;
                     //  }
                     break;
                   }
                   default:
                     break;
                   }
                 }
                 return e;
               },
               {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"})
          .Define("W", [](event &e) { return e.getW(); }, {"event"})
          .Define("protonP", [](event &e) { return e.getprimaryP().P(); },
                  {"event"});
          // .Filter([](double protonP) { return protonP < 0.5; }, {"protonP"});
  auto count = dataset.Count().GetValue();
  std::cout << "Number of events: " << count << std::endl;
  auto xsec = d.Mean("EvtWght").GetValue();
  std::cout << "Cross section: " << xsec << " +- " << xsec / sqrt(count)
            << std::endl;
  std::vector<ROOT::RDF::RResultPtr<TH1>> objs_list{};
  auto add_plot = [&](auto &&dataset_, const char *name, const char *varname) {
    objs_list
        .emplace_back(dataset_.Histo1D({name, varname, 100, 0., 0.5}, varname))
        ->Scale(xsec / count, "WIDTH");
  };
  // add_plot(dataset, "W", "W");
  add_plot(dataset, "protonP", "protonP");
  // add_plot(dataset.Filter("enu>1 && enu < 10"), "kaonE_cut", "kaonE");
  save(objs_list, file);
  return 0;
}