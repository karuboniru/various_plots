#pragma once
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <event.h>
#include <tools.h>

template <typename RDF> auto NuWro_RDF_setup_event(RDF &&df) {
  return df.Define(
      "event",
      [](int StdHepN, ROOT::RVec<int> &StdHepPdg, ROOT::RVec<int> &StdHepStatus,
         ROOT::RVec<double> &StdHepP4_, TObjString &EvtCode) {
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
        return e;
      },
      {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"});
}
