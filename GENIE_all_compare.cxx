#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TVector3.h>
#include <TVectorDfwd.h>
#include <event.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <ostream>
#include <string>
#include <tools.h>
#include <unordered_set>
#include <utility>

using std::string_literals::operator""s;

std::size_t replace_all(std::string &inout, std::string_view what,
                        std::string_view with) {
  std::size_t count{};
  for (std::string::size_type pos{};
       inout.npos != (pos = inout.find(what.data(), pos, what.length()));
       pos += with.length(), ++count) {
    inout.replace(pos, what.length(), with.data(), with.length());
  }
  return count;
}

auto define_particle_spectrum(ROOT::RDF::RNode df, int pdgid) {
  auto pdg = TDatabasePDG::Instance()->GetParticle(pdgid);
  std::string name = pdg->GetName();
  auto orig_name = name;
  // auto pretty_title = name;
  // normalize pdgname to C++ name ("+"->"_p", "-"->"_m")
  replace_all(name, "+", "_p");
  replace_all(name, "-", "_m");
  // replace_all(pretty_title, "+", "^{+}");
  // replace_all(pretty_title, "-", "^{-}");
  // replace_all(pretty_title, "0", "^{0}");

  auto newdf =
      df.Define(name,
                [&](int StdHepN, const ROOT::RVec<int> &StdHepPdg,
                    const ROOT::RVec<int> &StdHepStatus,
                    const ROOT::RVec<double> &StdHepP4_) {
                  double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                  ROOT::RVec<TLorentzVector> particles{};
                  for (int i = 0; i < StdHepN; ++i) {
                    if (StdHepPdg[i] == pdgid && StdHepStatus[i] == 1) {
                      particles.emplace_back(StdHepP4[i][0], StdHepP4[i][1],
                                             StdHepP4[i][2], StdHepP4[i][3]);
                    }
                  }
                  return particles;
                },
                {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4"})
          .Define(name + "_E"s,
                  [](const ROOT::RVec<TLorentzVector> &particles) {
                    ROOT::RVec<double> energies{};
                    for (auto &p : particles) {
                      energies.push_back(p.E() - p.M());
                    }
                    return energies;
                  },
                  {name})
          .Define(name + "_P"s,
                  [](const ROOT::RVec<TLorentzVector> &particles) {
                    ROOT::RVec<double> energies{};
                    for (auto &p : particles) {
                      energies.push_back(p.P());
                    }
                    return energies;
                  },
                  {name})
          .Define(name + "theta",
                  [](const ROOT::RVec<TLorentzVector> &particles,
                     const TVector3 &p3_nu) {
                    ROOT::RVec<double> thetas{};
                    for (auto &p : particles) {
                      thetas.push_back(p.Vect().Angle(p3_nu));
                    }
                    return thetas;
                  },
                  {name, "p3_nu"})
          .Define(name + "phi",
                  [](const ROOT::RVec<TLorentzVector> &particles) {
                    ROOT::RVec<double> thetas{};
                    for (auto &p : particles) {
                      thetas.push_back(p.Vect().Phi());
                    }
                    return thetas;
                  },
                  {name});

  auto hist = newdf.Histo1D(
      {(name + "_E"s).c_str(), (orig_name).c_str(), 100, 0, 1}, name + "_E"s);
  auto hist_P = newdf.Histo1D(
      {(name + "_P"s).c_str(), (orig_name).c_str(), 100, 0, 1}, name + "_P"s);
  auto hist_theta =
      newdf.Histo1D({(name + "_theta"s).c_str(),
                     (orig_name + " #theta").c_str(), 100, 0, M_PI},
                    name + "theta");
  auto hist_phi = newdf.Histo1D(
      {(name + "_phi"s).c_str(), (orig_name + " #phi").c_str(), 100, 0., 0.},
      name + "phi");
  return std::make_tuple((TH1D *)hist->Clone(), (TH1D *)hist_theta->Clone(),
                         (TH1D *)hist_phi->Clone(), (TH1D *)hist_P->Clone());
}

int main(int argc, char **argv) {
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
  auto r = d.Take<TObjString, std::vector<TObjString>>("EvtCode");
  std::set<std::string> m{}; // set of unique event codes
  for (auto j : *(r.GetPtr())) {
    m.insert(std::string{j.GetString().View()});
  }

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
                  {"StdHepP4"});

  // std::cout << "Found " << pdglist.size() << " particles" << std::endl;

  auto spline_file =
      get_object<TGraph>(config["spline_file"].get<std::string>(),
                         config["spline_path"].get<std::string>());
  auto enu = df.Histo1D("enu");
  auto [tot, xsec] = get_xsec(enu.GetPtr(), spline_file.get());
  // std::cout << "Total cross section: " << xsec << " mb" << std::endl;
  // std::cout << "Total number of events: " << tot << std::endl;
  auto makeplots = [&](ROOT::RDF::RNode dataset,
                       std::string filename = "out.root") {
    auto outfile = std::make_unique<TFile>(filename.c_str(), "RECREATE");

    auto pdglist = dataset
                       .Reduce(
                           [](const std::unordered_set<int> &a,
                              const std::unordered_set<int> &b) {
                             auto c = a;
                             c.insert(b.begin(), b.end());
                             return c;
                           },
                           {"pdglist_out"})
                       .GetValue();
    for (auto pdgid : pdglist) {
      if (pdgid > 1000000000) {
        continue;
      }
      auto &&[hist, hist_theta, hist_phi, hist_P] =
          define_particle_spectrum(dataset, pdgid);
      hist->Scale(xsec / tot, "WIDTH");
      hist_theta->Scale(xsec / tot, "WIDTH");
      hist_phi->Scale(xsec / tot, "WIDTH");
      hist_P->Scale(xsec / tot, "WIDTH");
      hist->SetDirectory(outfile.get());
      hist_theta->SetDirectory(outfile.get());
      hist_phi->SetDirectory(outfile.get());
      hist_P->SetDirectory(outfile.get());
      // std::cout << hist->Integral("WIDTH") << std::endl;
    }
    outfile->Write();
    outfile->Close();
  };
  for (auto &code : m) {
    auto newdf = df.Filter(
        [&](TObjString &str) { return str.String() == code; }, {"EvtCode"});
    makeplots(newdf, code + ".root");
  }
  makeplots(df, "out.root");
  return 0;
}