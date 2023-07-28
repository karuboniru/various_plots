#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TMatrixT.h>
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

template <typename A, typename B, typename T>
double do_chi2(A &&data, B &&theory, T &&error_matrix) {
  assert(data.GetNbinsX() == theory.GetNbinsX());
  assert(data.GetNbinsX() == error_matrix.GetNrows());
  assert(data.GetNbinsX() == error_matrix.GetNcols());
  auto Nbins = data.GetNbinsX();
  std::decay_t<T> diff_vector(1, Nbins);
  for (int i = 0; i < Nbins; ++i) {
    diff_vector[0][i] = data.GetBinContent(i + 1) - theory.GetBinContent(i + 1);
  }
  auto diff_vector_transpose = diff_vector;
  diff_vector_transpose.Transpose(diff_vector_transpose);
  auto error_matrix_inverse = error_matrix.Invert();
  auto chi2 = diff_vector * error_matrix_inverse * diff_vector_transpose;
  return chi2[0][0];
}

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
  // double count1 = d.Count().GetValue();
  auto dataset =
      d.Filter(
           [](TObjString &EvtCode) {
             return getmode_nuwro(EvtCode) != event::channel::Other;
           },
           {"EvtCode"})
          .Define(
              "event",
              [](int StdHepN, ROOT::RVec<int> &StdHepPdg,
                 ROOT::RVec<int> &StdHepStatus, ROOT::RVec<double> &StdHepP4_,
                 TObjString &EvtCode) {
                double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                event e{};
                // if (getmode_nuwro(EvtCode) == event::channel::Other) {
                //   return e;
                // }
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
                  case 2:
                    // break;
                    // case 2:
                    {
                      TLorentzVector p4(StdHepP4[i][0], StdHepP4[i][1],
                                        StdHepP4[i][2], StdHepP4[i][3]);
                      e.add_particle_nofsi(pdg, p4);
                      if (pdg == 2112 || pdg == 2112) {
                        if (proton_count == 0) {
                          e.setprimaryP(p4);
                        }
                        if (proton_count == 1) {
                          e.setspectatorP(p4);
                        }
                        proton_count++;
                      }
                    }
                    break;
                  default:
                    break;
                  }
                }
                return e;
              },
              {"StdHepN", "StdHepPdg", "StdHepStatus", "StdHepP4", "EvtCode"});
  // apply cut from
  // https://journals.aps.org/prd/abstract/1.1103/PhysRevD.92.092008
  auto xsec = dataset.Mean("EvtWght").GetValue();
  auto count = dataset.Count().GetValue();
  auto dataset_cut =
      dataset
          // .Filter(
          //     [](event &e) {
          //       auto E = e.get_enu();
          //       return E > 1.5 && E < 10;
          //     },
          //     {"event"})
          // Single pi+/pi- in final state
          .Filter(
              [](event &e) {
                return e.count_out(211) + e.count_out(-211) == 1;
              },
              {"event"})
          // W_rest < 1.4
          .Filter(
              [](event &e) {
                auto W_rest = e.W_rest();
                return W_rest < 1.4;
              },
              {"event"})
          .
      // Save pion momentum
      Define("pion_mom",
             [](event &e) {
               TLorentzVector pion{};
               if (e.count_out(211) == 1) {
                 pion = e.get_particle_out(211).begin()->second;
               } else {
                 pion = e.get_particle_out(-211).begin()->second;
               }
               return pion;
             },
             {"event"})
          .
      // save pion kinetic energy
      Define("pion_Tk",
             [](TLorentzVector &pion) { return (pion.E() - pion.M()) * 1e3; },
             {"pion_mom"})
          // .Filter([](double pion_Tk) { return pion_Tk < 350 && pion_Tk > 0;
          // },
          //         {"pion_Tk"})
          .
      // save pion angle
      Define("pion_angle",
             [](TLorentzVector &pion) {
               //  auto pion_dir = pion.Vect().Unit();
               auto pion_angle = pion.Vect().Angle(TVector3{0, 0, 1.});
               // convert to degrees
               pion_angle *= 180 / M_PI;
               return pion_angle;
             },
             {"pion_mom"})
          .Define("W", [](event &e) { return e.getW_nofsi(); }, {"event"})
          .Define("Q2", [](event &e) { return e.getQ2(); }, {"event"})
          .Define("mp", [](event &e) { return e.getQ2(); }, {"event"})
          .Define("xbj",
                  [](ROOT::RVec<double> &StdHepP4_, double W, double Q2) {
                    double(*StdHepP4)[4] = (double(*)[4]) & StdHepP4_[0];
                    TLorentzVector p4(StdHepP4[1][0], StdHepP4[1][1],
                                      StdHepP4[1][2], StdHepP4[1][3]);
                    return Q2 / ((W * W - p4.M2()) + Q2);
                    // return Q2 / 2 /(p4.Dot(e.get_lvq()));
                  },
                  {"StdHepP4", "W", "Q2"})
          .Define("proton",
                  [](event &e) {
                    ROOT::RVec<TLorentzVector> p;
                    if (e.count_out(2212)) {
                      p.push_back(e.get_leading_out(2212));
                    }
                    return p;
                  },
                  {"event"})
          .Define("protonangle",
                  [](ROOT::RVec<TLorentzVector> &p) {
                    ROOT::RVec<double> var;
                    for (auto &k : p) {
                      var.push_back(k.Vect().Angle(TVector3{0, 0, 1.}) *
                                    (180 / M_PI));
                    }
                    return var;
                  },
                  {"proton"})
          .Define("protonmomentum",
                  [](ROOT::RVec<TLorentzVector> &p) {
                    ROOT::RVec<double> var;
                    for (auto &k : p) {
                      var.push_back(k.P());
                    }
                    return var;
                  },
                  {"proton"})
          .Define("muon", [](event &e) { return e.getPrimaryLepton(); },
                  {"event"})
          .Define("muonangle",
                  [](TLorentzVector &p) {
                    return p.Vect().Angle(TVector3{0, 0, 1.}) * (180 / M_PI);
                  },
                  {"muon"})
          .Define("muonmomentum", [](TLorentzVector &p) { return p.P(); },
                  {"muon"});
  // .Define("pion", [](event &e) { return e.getQ2(); }, {"event"});
  std::cout << "There is " << dataset_cut.Count().GetValue() << " events left"
            << "from " << count << " events" << std::endl;
  // Tk bins: 35–55 55–75 75–100 100–125 125–150 150–200 200–350
  const double Tk_bin_edges[]{35, 55, 75, 100, 125, 150, 200, 350};
  const double Tk_bin_minerva[]{1.13, 1.16, 1.07, .85, .76, .66, .38};
  const double Tk_bin_minerva_error[]{.30, .25, .20, .15, .14, .11, .08};
  const double Tk_bin_minerva_error_shape[]{.20, .12, .09, .06, .05, .05, .04};
  // constexpr double bins_TK = 7;
  const double cov_TK[7][7]{
      {1, .74, .72, .68, .68, .59, .56}, {.74, 1, .87, .82, .81, .72, .70},
      {.72, .87, 1, .85, .84, .76, .71}, {.68, .82, .85, 1, .88, .83, .79},
      {.68, .81, .84, .88, 1, .84, .81}, {.59, .72, .76, .83, .84, 1, .89},
      {.56, .70, .71, .79, .81, .89, 1}};

  const double cov_TK_shape[7][7]{{1, .29, .20, .01, -.02, -.30, -.36},
                                  {.29, 1, .39, .09, .02, -.4, -.47},
                                  {.20, .39, 1, .21, .13, -.22, -.53},
                                  {.01, .09, .21, 1, .25, .01, -.31},
                                  {-.02, .02, .13, .25, 1, .05, -.21},
                                  {-.30, -.4, -.22, .01, .05, 1, .27},
                                  {-.36, -.47, -.53, -.31, -.21, .27, 1}};

  // double full_error_matrix_TK[7][7];
  TMatrixT<double> full_error_matrix_TK(7, 7);
  TMatrixT<double> full_error_matrix_TK_shape(7, 7);
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      full_error_matrix_TK[i][j] =
          cov_TK[i][j] * Tk_bin_minerva_error[i] * Tk_bin_minerva_error[j];
      full_error_matrix_TK_shape[i][j] = cov_TK_shape[i][j] *
                                         Tk_bin_minerva_error_shape[i] *
                                         Tk_bin_minerva_error_shape[j];
    }
  }

  // auto TkInvert = full_error_matrix_TK.Invert();
  // auto TkInvert_shape = full_error_matrix_TK_shape.Invert();
  // const double

  TH1D h_Tk_data("h_Tk_data", "h_Tk_data", 7, Tk_bin_edges);
  for (int i = 0; i < 7; ++i) {
    h_Tk_data.SetBinContent(i + 1, Tk_bin_minerva[i]);
    h_Tk_data.SetBinError(i + 1, Tk_bin_minerva_error[i]);
  }
  // theta bins: 0–15 15–22 22–29 29–36 36–43 43–50 50–57 57–72 72–108 108–130
  // 130–140 140–150 150–165
  const double theta_bin_edges[]{0,  15, 22,  29,  36,  43,  50,
                                 57, 72, 108, 130, 140, 150, 165};
  const double theta_bin_minerva[]{1.83, 2.87, 3.05, 3.87, 3.54, 2.91, 2.13,
                                   1.98, 1.55, .90,  .71,  .54,  .33};
  const double theta_bin_minerva_error[]{.40, .57, .60, .77, .74, .61, .45,
                                         .40, .29, .19, .14, .11, .07};

  const double theta_bin_minerva_error_shape[]{
      .23, .26, .25, .29, .26, .23, .20, .19, .10, .11, .08, .06, .05};

  const double cov_theta[13][13]{
      {1, .78, .71, .71, .73, .71, .66, .65, .78, .73, .72, .69, .60},
      {.78, 1, .82, .82, .83, .82, .78, .76, .84, .74, .74, .73, .64},
      {.71, .82, 1, .87, .86, .86, .83, .81, .85, .72, .73, .73, .65},
      {.71, .82, .87, 1, .89, .88, .85, .83, .85, .73, .75, .74, .66},
      {.73, .83, .86, .89, 1, .90, .86, .83, .87, .75, .76, .76, .67},
      {.71, .82, .86, .88, .90, 1, .86, .83, .86, .74, .75, .75, .67},
      {.66, .78, .83, .85, .86, .86, 1, .82, .84, .72, .73, .73, .66},
      {.65, .76, .81, .83, .83, .83, .82, 1, .82, .72, .73, .72, .65},
      {.78, .84, .85, .85, .87, .86, .84, .82, 1, .80, .80, .79, .71},
      {.73, .74, .72, .73, .75, .74, .72, .72, .80, 1, .75, .73, .66},
      {.72, .74, .73, .75, .76, .75, .73, .73, .80, .75, 1, .76, .67},
      {.69, .73, .73, .74, .76, .75, .73, .72, .79, .73, .76, 1, .66},
      {.60, .64, .65, .66, .67, .67, .66, .65, .71, .66, .67, .66, 1}};

  const double cov_theta_shape[13][13]{
      {1, 0.19, -0.16, -0.24, -0.20, -0.23, -0.31, -0.26, 0.02, 0.13, 0.08,
       0.04, -0.03},
      {0.19, 1, 0.03, -0.02, -0.04, -0.07, -0.13, -0.13, -0.03, -0.02, -0.04,
       -0.04, -0.06},
      {-0.16, 0.03, 1, 0.16, 0.06, 0.07, 0.08, 0.07, -0.06, -0.16, -0.13, -0.11,
       -0.07},
      {-0.24, -0.02, 0.16, 1, 0.19, 0.18, 0.13, 0.09, -0.12, -0.21, -0.14,
       -0.12, -0.11},
      {-0.20, -0.04, 0.06, 0.19, 1, 0.24, 0.16, 0.02, -0.10, -0.14, -0.16,
       -0.12, -0.11},
      {-0.23, -0.07, 0.07, 0.18, 0.24, 1, 0.17, 0.07, -0.09, -0.14, -0.14,
       -0.11, -0.09},
      {-0.31, -0.13, 0.08, 0.13, 0.16, 0.17, 1, 0.13, -0.03, -0.11, -0.11,
       -0.07, -0.02},
      {-0.26, -0.13, 0.07, 0.09, 0.02, 0.07, 0.13, 1, -0.04, -0.05, -0.01,
       -0.02, 0.02},
      {0.02, -0.03, -0.06, -0.12, -0.10, -0.09, -0.03, -0.04, 1, 0.07, 0.04,
       0.05, 0.06},
      {0.13, -0.02, -0.16, -0.21, -0.14, -0.14, -0.11, -0.05, 0.07, 1, 0.17,
       0.13, 0.11},
      {0.08, -0.04, -0.13, -0.14, -0.16, -0.14, -0.11, -0.01, 0.04, 0.17, 1,
       0.23, 0.16},
      {0.04, -0.04, -0.11, -0.12, -0.12, -0.11, -0.07, -0.02, 0.05, 0.13, 0.23,
       1, 0.17},
      {-0.03, -0.06, -0.07, -0.11, -0.11, -0.09, -0.02, 0.02, 0.06, 0.11, 0.16,
       0.17, 1}};
  // double full_error_matrix_theta[13][13];
  TMatrixT<double> full_error_matrix_theta(13, 13);
  TMatrixT<double> full_error_matrix_theta_shape(13, 13);
  for (int i = 0; i < 13; ++i) {
    for (int j = 0; j < 13; ++j) {
      full_error_matrix_theta[i][j] = cov_theta[i][j] *
                                      theta_bin_minerva_error[i] *
                                      theta_bin_minerva_error[j];
      full_error_matrix_theta_shape[i][j] = cov_theta_shape[i][j] *
                                            theta_bin_minerva_error_shape[i] *
                                            theta_bin_minerva_error_shape[j];
    }
  }

  TH1D h_theta_data("h_theta_data", "h_theta_data", 13, theta_bin_edges);
  for (int i = 0; i < 13; ++i) {
    h_theta_data.SetBinContent(i + 1, theta_bin_minerva[i]);
    h_theta_data.SetBinError(i + 1, theta_bin_minerva_error[i]);
  }
  std::vector<ROOT::RDF::RResultPtr<TH1>> objs_list{};
  objs_list
      .emplace_back(
          dataset_cut.Histo1D({"h_tk", "h_tk", 7, Tk_bin_edges}, "pion_Tk"))
      ->Scale(xsec / count * 1e3, "WIDTH");
  auto &&h_Tk = objs_list.back();
  // auto &&h_Tk_nobinning = objs_list.back();
  auto int_Tk = h_Tk->Integral("WIDTH");
  auto int_Tk_data = h_Tk_data.Integral("WIDTH");
  auto h_Tk_shape = (TH1D *)(h_Tk->Clone("TK_shape"));
  h_Tk_shape->Scale(int_Tk_data / int_Tk);

  auto chi2_Tk = do_chi2(h_Tk_data, *h_Tk, full_error_matrix_TK);
  auto chi2_Tk_shape =
      do_chi2(h_Tk_data, *h_Tk_shape, full_error_matrix_TK_shape);

  draw_same(h_Tk, &h_Tk_data,
            ";pion kinetic energy (MeV);d#sigma/dT_{#pi} (10^{-41} cm^{2}/MeV)",
            chi2_Tk);
  draw_same(h_Tk_shape, &h_Tk_data,
            ";pion kinetic energy (MeV);d#sigma/dT_{#pi} (10^{-41} cm^{2}/MeV) "
            "shape only",
            chi2_Tk_shape);
  objs_list
      .emplace_back(dataset_cut.Histo1D(
          {"h_theta", "h_theta", 13, theta_bin_edges}, "pion_angle"))
      ->Scale(xsec / count * 1e3, "WIDTH");
  auto &&h_theta = objs_list.back();
  // auto &&h_theta_nobinning = objs_list.back();
  auto int_theta = h_theta->Integral("WIDTH");
  auto int_theta_data = h_theta_data.Integral("WIDTH");
  auto theta_shape = (TH1D *)h_theta->Clone("theta_shape");
  theta_shape->Scale(int_theta_data / int_theta);

  auto chi2_theta = do_chi2(h_theta_data, *h_theta, full_error_matrix_theta);
  auto chi2_theta_shape =
      do_chi2(h_theta_data, *theta_shape, full_error_matrix_theta_shape);

  draw_same(h_theta, &h_theta_data,
            ";pion angle (deg);d#sigma/d#theta (10^{-41} cm^{2}/deg)",
            chi2_theta);
  draw_same(
      theta_shape, &h_theta_data,
      ";pion angle (deg);d#sigma/d#theta (10^{-41} cm^{2}/deg) shape only",
      chi2_theta_shape);

  objs_list
      .emplace_back(dataset_cut.Histo1D(
          {"h_theta_nobinning", "h_theta_nobinning", 500, 0, 180.},
          "pion_angle"))
      ->Scale(xsec / count * 1e3, "WIDTH");
  objs_list
      .emplace_back(dataset_cut.Histo1D(
          {"h_Tk_nobinning", "h_Tk_nobinning", 500, 0, 700.}, "pion_Tk"))
      ->Scale(xsec / count * 1e3, "WIDTH");
  std::cout << "Plotting for H" << std::endl;
  objs_list
      .emplace_back(dataset_cut.Filter("StdHepPdg[1] == 2212")
                        .Histo1D({"h_theta_nobinning_H", "h_theta_nobinning_H",
                                  500, 0, 180.},
                                 "pion_angle"))
      ->Scale(xsec / count * 1e3 * 13, "WIDTH");
  objs_list
      .emplace_back(
          dataset_cut.Filter("StdHepPdg[1] == 2212")
              .Histo1D({"h_Tk_nobinning_H", "h_Tk_nobinning_H", 500, 0, 700.},
                       "pion_Tk"))
      ->Scale(xsec / count * 1e3 * 13, "WIDTH");

  save(objs_list, file);
  std::cout << "chi2_Tk = " << chi2_Tk << std::endl;
  std::cout << "chi2_Tk_shape = " << chi2_Tk_shape << std::endl;
  std::cout << "chi2_theta = " << chi2_theta << std::endl;
  std::cout << "chi2_theta_shape = " << chi2_theta_shape << std::endl;

  // W
  auto W_hist = dataset_cut.Histo1D(
      {"W", ";W (GeV); d#sigma/dW (10^{-38} cm^{2}/MeV) ", 100, .9, 2.}, "W");
  W_hist->Scale(xsec / count, "WIDTH");
  W_hist->SetMaximum(3.0);
  draw(W_hist, nullptr, "HIST");

  // Q2
  auto Q2hist = dataset_cut.Histo1D(
      {"Q2", ";Q2 (GeV^{2}); d#sigma/dQ^{2} (10^{-38} cm^{2}/GeV^{2}) ", 100, 0,
       10.},
      "Q2");
  Q2hist->Scale(xsec / count, "WIDTH");
  Q2hist->SetMaximum(0.42);
  draw(Q2hist, nullptr, "HIST");

  // xBj
  auto xbjhist = dataset_cut.Histo1D(
      {"xbj", ";xbj; d#sigma/dxbj (10^{-38} cm^{2}/xbj) ", 100, 0, 1.1}, "xbj");
  xbjhist->Scale(xsec / count, "WIDTH");
  xbjhist->SetMaximum(1.0);
  draw(xbjhist, nullptr, "HIST");

  // protonmomentum
  auto pphist = dataset_cut.Histo1D(
      {"protonmomentum", ";p_{p} (GeV); d#sigma/dp_{p} (10^{-38} cm^{2}/GeV) ",
       100, 0, 2.},
      "protonmomentum");
  pphist->Scale(xsec / count, "WIDTH");
  pphist->SetMaximum(0.4);
  draw(pphist, nullptr, "HIST");

  // protonangle
  auto pahist = dataset_cut.Histo1D(
      {"protonangle",
       ";#theta_{p} (deg); d#sigma/d#theta_{p} (10^{-38} cm^{2}/deg) ", 100, 0,
       180.},
      "protonangle");
  pahist->Scale(xsec / count, "WIDTH");
  pahist->SetMaximum(7e-3);
  draw(pahist, nullptr, "HIST");

  // muonmomentum
  auto muhist = dataset_cut.Histo1D(
      {"muonmomentum",
       ";p_{#mu} (GeV); d#sigma/dp_{#mu} (10^{-38} cm^{2}/GeV) ", 100, 0, 10.},
      "muonmomentum");
  muhist->Scale(xsec / count, "WIDTH");
  muhist->SetMaximum(0.12);
  draw(muhist, nullptr, "HIST");

  // muonangle
  auto muahist = dataset_cut.Histo1D(
      {"muonangle",
       ";#theta_{#mu} (deg); d#sigma/d#theta_{#mu} (10^{-38} cm^{2}/deg) ", 100,
       0, 180.},
      "muonangle");
  muahist->Scale(xsec / count, "WIDTH");
  muahist->SetMaximum(18e-3);
  draw(muahist, nullptr, "HIST");

  auto report = dataset_cut.Report();
  report->Print();
  return 0;
}