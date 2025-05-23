#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDFHelpers.hxx>
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
  ROOT::RDF::Experimental::AddProgressBar(d);
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
  try {
    dataset = dataset.Redefine("W", [](double W) { return W / 1000.; }, {"W"})
                  .Redefine("Q2", [](double Q2) { return Q2 / 1e6; }, {"Q2"});
  } catch (...) {
    std::cerr << "failed to redefine W" << std::endl;
  }
  // try {
  //   dataset =
  //       dataset.Define("W", [](event &e) { return e.getW_nofsi(); },
  //       {"event"})
  //           .Define("Q2", [](event &e) { return e.getQ2(); }, {"event"});
  // } catch (...) {
  // }
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
          // .Define("W", [](event &e) { return e.getW_nofsi(); }, {"event"})
          // .Define("Q2", [](event &e) { return e.getQ2(); }, {"event"})
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
  const double Tk_bin_minerva[]{.12e1, .13e1, .12e1, .09e1,
                                .08e1, .07e1, .04e1};
  const double Tk_bin_minerva_error[]{.26, .21, .18, .17, .18, .17, .22};
  // constexpr double bins_TK = 7;
  const double cov_TK[7][7]{
      {1.0000, 0.7178, 0.6969, 0.6691, 0.6598, 0.5719, 0.5313},
      {0.7178, 1.0000, 0.8679, 0.8090, 0.7962, 0.6995, 0.6617},
      {0.6969, 0.8679, 1.0000, 0.8266, 0.8137, 0.7249, 0.6626},
      {0.6691, 0.8090, 0.8266, 1.0000, 0.8823, 0.8166, 0.7830},
      {0.6598, 0.7962, 0.8137, 0.8823, 1.0000, 0.8295, 0.8024},
      {0.5719, 0.6995, 0.7249, 0.8166, 0.8295, 1.0000, 0.8985},
      {0.5313, 0.6617, 0.6626, 0.7830, 0.8024, 0.8985, 1.0000}};

  // double full_error_matrix_TK[7][7];
  TMatrixT<double> full_error_matrix_TK(7, 7);
  // TMatrixT<double> full_error_matrix_TK_shape(7, 7);
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      full_error_matrix_TK[i][j] = cov_TK[i][j] * Tk_bin_minerva_error[i] *
                                   Tk_bin_minerva_error[j] * Tk_bin_minerva[j] *
                                   Tk_bin_minerva[j];
    }
  }

  // auto TkInvert = full_error_matrix_TK.Invert();
  // auto TkInvert_shape = full_error_matrix_TK_shape.Invert();
  // const double

  TH1D h_Tk_data("h_Tk_data", "h_Tk_data", 7, Tk_bin_edges);
  for (int i = 0; i < 7; ++i) {
    h_Tk_data.SetBinContent(i + 1, Tk_bin_minerva[i]);
    h_Tk_data.SetBinError(i + 1, Tk_bin_minerva_error[i] * Tk_bin_minerva[i]);
  }
  // theta bins: 0–15 15–22 22–29 29–36 36–43 43–50 50–57 57–72 72–108 108–130
  // 130–140 140–150 150–165
  const double theta_bin_edges[]{0,  15,  22,  29,  36,  43,  50, 57,
                                 72, 108, 130, 140, 150, 165, 180};
  const double theta_bin_minerva[]{.12e1, .23e1, .28e1, .38e1, .36e1,
                                   .30e1, .23e1, .22e1, .17e1, .10e1,
                                   .08e1, .06e1, .04e1, .02e1};
  const double theta_bin_minerva_error[]{.23, .21, .20, .20, .20, .20, .21,
                                         .20, .19, .21, .19, .19, .21, .26};

  // const double theta_bin_minerva_error_shape[]{
  //     .23, .26, .25, .29, .26, .23, .20, .19, .10, .11, .08, .06, .05};

  const double cov_theta[14][14]{
      {1.0000, 0.8411, 0.8344, 0.8561, 0.8561, 0.8437, 0.8104, 0.7787, 0.8562,
       0.7050, 0.7002, 0.6942, 0.5834, 0.4134},
      {0.8411, 1.0000, 0.8468, 0.8617, 0.8606, 0.8438, 0.8254, 0.7975, 0.8330,
       0.6916, 0.6979, 0.6959, 0.6005, 0.4242},
      {0.8344, 0.8468, 1.0000, 0.8735, 0.8608, 0.8527, 0.8377, 0.8141, 0.8475,
       0.6781, 0.6873, 0.6856, 0.5911, 0.4354},
      {0.8561, 0.8617, 0.8735, 1.0000, 0.8897, 0.8810, 0.8574, 0.8376, 0.8573,
       0.7041, 0.7238, 0.7180, 0.6152, 0.4589},
      {0.8561, 0.8606, 0.8608, 0.8897, 1.0000, 0.8924, 0.8637, 0.8317, 0.8607,
       0.7207, 0.7273, 0.7261, 0.6309, 0.4594},
      {0.8437, 0.8438, 0.8527, 0.8810, 0.8924, 1.0000, 0.8564, 0.8264, 0.8589,
       0.7070, 0.7171, 0.7169, 0.6203, 0.4633},
      {0.8104, 0.8254, 0.8377, 0.8574, 0.8637, 0.8564, 1.0000, 0.8160, 0.8358,
       0.7020, 0.7147, 0.7176, 0.6325, 0.4621},
      {0.7787, 0.7975, 0.8141, 0.8376, 0.8317, 0.8264, 0.8160, 1.0000, 0.8095,
       0.7047, 0.7304, 0.7186, 0.6330, 0.5031},
      {0.8562, 0.8330, 0.8475, 0.8573, 0.8607, 0.8589, 0.8358, 0.8095, 1.0000,
       0.7619, 0.7509, 0.7447, 0.6452, 0.4745},
      {0.7050, 0.6916, 0.6781, 0.7041, 0.7207, 0.7070, 0.7020, 0.7047, 0.7619,
       1.0000, 0.7444, 0.7208, 0.6442, 0.4579},
      {0.7002, 0.6979, 0.6873, 0.7238, 0.7273, 0.7171, 0.7147, 0.7304, 0.7509,
       0.7444, 1.0000, 0.7534, 0.6670, 0.4832},
      {0.6942, 0.6959, 0.6856, 0.7180, 0.7261, 0.7169, 0.7176, 0.7186, 0.7447,
       0.7208, 0.7534, 1.0000, 0.6555, 0.4646},
      {0.5834, 0.6005, 0.5911, 0.6152, 0.6309, 0.6203, 0.6325, 0.6330, 0.6452,
       0.6442, 0.6670, 0.6555, 1.0000, 0.4535},
      {0.4134, 0.4242, 0.4354, 0.4589, 0.4594, 0.4633, 0.4621, 0.5031, 0.4745,
       0.4579, 0.4832, 0.4646, 0.4535, 1.0000}};

  TMatrixT<double> full_error_matrix_theta(14,14);
  for (int i = 0; i < 14; ++i) {
    for (int j = 0; j < 14; ++j) {
      full_error_matrix_theta[i][j] =
          cov_theta[i][j] * theta_bin_minerva_error[i] *
          theta_bin_minerva_error[j] * theta_bin_minerva[j] *
          theta_bin_minerva[j];
    }
  }

  TH1D h_theta_data("h_theta_data", "h_theta_data", 14, theta_bin_edges);
  for (int i = 0; i < 14; ++i) {
    h_theta_data.SetBinContent(i + 1, theta_bin_minerva[i]);
    h_theta_data.SetBinError(i + 1,
                             theta_bin_minerva_error[i] * theta_bin_minerva[i]);
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

  draw_same(h_Tk, &h_Tk_data,
            ";pion kinetic energy (MeV);d#sigma/dT_{#pi} (10^{-41} cm^{2}/MeV)",
            chi2_Tk);
  objs_list
      .emplace_back(dataset_cut.Histo1D(
          {"h_theta", "h_theta", 14, theta_bin_edges}, "pion_angle"))
      ->Scale(xsec / count * 1e3, "WIDTH");
  auto &&h_theta = objs_list.back();
  // auto &&h_theta_nobinning = objs_list.back();
  auto int_theta = h_theta->Integral("WIDTH");
  auto int_theta_data = h_theta_data.Integral("WIDTH");
  auto theta_shape = (TH1D *)h_theta->Clone("theta_shape");
  theta_shape->Scale(int_theta_data / int_theta);

  auto chi2_theta = do_chi2(h_theta_data, *h_theta, full_error_matrix_theta);

  draw_same(h_theta, &h_theta_data,
            ";pion angle (deg);d#sigma/d#theta (10^{-41} cm^{2}/deg)",
            chi2_theta);

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
      ->Scale(xsec / count * 1e4, "WIDTH");
  objs_list
      .emplace_back(
          dataset_cut.Filter("StdHepPdg[1] == 2212")
              .Histo1D({"h_Tk_nobinning_H", "h_Tk_nobinning_H", 500, 0, 700.},
                       "pion_Tk"))
      ->Scale(xsec / count * 1e4, "WIDTH");

  const double Wfactor = 1e-38;
  objs_list
      .emplace_back(
          dataset_cut.Histo1D({"Wtrue_all", "Wtrue_all", 500, .8, 3.}, "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  auto dataset_hydrogen = dataset_cut.Filter(
      [](ROOT::RVec<int> &StdHepPdg) { return StdHepPdg[1] == 2212; },
      {"StdHepPdg"});
  objs_list
      .emplace_back(dataset_hydrogen.Histo1D(
          {"Wtrue_hydrogen", "Wtrue_hydrogen", 500, .8, 3.}, "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  objs_list
      .emplace_back(dataset_hydrogen
                        .Filter(
                            [](const event &e, int flag_delta, double W) {
                              return e.get_mode() == event::channel::RES &&
                                     (flag_delta || W <= 1.210);
                            },
                            {"event", "flag_delta", "W"})
                        .Histo1D({"Wtrue_hydrogen_res", "Wtrue_hydrogen_res",
                                  500, .8, 3.},
                                 "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  objs_list
      .emplace_back(
          dataset_cut
              .Filter(
                  [](const event &e) {
                    return e.get_mode() != event::channel::RES &&
                           e.get_mode() != event::channel::DIS;
                  },
                  {"event"})
              .Histo1D({"Wtrue_nonres", "Wtrue_nonres", 500, .8, 3.}, "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  objs_list
      .emplace_back(
          dataset_cut
              .Filter(
                  [](const event &e) {
                    return e.get_mode() == event::channel::RES ||
                           e.get_mode() == event::channel::DIS;
                  },
                  {"event"})
              .Histo1D({"Wtrue_resdis", "Wtrue_resdis", 500, .8, 3.}, "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  objs_list
      .emplace_back(dataset_cut
                        .Filter(
                            [](const event &e, int flag_delta, double W) {
                              return e.get_mode() == event::channel::RES &&
                                     (flag_delta || W <= 1.210);
                            },
                            {"event", "flag_delta", "W"})
                        .Histo1D({"Wtrue_res", "Wtrue_res", 500, .8, 3.}, "W"))
      ->Scale(xsec / count * Wfactor, "WIDTH");

  save(objs_list, file);
  std::cout << "chi2_Tk = " << chi2_Tk << std::endl;
  // std::cout << "chi2_Tk_shape = " << chi2_Tk_shape << std::endl;
  std::cout << "chi2_theta = " << chi2_theta << std::endl;
  // std::cout << "chi2_theta_shape = " << chi2_theta_shape << std::endl;

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