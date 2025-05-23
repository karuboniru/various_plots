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
#include <iostream>
#include <memory>
#include <ostream>
#include <tools.h>

template <typename A, typename B, typename T>
double do_chi2(A &&data, B &&theory, T &&error_matrix) {
  assert(data.GetNbinsX() == theory.GetNbinsX());
  assert(data.GetNbinsX() == error_matrix.GetNrows());
  assert(data.GetNbinsX() == error_matrix.GetNcols());
  auto Nbins = data.GetNbinsX();
  std::decay_t<T> data_vec(1, Nbins), prediction_vec(1, Nbins);
  for (int i = 0; i < Nbins; ++i) {
    // diff_vector[0][i] = data.GetBinContent(i + 1) - theory.GetBinContent(i +
    // 1); std::cout << "data " << data.GetBinContent(i + 1) << " theory "
    //           << theory.GetBinContent(i + 1) << " diff " << diff_vector[0][i]
    //           << std::endl;
    data_vec[0][i] = data.GetBinContent(i + 1);
    prediction_vec[0][i] = theory.GetBinContent(i + 1);
  }
  // auto diff_vector_smeared = diff_vector * smear_matrix ;
  auto diff_vector_smeared = data_vec - prediction_vec;

  auto diff_vector_transpose = diff_vector_smeared;

  diff_vector_transpose.Transpose(diff_vector_transpose);
  auto error_matrix_inverse = error_matrix.Invert();
  auto chi2 =
      diff_vector_smeared * error_matrix_inverse * diff_vector_transpose;
  return chi2[0][0];
}

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
          // 1pi0, 0pi+-

          // .Filter(
          //     [](const event &e) {
          //       return e.count_out(13) == 1 && e.count_out(111) == 1;
          //     },
          //     {"event"}, "1mu1pi00pip")
          // // 0 other particle than proton and neutron
          .Filter(
              [](const event &e) {
                if (e.count_out(13) != 1 || e.count_out(111) != 1) {
                  return false;
                }
                for (auto &[pdg, count] : e.get_pdg_list_out()) {
                  switch (abs(pdg)) {
                  case 13:
                  case 111:
                    if (count != 1) {
                      return false;
                    }
                    break;

                  case 2212:
                  case 2112:
                    break;

                  default:
                    return false;
                  }
                }
                return true;
              },
              {"event"}, "topocut")
          .
      // Save pion momentum
      Define("pion_mom",
             [](const event &e) {
               return e.get_particle_out(111).begin()->second;
             },
             {"event"})
          .
      // save pion kinetic energy
      Define("pion_Tk",
             [](TLorentzVector &pion) { return (pion.E() - pion.M()) * 1e3; },
             {"pion_mom"})
          .Define("pion_p", [](TLorentzVector &pion) { return pion.P(); },
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

  const double pion_mom_bin_edges[]{0, .1, .15, .2, .3, .4, .5, .6, .799};
  const double pion_mom_bin_edges_with_overflow[]{0,  .1, .15, .2,  .3,
                                                  .4, .5, .6,  1e10};
  const double pion_mom_bin_content[]{0.66, 1.92, 2.44, 1.15,
                                      0.73, 0.34, 0.18, 0.09};
  const double pion_mom_cov_matrix[8][8]{
      {0.032, 0.015, 0.005, 0.006, 0.008, 0.007, 0.006, 0.004},
      {0.015, 0.023, 0.021, 0.011, 0.008, 0.009, 0.007, 0.004},
      {0.005, 0.021, 0.031, 0.024, 0.013, 0.009, 0.008, 0.005},
      {0.006, 0.011, 0.024, 0.035, 0.023, 0.009, 0.006, 0.006},
      {0.008, 0.008, 0.013, 0.023, 0.027, 0.014, 0.003, -0.000},
      {0.007, 0.009, 0.009, 0.009, 0.014, 0.016, 0.007, -0.002},
      {0.006, 0.007, 0.008, 0.006, 0.003, 0.007, 0.012, 0.009},
      {0.004, 0.004, 0.005, 0.006, -0.000, -0.002, 0.009, 0.017}};
  TMatrixT<double> full_error_matrix_pion_mom_raw(8, 8);
  full_error_matrix_pion_mom_raw.SetTol(1e-32);
  // TMatrixT<double> full_error_matrix_TK_shape(7, 7);
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      full_error_matrix_pion_mom_raw[i][j] =
          pion_mom_cov_matrix[i][j] / 1e4 /
          (pion_mom_bin_edges[j + 1] - pion_mom_bin_edges[j]) /
          (pion_mom_bin_edges[i + 1] - pion_mom_bin_edges[i]);
    }
  }

  const double pion_mom_smear_matrix[8][8]{
      {0.750, 0.320, 0.010, -0.027, -0.002, 0.014, 0.026, 0.024},
      {0.273, 0.393, 0.273, 0.025, -0.019, 0.012, 0.016, 0.012},
      {0.009, 0.282, 0.427, 0.229, 0.024, -0.014, 0.003, -0.004},
      {-0.061, 0.001, 0.266, 0.525, 0.284, -0.019, -0.084, -0.023},
      {-0.030, -0.040, 0.026, 0.259, 0.429, 0.239, -0.008, -0.069},
      {0.014, 0.011, -0.019, 0.006, 0.270, 0.477, 0.287, 0.016},
      {0.014, -0.001, -0.023, -0.070, 0.008, 0.273, 0.455, 0.313},
      {-0.004, -0.026, -0.057, -0.076, -0.132, -0.058, 0.270, 0.429}};
  TMatrixT<double> smear_matrix_pion_mom(8, 8);
  // full_error_matrix_pion_mom.SetTol(1e-32);
  // TMatrixT<double> full_error_matrix_TK_shape(7, 7);
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      smear_matrix_pion_mom[i][j] = pion_mom_smear_matrix[i][j];
    }
  }
  auto smear_matrix_pion_mom_T = smear_matrix_pion_mom;
  smear_matrix_pion_mom_T.Transpose(smear_matrix_pion_mom_T);
  // auto full_error_matrix_pion_mom_smear = smear_matrix_pion_mom_T *
  // full_error_matrix_pion_mom_raw * smear_matrix_pion_mom; auto
  // full_error_matrix_pion_mom_smear = smear_matrix_pion_mom *
  // full_error_matrix_pion_mom_raw * smear_matrix_pion_mom_T;
  auto full_error_matrix_pion_mom_smear = full_error_matrix_pion_mom_raw;
  full_error_matrix_pion_mom_smear.Print();
  TH1D h_pion_mom_data("h_pion_p_data", "h_pion_p_data", 8, pion_mom_bin_edges);
  for (int i = 0; i < 8; ++i) {
    h_pion_mom_data.SetBinContent(i + 1, pion_mom_bin_content[i] * 1e-1);
    h_pion_mom_data.SetBinError(i + 1,
                                sqrt(full_error_matrix_pion_mom_smear[i][i]));
  }

  std::vector<ROOT::RDF::RResultPtr<TH1>> objs_list{};

  objs_list.emplace_back(dataset_cut.Histo1D(
      {"h_pion_p_raw", "h_pion_p_raw", 8, pion_mom_bin_edges_with_overflow},
      "pion_p"));

  // ->Scale(xsec / count, "WIDTH");
  auto &&h_pion_p_with_overflow = objs_list.back();

  auto h_pion_p = new TH1D("h_pion_p", "h_pion_p", 8, pion_mom_bin_edges);
  for (int h_pion_p_bin_id = 0; h_pion_p_bin_id < 8; h_pion_p_bin_id++) {
    double sum = 0;
    for (int true_bin_id = 0; true_bin_id < 8; true_bin_id++) {
      sum += h_pion_p_with_overflow->GetBinContent(true_bin_id + 1) *
             //  pion_mom_smear_matrix[h_pion_p_bin_id][true_bin_id];
             pion_mom_smear_matrix[true_bin_id][h_pion_p_bin_id];
    }
    h_pion_p->SetBinContent(h_pion_p_bin_id + 1, sum);
    std::cout << "bin " << h_pion_p_bin_id << " " << sum << std::endl;
  }
  h_pion_p->Scale(xsec / count, "WIDTH");
  auto chi2_pion_p =
      do_chi2(*h_pion_p, h_pion_mom_data, full_error_matrix_pion_mom_smear);
  std::cout << "chi2 pion p " << chi2_pion_p << std::endl;

  const double Wfactor = 1e-38;
  objs_list
      .emplace_back(
          dataset_cut.Histo1D({"Wtrue_all", "Wtrue_all", 500, .8, 3.}, "W"))
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

  // W
  // auto W_hist = dataset_cut.Histo1D(
  //     {"W", ";W (GeV); d#sigma/dW (10^{-38} cm^{2}/MeV) ", 100, .9, 2.},
  //     "W");
  // W_hist->Scale(xsec / count, "WIDTH");
  // W_hist->SetMaximum(3.0);
  // draw(W_hist, nullptr, "HIST");

  // // Q2
  // auto Q2hist = dataset_cut.Histo1D(
  //     {"Q2", ";Q2 (GeV^{2}); d#sigma/dQ^{2} (10^{-38} cm^{2}/GeV^{2}) ", 100,
  //     0,
  //      10.},
  //     "Q2");
  // Q2hist->Scale(xsec / count, "WIDTH");
  // Q2hist->SetMaximum(0.42);
  // draw(Q2hist, nullptr, "HIST");

  // // xBj
  // auto xbjhist = dataset_cut.Histo1D(
  //     {"xbj", ";xbj; d#sigma/dxbj (10^{-38} cm^{2}/xbj) ", 100, 0, 1.1},
  //     "xbj");
  // xbjhist->Scale(xsec / count, "WIDTH");
  // xbjhist->SetMaximum(1.0);
  // draw(xbjhist, nullptr, "HIST");

  // // protonmomentum
  // auto pphist = dataset_cut.Histo1D(
  //     {"protonmomentum", ";p_{p} (GeV); d#sigma/dp_{p} (10^{-38} cm^{2}/GeV)
  //     ",
  //      100, 0, 2.},
  //     "protonmomentum");
  // pphist->Scale(xsec / count, "WIDTH");
  // pphist->SetMaximum(0.4);
  // draw(pphist, nullptr, "HIST");

  // // protonangle
  // auto pahist = dataset_cut.Histo1D(
  //     {"protonangle",
  //      ";#theta_{p} (deg); d#sigma/d#theta_{p} (10^{-38} cm^{2}/deg) ", 100,
  //      0, 180.},
  //     "protonangle");
  // pahist->Scale(xsec / count, "WIDTH");
  // pahist->SetMaximum(7e-3);
  // draw(pahist, nullptr, "HIST");

  // // muonmomentum
  // auto muhist = dataset_cut.Histo1D(
  //     {"muonmomentum",
  //      ";p_{#mu} (GeV); d#sigma/dp_{#mu} (10^{-38} cm^{2}/GeV) ", 100,
  //      0, 10.},
  //     "muonmomentum");
  // muhist->Scale(xsec / count, "WIDTH");
  // muhist->SetMaximum(0.12);
  // draw(muhist, nullptr, "HIST");

  // // muonangle
  // auto muahist = dataset_cut.Histo1D(
  //     {"muonangle",
  //      ";#theta_{#mu} (deg); d#sigma/d#theta_{#mu} (10^{-38} cm^{2}/deg) ",
  //      100, 0, 180.},
  //     "muonangle");
  // muahist->Scale(xsec / count, "WIDTH");
  // muahist->SetMaximum(18e-3);
  // draw(muahist, nullptr, "HIST");

  auto report = dataset_cut.Report();
  report->Print();
  return 0;
}