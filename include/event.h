#pragma once

// #include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <array>
#include <map>
#include <set>
#include <type_traits>
#include <vector>
template <typename T> class eq_range {
private:
  T pair;

public:
  template <typename T_> eq_range(T_ &&p) : pair(p) {}
  decltype(auto) begin() { return pair.first; };
  decltype(auto) end() { return pair.second; };
};
template <typename T_> eq_range(T_ &&) -> eq_range<std::remove_reference_t<T_>>;
class event {
public:
  enum class channel { QE, RES, DIS, MEC, Other };

private:
  std::unordered_multimap<int, TLorentzVector> in_particles{};
  std::unordered_multimap<int, TLorentzVector> out_particles{};
  size_t out_count{};
  // std::set<int> pdg_list_out{};
  std::map<int, size_t> pdg_list_out{};
  std::unordered_multimap<int, TLorentzVector> nofsi_particles{};
  double weight{1};
  channel mode;
  // TDatabasePDG db{};
  std::string channelname{};
  std::string channelname_nonucleon{};
  std::array<int, 18> m_nod{}; ///< number of rescattering interactions of given
                               ///< type:
                               ///< 0 - nucleon elastic,
                               ///< 1 - nucleon ce,
                               ///< 2 - nucleon spp,
                               ///< 3 - nucleon dpp,
                               ///< 4 - pion elastic,
                               ///< 5 - pion ce,
                               ///< 6 - pion spp,
                               ///< 7 - pion dpp,
                               ///< 8 - pion abs,
                               ///< 9 - jailed,
                               ///< 10 - escape,
                               ///< 11 - pion tpp,
                               ///< 12 - pion no interaction,
                               ///< 13 - nucleon no interaction,
                               ///< 14 - hyperon no interaction,
                               ///< 15 - hyperon elastic scatter,
                               ///< 16 - hyperon Lambda -> Sigma conversion,
                               ///< 17 - hyperon Sigma -> Lambda conversion,
  TLorentzVector primaryP{}, spectatorP{};
  TLorentzVector primarylepton{};
  bool is_spectator{false};
  bool found_lepton{false};

public:
  event() {}
  ~event();
  void setprimaryP(const TLorentzVector &p) { primaryP = p; }
  void setspectatorP(const TLorentzVector &p) {
    spectatorP = p;
    is_spectator = true;
  }
  size_t get_out_count() const { return out_count; }
  const TLorentzVector &getprimaryP() const { return primaryP; }
  const TLorentzVector &getspectatorP() const { return spectatorP; }
  bool is_spectator_event() const { return is_spectator; }
  double getQ2() const;
  double getW() const;
  void add_particle_in(int id, const TLorentzVector &p4);
  void add_particle_out(int id, const TLorentzVector &p4);
  void add_particle_nofsi(int id, const TLorentzVector &p4);
  bool is_good_event() const;
  std::pair<double, double> get_q2_w() const;
  double get_enu() const;
  TLorentzVector get_leading_proton() const;
  TLorentzVector get_leading_out(int) const;
  TLorentzVector get_leading_nofsi(int) const;
  auto get_particle_out(int pdgid) {
    return eq_range{out_particles.equal_range(pdgid)};
  }
  auto get_particle_in(int pdgid) {
    return eq_range{in_particles.equal_range(pdgid)};
  }
  auto get_particle_nofsi(int pdgid) {
    return eq_range{nofsi_particles.equal_range(pdgid)};
  }
  size_t count_particle_out(int pdgid) const noexcept;
  size_t count_particle_nofsi(int pdgid) const noexcept;
  bool TKI_mu_p_cut() const;
  bool TKI_mu_cut() const;
  void set_mode(channel);
  channel get_mode() const;
  void set_weight(double);
  double get_weight() const;
  const std::unordered_multimap<int, TLorentzVector> &get_particle_out() const {
    return out_particles;
  };
  const std::unordered_multimap<int, TLorentzVector> &get_particle_in() const {
    return in_particles;
  };
  const std::unordered_multimap<int, TLorentzVector> &
  get_particle_nofsi() const {
    return nofsi_particles;
  };
  size_t count_out(int) const noexcept;
  const std::string &get_channelname();
  const std::string &get_channelname_no_nucleon();
  void set_nod(const std::array<int, 18> &nod);
  bool is_no_pion_interaction() const;
  bool is_true_elastic() const;
  const std::array<int, 18> &get_nod() const;
  int get_pion_interaction_count() const;
  double getW_nofsi() const;
  void setPrimaryLepton(const TLorentzVector &p) { primarylepton = p; }
  const TLorentzVector & getPrimaryLepton() { return primarylepton; }
  double W_rest() const;
  TLorentzVector get_lvq() const;
};
