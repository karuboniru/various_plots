#include <cassert>
#include <event.h>
#include <sstream>

template <typename T> inline TLorentzVector get_sum(T &&cont, int pdgid) {
  TLorentzVector sum{};
  for (const auto &[id, p4] : eq_range{cont.equal_range(pdgid)}) {
    sum += p4;
  }
  return sum;
}

template <typename T> inline TLorentzVector get_leading(T &&cont, int pdgid) {
  const TLorentzVector *max{nullptr};
  double max_E{-INFINITY};
  for (const auto &[id, p4] : eq_range{cont.equal_range(pdgid)}) {
    if (p4.E() > max_E) {
      max_E = p4.E();
      max = &p4;
    }
  }
  if(max){
    return *max;
  } else {
    return TLorentzVector{};
  }
}

event::~event() {}

double event::getQ2() const {
  // const auto &p4_muon = out_particles.at(13);
  const auto &p4_muon = primarylepton;
  const auto &p4_neutrino = in_particles.find(14)->second;
  auto q = p4_muon - p4_neutrino;
  return -q.M2();
}

double event::getW() const {
  TLorentzVector p_final_state{};
  for (const auto &[pdgid, p4v] : in_particles) {
    p_final_state += p4v;
  }
  p_final_state -= primarylepton;
  return p_final_state.M();
}

double event::getW_nofsi() const {
  TLorentzVector p_final_state{};
  for (const auto &[pdgid, p4v] : nofsi_particles) {
    p_final_state += p4v;
  }
  p_final_state -= primarylepton;
  return p_final_state.M();
}

void event::add_particle_in(int id, const TLorentzVector &p4) {
  in_particles.emplace(id, p4);
}

void event::add_particle_out(int id, const TLorentzVector &p4) {
  out_particles.emplace(id, p4);
  // pdg_list_out.insert(id);
  pdg_list_out[id]++;
  if ((!found_lepton) && id == 13) {
    primarylepton = p4;
    found_lepton = true;
  }
  out_count++;
}

void event::add_particle_nofsi(int id, const TLorentzVector &p4) {
  nofsi_particles.emplace(id, p4);
}

std::pair<double, double> event::get_q2_w() const {
  return std::make_pair(getQ2(), getW());
}

double event::get_enu() const {
  assert(in_particles.count(14) == 1);
  const auto &p4_neutrino = in_particles.find(14)->second;
  return p4_neutrino.E();
}

bool event::TKI_mu_p_cut() const {
  if (count_out(13) == 0 || count_out(2212) == 0) {
    return false;
  }
  const auto leading_p4_muon = get_leading(out_particles, 13);
  const auto p_mu = leading_p4_muon.P();
  const auto angle_mu = leading_p4_muon.Theta() * 180 / M_PI;
  const auto leading_p4_proton = get_leading(out_particles, 2212);
  // ignore any proton cut
  // const auto p_proton = leading_p4_proton.P();
  const auto angle_proton = leading_p4_proton.Theta() * 180 / M_PI;
  // return 1.5 < p_mu && p_mu < 20. && angle_mu < 25. && 0.45 < p_proton;
  return 1.5 < p_mu && p_mu < 20. && angle_mu < 25. && angle_proton < 70;
}

TLorentzVector event::get_leading_proton() const {
  return get_leading(out_particles, 2212);
}

size_t event::count_particle_out(int pdgid) const noexcept {
  return count_out(pdgid);
}

size_t event::count_particle_nofsi(int pdgid) const noexcept {
  return nofsi_particles.count(pdgid);
}

TLorentzVector event::get_leading_out(int pdgid) const {
  return get_leading(out_particles, pdgid);
}

TLorentzVector event::get_leading_nofsi(int pdgid) const {
  return get_leading(nofsi_particles, pdgid);
}

bool event::TKI_mu_cut() const {
  if (count_out(13) == 0 || count_out(2212) == 0) {
    return false;
  }
  const auto leading_p4_muon = get_leading(out_particles, 13);
  const auto p_mu = leading_p4_muon.P();
  const auto angle_mu = leading_p4_muon.Theta() * 180 / M_PI;
  return 1.5 < p_mu && p_mu < 20. && angle_mu < 25.;
}

void event::set_mode(event::channel m) { mode = m; }

event::channel event::get_mode() const { return mode; }

void event::set_weight(double w) { weight = w; }

double event::get_weight() const { return weight; }

size_t event::count_out(int id) const noexcept {
  auto it = pdg_list_out.find(id);
  if (it == pdg_list_out.end()) {
    return 0;
  }
  return it->second;
}

void event::set_nod(const std::array<int, 18> &nod) { m_nod = nod; }

bool event::is_no_pion_interaction() const {
  return !(m_nod[4] || m_nod[5] || m_nod[6] || m_nod[7] || m_nod[8]);
}

bool event::is_true_elastic() const {
  return !(m_nod[5] || m_nod[6] || m_nod[7] || m_nod[8]);
}

const std::array<int, 18> &event::get_nod() const { return m_nod; }

int event::get_pion_interaction_count() const { return m_nod[4] + m_nod[5]; }

#include <TDatabasePDG.h>
const std::string &event::get_channelname() {
  if (!channelname.empty())
    return channelname;
  std::stringstream ss{};
  thread_local TDatabasePDG db{};
  for (const auto &[id, count] : pdg_list_out) {
    if (count > 1) {
      ss << count;
    }
    ss << db.GetParticle(id)->GetName() << " ";
  }
  channelname = ss.str();
  channelname = channelname.substr(0, channelname.size() - 1);
  return channelname;
}

const std::string &event::get_channelname_no_nucleon() {
  if (!channelname_nonucleon.empty())
    return channelname_nonucleon;
  std::stringstream ss{};
  thread_local TDatabasePDG db{};
  // auto & db = * TDatabasePDG::Instance();
  for (const auto &[id, count] : pdg_list_out) {
    if (id == 2212 || id == 2112) {
      continue;
    }
    if (count > 1) {
      ss << count;
    }
    ss << db.GetParticle(id)->GetName() << " ";
  }
  channelname_nonucleon = ss.str();
  channelname_nonucleon =
      channelname_nonucleon.substr(0, channelname_nonucleon.size() - 1);
  return channelname_nonucleon;
}

double event::W_rest() const {
  constexpr double m_p = 0.938272;
  TLorentzVector dummyp{0., 0, 0, m_p};
  auto p4_neutrino = in_particles.find(14)->second;
  auto lvq = p4_neutrino - primarylepton;
  return (lvq + dummyp).M();

  // return std::sqrt(m_p2 + 2 * m_p * (get_enu() - e_mu) - getQ2());
}

TLorentzVector event::get_lvq() const {
  auto p4_neutrino = in_particles.find(14)->second;
  return p4_neutrino - primarylepton;
}