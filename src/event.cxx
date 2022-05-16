#include <event.h>
#include <functional>
#include <set>
#include <sstream>

template <typename T, template <typename...> typename V> inline T add(const V<T> &vec) {
    T sum{};
    for (const T &p4 : vec) {
        sum += p4;
    }
    return sum;
}

template <typename T, template <typename...> typename V> inline T leading(const V<T> &vec) {
    T max{};
    double max_p{};
    for (const T &p4 : vec) {
        if (p4.P() > max_p) {
            max_p = p4.P();
            max = p4;
        }
    }
    return max;
}

event::~event() {}

double event::getQ2() const {
    const auto &p4_muon = out_particles.at(13);
    const auto &p4_neutrino = in_particles.at(14);
    auto q = add(p4_muon) - add(p4_neutrino);
    return -q.M2();
}

double event::getW() const {
    // const auto p_final_state = add(in_particles.at(init_nucleon_pdgid)) - add(out_particles.at(13)) + add(in_particles.at(14));
    TLorentzVector p_final_state{};
    for (const auto &[pdgid, p4v] : out_particles) {
        if (pdgid != 13) {
            p_final_state += add(p4v);
        }
    }
    return p_final_state.M();
}

void event::add_particle_in(int id, const TLorentzVector &p4) {
    if (id > 1000) {
        init_nucleon_pdgid = id;
        // init_nucleon_charge = pdgdb.GetParticle(id)->Charge();
    }
    in_particles[id].push_back(p4);
}

void event::add_particle_out(int id, const TLorentzVector &p4) { out_particles[id].push_back(p4); }

// const std::string &event::get_event_info() {
//     if (!channelname.empty()) {
//         return channelname;
//     }
//     std::stringstream ss;
//     // for (const auto &[id, p4] : in_particles)
//     // {
//     //     ss << pdgdb.GetParticle(id)->GetName();
//     // }
//     // ss << ">";
//     for (const auto &[id, p4] : out_particles) {
//         if (out_particles.at(id).size() > 1) {
//             ss << out_particles.at(id).size();
//         }
//         ss << pdgdb.GetParticle(id)->GetName();
//     }
//     channelname = ss.str();
//     return channelname;
// }

std::size_t event::get_count(int id) const {
    if (out_particles.find(id) == out_particles.end()) {
        return 0;
    }
    return out_particles.at(id).size();
}

std::size_t event::get_pi0_count() const { return get_count(111); }

std::size_t event::get_pip_count() const { return get_count(211); }

std::size_t event::get_pim_count() const { return get_count(-211); }

std::size_t event::get_pions() const { return get_pip_count() + get_pim_count() + get_pi0_count(); }

std::pair<double, double> event::get_q2_w() const { return std::make_pair(getQ2(), getW()); }

bool event::is_good_event() const {
    if (get_pions() > 2 || get_pions() < 1) // 1 or 2 pi
    {
        return false;
    }
    for (const auto &[id, p4] : out_particles) // only allowed particles
    {
        if (std::set{2212, 2112, 13, 211, 111}.count(abs(id)) == 0) // not p, n, mu, pi0, pi+, pi-
        {
            return false;
        }
    }
    if (out_particles.at(13).size() != 1) // one and only one muon
    {
        return false;
    }

    return true;
}

double event::get_enu() const {
    const auto &p4_neutrino = in_particles.at(14);
    assert(p4_neutrino.size() == 1);
    return p4_neutrino[0].E();
}

double event::get_p_mu() const {
    const auto &p4_muon = out_particles.at(13);
    assert(p4_muon.size() == 1);
    return p4_muon[0].P();
}

double event::get_pt_mu() const {
    const auto &p4_muon = out_particles.at(13);
    assert(p4_muon.size() == 1);
    return p4_muon[0].Pt();
}

double event::get_pl_mu() const {
    const auto &p4_muon = out_particles.at(13);
    assert(p4_muon.size() == 1);
    return p4_muon[0].Pz();
}

double event::get_angle_mu() const {
    const auto &p4_muon = out_particles.at(13);
    assert(p4_muon.size() == 1);
    return p4_muon[0].Theta();
}

bool event::TKI_phase_cut() const {
    if (out_particles.count(13) == 0 || out_particles.count(2212) == 0) {
        return false;
    }
    const auto leading_p4_muon = leading(out_particles.at(13));
    const auto p_mu = leading_p4_muon.P();
    const auto angle_mu = leading_p4_muon.Theta() * 180 / M_PI;
    const auto leading_p4_proton = leading(out_particles.at(2212));
    const auto p_proton = leading_p4_proton.P();
    // const auto angle_proton = leading_p4_proton.Theta() * 180 / M_PI;
    return 1.5 < p_mu && p_mu < 20. && angle_mu < 25. && 0.45 < p_proton;
}

TLorentzVector event::get_leading_proton() const { return leading(out_particles.at(2212)); }

std::vector<TLorentzVector> &event::get_particle(int pdgid){
    return out_particles[pdgid];
}