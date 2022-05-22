#pragma once

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <map>
#include <vector>
template <typename T> class eq_range {
private:
    T pair;

public:
    template <typename T_> eq_range(T_ &&p) : pair(p){}
    decltype(auto) begin() { return pair.first; };
    decltype(auto) end() { return pair.second; };
};
template <typename T_> eq_range(T_ &&) -> eq_range<std::remove_cvref_t<T_>>;
class event {
private:
    std::unordered_multimap<int, TLorentzVector> in_particles{};
    std::unordered_multimap<int, TLorentzVector> out_particles{};

public:
    event() {}
    ~event();
    double getQ2() const;
    double getW() const;
    void add_particle_in(int id, const TLorentzVector &p4);
    void add_particle_out(int id, const TLorentzVector &p4);
    bool is_good_event() const;
    std::pair<double, double> get_q2_w() const;
    double get_enu() const;
    bool TKI_phase_cut() const;
    TLorentzVector get_leading_proton() const;
    auto get_particle_out(int pdgid) { return eq_range{out_particles.equal_range(pdgid)}; }
    auto get_particle_in(int pdgid) { return eq_range{in_particles.equal_range(pdgid)}; }
    size_t count_particle_out (int pdgid) const noexcept;
};
