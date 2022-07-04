#pragma once

// #include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <map>
#include <set>
#include <vector>
template <typename T> class eq_range {
private:
    T pair;

public:
    template <typename T_> eq_range(T_ &&p) : pair(p) {}
    decltype(auto) begin() { return pair.first; };
    decltype(auto) end() { return pair.second; };
};
template <typename T_> eq_range(T_ &&) -> eq_range<std::remove_cvref_t<T_>>;
class event {
public:
    enum class channel { QE, RES, DIS, MEC, Other };

private:
    std::unordered_multimap<int, TLorentzVector> in_particles{};
    std::unordered_multimap<int, TLorentzVector> out_particles{};
    // std::set<int> pdg_list_out{};
    std::map<int, size_t> pdg_list_out{};
    std::unordered_multimap<int, TLorentzVector> nofsi_particles{};
    double weight{1};
    channel mode;
    // TDatabasePDG db{};
    std::string channelname{};
    std::string channelname_nonucleon{};

public:
    event() {}
    ~event();
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
    auto get_particle_out(int pdgid) { return eq_range{out_particles.equal_range(pdgid)}; }
    auto get_particle_in(int pdgid) { return eq_range{in_particles.equal_range(pdgid)}; }
    auto get_particle_nofsi(int pdgid) { return eq_range{nofsi_particles.equal_range(pdgid)}; }
    size_t count_particle_out(int pdgid) const noexcept;
    size_t count_particle_nofsi(int pdgid) const noexcept;
    bool TKI_mu_p_cut() const;
    bool TKI_mu_cut() const;
    void set_mode(channel);
    channel get_mode() const;
    void set_weight(double);
    double get_weight() const;
    const std::unordered_multimap<int, TLorentzVector> &get_particle_out() const { return out_particles; };
    const std::unordered_multimap<int, TLorentzVector> &get_particle_in() const { return in_particles; };
    const std::unordered_multimap<int, TLorentzVector> &get_particle_nofsi() const { return nofsi_particles; };
    size_t count_out(int) const noexcept;
    const std::string & get_channelname();
    const std::string & get_channelname_no_nucleon();
};
