#pragma once

#include <TLorentzVector.h>
#include <map>
#include <vector>

class event {
private:
    std::map<int, std::vector<TLorentzVector>> in_particles{};
    std::map<int, std::vector<TLorentzVector>> out_particles{};
    // TLorentzVector init_nucleon{};
    int init_nucleon_pdgid{};
    int init_nucleon_charge{};
    double wgt{};
    std::string channelname{};

public:
    event();
    ~event();
    double getQ2() const;
    double getW() const;
    void add_particle_in(int id, const TLorentzVector &p4);
    void add_particle_out(int id, const TLorentzVector &p4);
    const std::string &get_event_info();
    std::size_t get_pi0_count() const;
    std::size_t get_pip_count() const;
    std::size_t get_pim_count() const;
    std::size_t get_pions() const;
    std::size_t get_count(int) const;
    bool is_good_event() const;
    std::pair<double, double> get_q2_w() const;
    double get_enu() const;
    double get_p_mu() const;
    double get_pt_mu() const;
    double get_pl_mu() const;
    double get_angle_mu() const;
};
