#ifndef __FK_MC_MOVES_HPP_
#define __FK_MC_MOVES_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include "triqs_extra.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {

// flip move
struct move_flip {
    typedef double mc_weight_type;
    typedef typename configuration_t::real_array_t  real_array_t;

    double beta;
    configuration_t& config;
    configuration_t new_config;
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt();
    mc_weight_type accept();
    void reject();
 };

//************************************************************************************

struct move_randomize : move_flip {
    move_randomize(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_) {}

    mc_weight_type attempt();
};

//************************************************************************************

struct move_addremove : move_flip {
    double exp_beta_mu_f;
    move_addremove(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_),exp_beta_mu_f(exp(beta*config.params_.mu_f)) {}

    mc_weight_type attempt();
};
//************************************************************************************
/*
template <typename Lattice>
struct move_cluster : move_flip {
    typedef Lattice lattice_t;
    const lattice_t& lattice_;
    move_cluster(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_), lattice_(static_cast<Lattice>(config_.lattice_)){}
    mc_weight_type attempt();
};
*/


}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
