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
    static double __calc_weight_ratio(const configuration_t &old_config, const configuration_t &new_config);
    configuration_t& config;
    configuration_t new_config;
    bool calc_eigenvectors_ = false;
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt();
    mc_weight_type accept();
    void reject();
 };

//************************************************************************************

struct move_randomize : move_flip {
    using typename move_flip::mc_weight_type;
    using move_flip::beta;
    using move_flip::config;
    using move_flip::new_config;
    using move_flip::RND;
    using move_flip::__calc_weight_ratio;
    move_randomize(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_) {}

    mc_weight_type attempt();
};

//************************************************************************************

struct move_addremove : move_flip {
    using typename move_flip::mc_weight_type;
    using move_flip::beta;
    using move_flip::config;
    using move_flip::new_config;
    using move_flip::RND;
    using move_flip::__calc_weight_ratio;
    double exp_beta_mu_f;
    move_addremove(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_),exp_beta_mu_f(exp(beta*config.mu_f)) {}

    mc_weight_type attempt();
};

}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
