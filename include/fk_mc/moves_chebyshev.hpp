#ifndef __FK_MC_MOVES_CHEBYSHEV_HPP_
#define __FK_MC_MOVES_CHEBYSHEV_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include "chebyshev.hpp" 
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {

namespace chebyshev { 
// flip move
struct move_flip {
    typedef double mc_weight_type;
    typedef typename configuration_t::real_array_t  real_array_t;

    double beta;
    configuration_t& config;
    configuration_t new_config;
    const chebyshev_eval& cheb_;
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, configuration_t& current_config, const chebyshev_eval& cheb, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), cheb_(cheb), RND(RND_) {}

    mc_weight_type attempt();
    mc_weight_type accept();
    void reject();
 };

//************************************************************************************

struct move_randomize : move_flip {
    move_randomize(double beta, configuration_t& current_config, const chebyshev_eval& cheb, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, cheb, RND_) {}

    mc_weight_type attempt();
};

//************************************************************************************

struct move_addremove : move_flip {
    double exp_beta_mu_f;
    move_addremove(double beta, configuration_t& current_config, const chebyshev_eval& cheb, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, cheb, RND_),exp_beta_mu_f(exp(beta*config.params_.mu_f)) {}

    mc_weight_type attempt();
};

} // end of namespace chebyshev
}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
