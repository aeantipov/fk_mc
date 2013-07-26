#ifndef __FK_MC_MOVES_HPP_
#define __FK_MC_MOVES_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {

// flip move
template <class lattice>
struct move_flip {
    typedef configuration<lattice> config_t;

    typedef double mc_weight_type;

    double beta;

    config_t& config;
    config_t new_config;
     
    triqs::mc_tools::random_generator &RND;

    // constructor
    move_flip<config_t>(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
    beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    // trial move
    mc_weight_type attempt(){
        new_config = config;
        size_t m_size = config.f_config.shape()[0];
        size_t from = RND(m_size); while (new_config.f_config(from)==0) from = RND(m_size);
        size_t to = RND(m_size); while (new_config.f_config(to)==1) to = RND(m_size);

        new_config.f_config(from) = 0;
        new_config.f_config(to) = 1;

        auto evals_old = config.cached_spectrum;
        auto evals_new = new_config.get_spectrum();

        double e_min = std::min(evals_new(0),evals_old(0)); 
        double exp_emin = exp(beta*e_min);

        // calculate ( \prod exp(beta*E_0) + \exp(-beta*(E-E_0). The factor exp(beta*E_0) is multiplied on both sides of a fraction P_{n+1}/P_{n}.
        auto dm_function = triqs::arrays::map(std::function<double(double)>( [this,exp_emin,e_min](double E){return exp_emin+exp(-beta*(E-e_min));} ));
        evals_old = dm_function(evals_old);
        evals_new = dm_function(evals_new);

        mc_weight_type out = __prod(evals_new) / __prod(evals_old);
        return out;
    }

    mc_weight_type accept() {
    #ifdef FK_MC_DEBUG
        DEBUG(config.f_config << "->" << new_config.f_config);
    #endif
        config = new_config;
        return 1.0; 
    }

    void reject() {
    }
 };

 //************************************************************************************

}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
