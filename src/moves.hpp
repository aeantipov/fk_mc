#ifndef __FK_MC_MOVES_HPP_
#define __FK_MC_MOVES_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {
 /*
  * This is the insertion move
  */
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
        double ground_energy_old = evals_old(0);
        evals_old-=ground_energy_old;

        auto evals_new = new_config.get_spectrum();
        double ground_energy_new = evals_new(0);
        evals_new-=ground_energy_new;

        auto F = triqs::arrays::map(std::function<double(double)>( [this](double E){return 1.0+exp(-beta*E);} ));
        evals_old = F(evals_old);
        evals_new = F(evals_new);

        mc_weight_type out = std::accumulate(evals_new.begin(), evals_new.end(), 1.0, std::multiplies<double>()); 
        out /= std::accumulate(evals_old.begin(), evals_old.end(), 1.0, std::multiplies<double>());
        out *= exp(-beta*(ground_energy_new - ground_energy_old));
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
