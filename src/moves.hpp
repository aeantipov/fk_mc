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
    static double __calc_weight_ratio(double beta, const real_array_t &evals_old, real_array_t &evals_new);
     config_t& config;
    config_t new_config;
     
    triqs::mc_tools::random_generator &RND;

    move_flip<config_t>(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
    beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt(){
        new_config = config;
        size_t m_size = config.f_config.shape()[0];
        size_t from = RND(m_size); while (new_config.f_config(from)==0) from = RND(m_size);
        size_t to = RND(m_size); while (new_config.f_config(to)==1) to = RND(m_size);

        new_config.f_config(from) = 0;
        new_config.f_config(to) = 1;

        auto evals_old = config.cached_spectrum;
        auto evals_new = new_config.get_spectrum();

        return __calc_weight_ratio(beta, evals_new, evals_old);
    }

    mc_weight_type accept() {
        MY_DEBUG(config.f_config << "->" << new_config.f_config);
        config = new_config; return 1.0; 
    }

    void reject() {}
 };

 //************************************************************************************

template <class lattice>
struct move_randomize : move_flip<lattice> {
    using typename move_flip<lattice>::config_t;
    using typename move_flip<lattice>::mc_weight_type;
    using move_flip<lattice>::beta;
    using move_flip<lattice>::config;
    using move_flip<lattice>::new_config;
    using move_flip<lattice>::RND;
    using move_flip<lattice>::__calc_weight_ratio;
    move_randomize<config_t>(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip<config_t>::move_flip(beta, current_config, RND_) {}

    mc_weight_type attempt(){
        new_config = config;
        new_config.randomize_f(RND, config.get_nf());
        auto evals_old = config.cached_spectrum;
        auto evals_new = new_config.get_spectrum();
        return __calc_weight_ratio(beta, evals_new, evals_old);
    }

};

 //************************************************************************************
template <class lattice>
inline double move_flip<lattice>::__calc_weight_ratio(double beta, const real_array_t &evals_old, real_array_t &evals_new)
{
    double e_min = std::min(evals_new(0),evals_old(0)); 
    double exp_emin = exp(beta*e_min);

    // calculate ( \prod exp(beta*E_0) + \exp(-beta*(E-E_0). The factor exp(beta*E_0) is multiplied on both sides of a fraction P_{n+1}/P_{n}.
    real_array_t evals_rate(evals_new.shape()[0]);// = evals_new / evals_old;
    triqs::clef::placeholder<0> i_;
    evals_rate(i_) << (exp_emin+exp(-beta*(evals_new(i_)-e_min))) / (exp_emin+exp(-beta*(evals_old(i_)-e_min)));
    //MY_DEBUG(evals_rate);

    return __prod(evals_rate);
};



}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
