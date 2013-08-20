#ifndef __FK_MC_MOVES_HPP_
#define __FK_MC_MOVES_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include "triqs_extra.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {

// flip move
template <class config_t>
struct move_flip {
    typedef double mc_weight_type;
    typedef typename config_t::real_array_t  real_array_t;

    double beta;
    static double __calc_weight_ratio(const config_t &old_config, const config_t &new_config);
    config_t& config;
    config_t new_config;
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt(){
        if (config.get_nf() == 0 || config.get_nf() == config.lattice.m_size) return 0; // this move won't work when the configuration is completely full or empty
        new_config = config;
        size_t m_size = config.lattice.m_size;
        size_t from = RND(m_size); while (new_config.f_config(from)==0) from = RND(m_size);
        size_t to = RND(m_size); while (new_config.f_config(to)==1) to = RND(m_size);

        new_config.f_config(from) = 0;
        new_config.f_config(to) = 1;

        new_config.calc_hamiltonian();
        new_config.calc_spectrum();
        return __calc_weight_ratio(config, new_config);
    }

    mc_weight_type accept() {
        config = new_config; 
        return 1.0; 
    }

    void reject() {}
 };

//************************************************************************************

template <class config_t>
struct move_randomize : move_flip<config_t> {
    using typename move_flip<config_t>::mc_weight_type;
    using move_flip<config_t>::beta;
    using move_flip<config_t>::config;
    using move_flip<config_t>::new_config;
    using move_flip<config_t>::RND;
    using move_flip<config_t>::__calc_weight_ratio;
    move_randomize<config_t>(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip<config_t>::move_flip(beta, current_config, RND_) {}

    mc_weight_type attempt(){
        new_config = config;
        //new_config.randomize_f(RND, config.get_nf());
        new_config.randomize_f(RND);
        new_config.calc_hamiltonian();
        new_config.calc_spectrum();
        auto ratio = __calc_weight_ratio(config, new_config);
        //MY_DEBUG(ratio << "*" << exp(beta*config.mu_f*(new_config.get_nf()-config.get_nf())) << "=" << exp(beta*config.mu_f*(new_config.get_nf()-config.get_nf())));
        if (beta*config.mu_f*(new_config.get_nf()-config.get_nf()) > 2.7182818 - log(ratio)) { return 1;}
        else return ratio*exp(beta*config.mu_f*(new_config.get_nf()-config.get_nf())); 
    }

};

//************************************************************************************

template <class config_t>
struct move_addremove : move_flip<config_t> {
    using typename move_flip<config_t>::mc_weight_type;
    using move_flip<config_t>::beta;
    using move_flip<config_t>::config;
    using move_flip<config_t>::new_config;
    using move_flip<config_t>::RND;
    using move_flip<config_t>::__calc_weight_ratio;
    double exp_beta_mu_f;
    move_addremove<config_t>(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip<config_t>::move_flip(beta, current_config, RND_),exp_beta_mu_f(exp(beta*config.mu_f)) {}

    mc_weight_type attempt(){
        new_config = config;
        size_t m_size = config.lattice.m_size;
        size_t to = RND(m_size);
        new_config.f_config(to) = 1 - config.f_config(to);

        new_config.calc_hamiltonian();
        new_config.calc_spectrum(config_t::calc_eval::arpack);
        auto ratio = __calc_weight_ratio(config, new_config);
        //MY_DEBUG("Exp weight: " << t1);
        auto out = (new_config.f_config(to)?ratio*exp_beta_mu_f:ratio/exp_beta_mu_f);
        //MY_DEBUG("weight: " << out);
        return out;
    }
};


 //************************************************************************************
template <class config_t>
inline double move_flip<config_t>::__calc_weight_ratio(const config_t &old_config, const config_t &new_config)
{
    size_t size = std::max(old_config.cached_weights.size(), new_config.cached_weights.size());

    real_array_t weights_old(size), weights_new(size);
    weights_old.setZero(); weights_new.setZero();
    weights_old.head(old_config.cached_weights.size()) = old_config.cached_weights;
    weights_new.head(new_config.cached_weights.size()) = new_config.cached_weights;

    double exp_e0_old = exp(old_config.beta*old_config.cached_spectrum[0]);
    double exp_e0_new = exp(new_config.beta*new_config.cached_spectrum[0]);

    real_array_t evals_rate = (exp_e0_new + weights_new) / (exp_e0_old + weights_old ) * (exp_e0_old / exp_e0_new);
   /* 
    double exp_e0_old  = std::min(evals_new(0),evals_old(0)); 
    double exp_emin = exp(beta*e_min);

    // calculate ( \prod exp(beta*E_0) + \exp(-beta*(E-E_0). The factor exp(beta*E_0) is multiplied on both sides of a fraction P_{n+1}/P_{n}.
    real_array_t evals_rate(evals_new.size());// = evals_new / evals_old;
    for (size_t i=0; i<evals_rate.size(); ++i) 
        evals_rate(i) = (exp_emin+exp(-beta*(evals_new(i)-e_min))) / (exp_emin+exp(-beta*(evals_old(i)-e_min)));
    */

    return evals_rate.prod();
};



}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
