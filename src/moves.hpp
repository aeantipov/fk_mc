#ifndef __FK_MC_MOVES_HPP_
#define __FK_MC_MOVES_HPP_

#include "common.hpp"
#include "configuration.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace fk {

// flip move
template <class config_t>
struct move_flip {
    typedef double mc_weight_type;

    double beta;
    static double __calc_weight_ratio(double beta, const real_array_t &evals_old, real_array_t &evals_new);
     config_t& config;
    config_t new_config;
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, config_t& current_config, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt(){
        if (config.get_nf() == 0 || config.get_nf() == config.lattice.m_size) return 0; // this move won't work when the configuration is completely full or empty
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
        #ifdef FK_MC_DEBUG
        //auto x=triqs::arrays::immutable_diagonal_matrix_view<int>(new_config.f_config);
        auto x = reinterpret_array_view(new_config.f_config,new_config.lattice.dims[0],new_config.lattice.dims[1]);
        MY_DEBUG(config.f_config << "->" << new_config.f_config);
        #endif
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
        auto evals_old = config.cached_spectrum;
        auto evals_new = new_config.get_spectrum();
        auto ratio = __calc_weight_ratio(beta, evals_new, evals_old);
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
        size_t m_size = config.f_config.shape()[0];
        size_t to = RND(m_size);
        new_config.f_config(to) = 1 - config.f_config(to);
        auto evals_old = config.cached_spectrum;
        auto evals_new = new_config.get_spectrum();
        auto t1 = __calc_weight_ratio(beta, evals_new, evals_old);
        //MY_DEBUG("Exp weight: " << t1);
        auto out = (new_config.f_config(to)?t1*exp_beta_mu_f:t1/exp_beta_mu_f);
        //MY_DEBUG("weight: " << out);
        return out;
    }
};


 //************************************************************************************
template <class config_t>
inline double move_flip<config_t>::__calc_weight_ratio(double beta, const real_array_t &evals_new, real_array_t &evals_old)
{
    double e_min = std::min(evals_new(0),evals_old(0)); 
    double exp_emin = exp(beta*e_min);

    // calculate ( \prod exp(beta*E_0) + \exp(-beta*(E-E_0). The factor exp(beta*E_0) is multiplied on both sides of a fraction P_{n+1}/P_{n}.
    real_array_t evals_rate(evals_new.shape()[0]);// = evals_new / evals_old;
    triqs::clef::placeholder<0> i_;
    evals_rate(i_) << (exp_emin+exp(-beta*(evals_new(i_)-e_min))) / (exp_emin+exp(-beta*(evals_old(i_)-e_min)));
    //MY_DEBUG("Evals ratio: " << evals_rate);
//    evals_rate(i_) << (1.0+exp(-beta*(evals_new(i_)))) / (1.0+exp(-beta*(evals_old(i_))));

    return __prod(evals_rate);
};



}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
