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
     
    triqs::mc_tools::random_generator &RND;

    move_flip(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        beta(beta), config(current_config), new_config(current_config), RND(RND_) {}

    mc_weight_type attempt(){
        if (config.get_nf() == 0 || config.get_nf() == config.lattice.get_msize()) return 0; // this move won't work when the configuration is completely full or empty
        new_config = config;
        size_t m_size = config.lattice.get_msize();
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

struct move_randomize : move_flip {
    using typename move_flip::mc_weight_type;
    using move_flip::beta;
    using move_flip::config;
    using move_flip::new_config;
    using move_flip::RND;
    using move_flip::__calc_weight_ratio;
    move_randomize(double beta, configuration_t& current_config, triqs::mc_tools::random_generator &RND_): 
        move_flip::move_flip(beta, current_config, RND_) {}

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

    mc_weight_type attempt(){
        new_config = config;
        size_t m_size = config.lattice.get_msize();
        size_t to = RND(m_size);
        new_config.f_config(to) = 1 - config.f_config(to);

        new_config.calc_hamiltonian();
        new_config.calc_spectrum();//configuration_t::calc_eval::arpack);
        auto ratio = __calc_weight_ratio(config, new_config);
        //MY_DEBUG("Exp weight: " << t1);
        auto out = (new_config.f_config(to)?ratio*exp_beta_mu_f:ratio/exp_beta_mu_f);
        //MY_DEBUG("weight: " << out);
        return out;
    }
};



}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
