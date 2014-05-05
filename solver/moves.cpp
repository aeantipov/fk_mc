#include "moves.hpp"

namespace fk {

double calc_weight_ratio(const configuration_t &old_config, const configuration_t &new_config);

typename move_flip::mc_weight_type move_flip::attempt()
{
    if (config.get_nf() == 0 || config.get_nf() == config.lattice_.get_msize()) return 0; // this move won't work when the configuration is completely full or empty
    new_config = config;
    size_t m_size = config.lattice_.get_msize();
    size_t from = RND(m_size); while (new_config.f_config_(from)==0) from = RND(m_size);
    size_t to = RND(m_size); while (new_config.f_config_(to)==1) to = RND(m_size);

    new_config.f_config_(from) = 0;
    new_config.f_config_(to) = 1;

    new_config.calc_hamiltonian();
    new_config.calc_ed(false);//calc_eigenvectors_);
    return calc_weight_ratio(config, new_config);
}

typename move_flip::mc_weight_type move_flip::accept() 
{

    config = new_config; 
    return 1.0; 
}

void move_flip::reject() 
{
}

// move_randomize
typename move_randomize::mc_weight_type move_randomize::attempt()
{
    new_config = config;
    //new_config.randomize_f(RND, config.get_nf());
    new_config.randomize_f(RND);
    new_config.calc_hamiltonian();
    new_config.calc_ed(false);
    auto ratio = calc_weight_ratio(config, new_config);
    if (beta*config.params_.mu_f*(new_config.get_nf()-config.get_nf()) > 2.7182818 - log(ratio)) { return 1;}
    else return ratio*exp(beta*config.params_.mu_f*(new_config.get_nf()-config.get_nf())); 
}

// move_addremove
typename move_addremove::mc_weight_type move_addremove::attempt()
{
    new_config = config;
    size_t m_size = config.lattice_.get_msize();
    size_t to = RND(m_size);
    new_config.f_config_(to) = 1 - config.f_config_(to);

    new_config.calc_hamiltonian();
    new_config.calc_ed(false);//calc_eigenvectors_);//configuration_t::calc_eval::arpack);
    auto ratio = calc_weight_ratio(config, new_config);
    auto out = (new_config.f_config_(to)?ratio*exp_beta_mu_f:ratio/exp_beta_mu_f);
    return out;
}



 //************************************************************************************
double calc_weight_ratio(const configuration_t &old_config, const configuration_t &new_config)
{
    size_t size = std::max(old_config.ed_data_.cached_weights.size(), new_config.ed_data_.cached_weights.size());

    typename move_flip::real_array_t weights_old(size), weights_new(size);
    weights_old.setZero(); weights_new.setZero();
    weights_old.head(old_config.ed_data_.cached_weights.size()) = old_config.ed_data_.cached_weights;
    weights_new.head(new_config.ed_data_.cached_weights.size()) = new_config.ed_data_.cached_weights;

    double beta = new_config.params_.beta;
    double exp_e0_old = exp(beta*old_config.ed_data_.cached_spectrum[0]);
    double exp_e0_new = exp(beta*new_config.ed_data_.cached_spectrum[0]);

    auto evals_rate = (exp_e0_new + weights_new) / (exp_e0_old + weights_old ) * (exp_e0_old / exp_e0_new);
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




} // end of namespace fk
