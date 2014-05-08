#include "moves.hpp"

namespace fk {

typename move_flip::mc_weight_type move_flip::attempt()
{
    if (config.get_nf() == 0 || config.get_nf() == config.lattice_.get_msize()) return 0; // this move won't work when the configuration is completely full or empty
    new_config = config;
    size_t m_size = config.lattice_.get_msize();
    size_t from = RND(m_size); while (new_config.f_config_(from)==0) from = RND(m_size);
    size_t to = RND(m_size); while (new_config.f_config_(to)==1) to = RND(m_size);
    config.calc_ed(false);

    new_config.f_config_(from) = 0;
    new_config.f_config_(to) = 1;

    new_config.calc_hamiltonian();
    new_config.calc_ed(false);//calc_eigenvectors_);
    auto ratio = std::exp(new_config.ed_data_.logZ - config.ed_data_.logZ );
    return ratio;
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
    config.calc_ed(false);
    new_config.calc_ed(false);
    auto log_ratio = new_config.ed_data_.logZ - config.ed_data_.logZ;
    //FKDEBUG(log_ratio);
    if (beta*config.params_.mu_f*(new_config.get_nf()-config.get_nf()) > 2.7182818 - log_ratio) { return 1;}
    else if (beta*config.params_.mu_f*(new_config.get_nf()-config.get_nf()) + log_ratio < 0) {return 0;}
    else return std::exp(log_ratio)*exp(beta*config.params_.mu_f*(new_config.get_nf()-config.get_nf())); 
}

// move_addremove
typename move_addremove::mc_weight_type move_addremove::attempt()
{
    new_config = config;
    size_t m_size = config.lattice_.get_msize();
    size_t to = RND(m_size);
    new_config.f_config_(to) = 1 - config.f_config_(to);

    config.calc_ed(false);
    new_config.calc_hamiltonian();
    new_config.calc_ed(false);//calc_eigenvectors_);//configuration_t::calc_eval::arpack);
    auto ratio = std::exp(new_config.ed_data_.logZ - config.ed_data_.logZ );
    auto out = (new_config.f_config_(to)?ratio*exp_beta_mu_f:ratio/exp_beta_mu_f);
    return out;
}



} // end of namespace fk
