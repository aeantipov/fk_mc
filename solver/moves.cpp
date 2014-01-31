
#include "moves.hpp"

namespace fk {

 //************************************************************************************
double move_flip::__calc_weight_ratio(const configuration_t &old_config, const configuration_t &new_config)
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




} // end of namespace fk
