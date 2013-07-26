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
    config_t * config;
     
    triqs::mc_tools::random_generator &RND;

    // constructor
    move_flip<config_t>(config_t& current_config, triqs::mc_tools::random_generator &RND_): 
    config(&current_config), RND(RND_) {}

    // trial move
    mc_weight_type attempt(){
    /*
        size_t m_size = old_f.shape()[0];
        size_t from = RND(m_size); while (old_f(from)==0) from = RND(m_size);
        size_t to = RND(m_size); while (old_f(to)==1) to = RND(m_size);

        real_array_t new_f(old_f); 
        new_f(from)=0; new_f(to)=1;
        DEBUG(old_f << "->" << new_f);
        DEBUG(l.get_hopping_matrix());

        auto evals_old = triqs::arrays::linalg::eigenvalues(l.get_hopping_matrix());
        DEBUG(evals_old);
        //auto evals_new = triqs::arrays::linalg
    */
    return 1.0;
    }

    mc_weight_type accept() {
    return 1.0; 
    }

    void reject() {
    }
 };

 //************************************************************************************

}

#endif // endif :: ifndef __FK_MC_MOVES_HPP_
