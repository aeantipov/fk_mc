#include "fk_mc/lattice/chain.hpp"

namespace fk { 

void chain_lattice::fill(double t, double eta, double delta)
{
    double t1 = t-eta;
    double t2 = t+eta;

    if (m_size_%2 != 0) { FK_ERROR << "Can't have a chain of even elements"; exit(1); };
    
    for (size_t i=0; i<m_size_; i+=2) {
        hopping_m_.insert(i,i) = eta;
        hopping_m_.insert(i+1,i+1) = -eta;
        hopping_m_.insert(i,i+1)   = -1.0*t1; hopping_m_.insert(i+1,i)   = -1.0*t1;
        hopping_m_.insert(i+1,(i+2)%(m_size_)) = -1.0*t2; hopping_m_.insert((i+2)%m_size_,i+1) = -1.0*t2;
        };
}



} // end of namespace fk
