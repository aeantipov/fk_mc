#pragma once 

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>
#include <Eigen/Dense>

namespace fk {

/** A stiffness measurement. 
 * In order to obtain the stiffness for a given configuration, we need it's eigenvectors. 
 * The stiffness is obtained by obtaining kinetic and current operators.
 * We first define the matrices of operators and then obtain averages by matrix-matrix multiplication
 */
template <typename lattice_t>
struct measure_stiffness {
    typedef typename configuration_t::real_array_t  real_array_t;


    measure_stiffness(configuration_t& in, const lattice_t& lattice, 
                  std::vector<double>& stiffness_vals); 
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

    // lattice
    const lattice_t& lattice_;
    // configuration
    configuration_t& config_;

    std::vector<double>& stiffness_vals_;
    double t_ = 1.; // fix for nontrivial lattices
    static constexpr int DEBUG_LEVEL = 2;
    typedef Eigen::MatrixXcd matrix_type;
protected:
    matrix_type Tm_, Jm_;

};


template <typename lattice_t>
measure_stiffness<lattice_t>::measure_stiffness(configuration_t& in, const lattice_t& lattice, std::vector<double>& stiffness_vals): 
    config_(in), 
    lattice_(lattice),
    stiffness_vals_(stiffness_vals),
    Tm_(matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize())),
    Jm_(matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize()))
{
    auto dims = lattice_.dims;
    size_t Lx = dims[0];
    auto pos0(dims),pos_left(dims), pos_right(dims);
    for (size_t i=0; i<lattice_.get_msize();i++) { 
        pos0 = lattice_.index_to_pos(i);
        pos_left = pos0; pos_right = pos0;
        pos_left[0] = ((pos0[0] - 1)+Lx)%Lx; // pos_left = |r-x>
        pos_right[0] = ((pos0[0] + 1)+Lx)%Lx; // pos_right = |r+x>
        auto index0 = lattice_.pos_to_index(pos0);
        assert(lattice_.pos_to_index(pos0) == i);
        auto index_left = lattice_.pos_to_index(pos_left);
        auto index_right = lattice_.pos_to_index(pos_right);
        //FKDEBUG(pos_left << " <- " << pos0 << " -> " << pos_right); 
        //FKDEBUG(index_left << " <- " << index0 << " -> " << index_right);
        Tm_(i, index_left) = -t_;
        Tm_(i, index_right) = -t_;

        Jm_(i, index_left) = -I;
        Jm_(i, index_right) = I;
    }

    FKDEBUG("T = " << Tm_, DEBUG_LEVEL, 3);
    FKDEBUG("J = " << Jm_, DEBUG_LEVEL, 3);
}

template <typename lattice_t>
void measure_stiffness<lattice_t>::accumulate(double sign)
{
    config_.calc_ed(true);

    const auto& evecs = config_.ed_data().cached_evecs;
    const auto& evals = config_.ed_data().cached_spectrum;
    const auto& dims = lattice_.dims;
    const auto& fermi = config_.ed_data().cached_fermi;
    int m_size = evals.size(); 
    double Volume = m_size;
    assert(m_size == lattice_.get_msize());

    assert(dims.size() >= 2); // no stiffness in 

    //auto fermi[](double e){return 1.0/(1.0 + e); }

    //size_t Lx = dims[0];
    
    double T = 0, V = 0;

    for (size_t i=0; i<evals.size(); i++) { 
        double T_val = (-M_PI) * (evecs.col(i).transpose()*Tm_*evecs.col(i) * fermi(i)).real()(0,0);
        T+=T_val;
        for (size_t j=0; j<evals.size(); j++) { 
            bool neq = std::abs(evals[i] - evals[j]) > 1e-12;
            double V_val = neq ? M_PI * (fermi(i) - fermi(j)) / (evals(i) - evals(j)) * (evecs.col(j).transpose() * Jm_ * evecs.col(i)).squaredNorm() : 0.0;
            V+=V_val;
            }
        }

    double stiffness = (V + T)/Volume;
    FKDEBUG("T = " << T << "; V = " << V << "; stiffness = " << stiffness, DEBUG_LEVEL, 2);
    stiffness_vals_.push_back(stiffness);
}

template <typename lattice_t>
void measure_stiffness<lattice_t>::collect_results(boost::mpi::communicator const &c)
{
    std::vector<double> stiffness_vals(stiffness_vals_.size()*c.size());
    boost::mpi::gather(c, stiffness_vals_.data(), stiffness_vals_.size(), stiffness_vals, 0);
    stiffness_vals_.swap(stiffness_vals);
}

} // end of namespace fk

