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

};


template <typename lattice_t>
measure_stiffness<lattice_t>::measure_stiffness(configuration_t& in, const lattice_t& lattice, std::vector<double>& stiffness_vals): 
    config_(in), 
    lattice_(lattice),
    stiffness_vals_(stiffness_vals)
{
}

template <typename lattice_t>
void measure_stiffness<lattice_t>::accumulate(double sign)
{
    config_.calc_ed(true);

    const auto& evecs = config_.ed_data().cached_evecs;
    const auto& evals = config_.ed_data().cached_spectrum;
    const auto& dims = lattice_.dims;

    assert(dims.size() >= 2); // no stiffness in 
    
    auto pos0(dims);
    pos0.fill(0);
    auto pos1(pos0);
    pos0[0] = 1;

    double Lx = dims[0];

    double T = 0, V = 0;

    //Eigen::MatrixXd Vvals(evals.size(),evals.size());
    //Vvals.setZero();

    // loop over eigenvectors
    for (size_t i=0; i<evals.size(); i++) {
        const auto &alpha = evecs.col(i);
        double Te = 0;
        double Ve = 0;
        for (size_t j=0; j<evals.size(); j++) {
            if (std::abs(evals[i] - evals[j]) < 1e-12 && i!=j) continue; 
            const auto &alpha_p = evecs.col(j);

            double T_val = 0;
            std::complex<double> curr_val = 0;
            // loop over y
            for (int y=0; y<dims[1]; y++) { 
                pos0[1] = y; // pos0 = <r|
                pos1[1] = y; // pos1 = |r-x>
                // index0 - r
                auto index0 = lattice_.pos_to_index(pos0);
                // index1 - r-x
                auto index1 = lattice_.pos_to_index(pos1);
                //FKDEBUG(index0 << " " << index1);
                T_val += -t_*std::real(std::conj(alpha[index0]) * alpha_p[index1] + std::conj(alpha[index1]) * alpha_p[index0]);
                curr_val += t_*(std::conj(alpha[index0])*alpha_p[index1] - std::conj(alpha[index1])*alpha_p[index0]);
                };

            double V_val = (i!=j)?2.0*0.5*std::abs(curr_val * curr_val) / (evals[i] - evals[j]):0.0;
            if (i==j) Te = 0.5*T_val;
            Ve+=V_val; 
            //Vvals(i,j)=Ve;
            };
        if (DEBUG_LEVEL > 1) FKDEBUG("-T = " << -Te << " V = " << Ve << " weight = " << 1./(1.0+config_.ed_data().cached_exp(i)) );
        //FKDEBUG(Vvals);
        auto T_val = 1.0/(1.0+config_.ed_data().cached_exp(i)) * Te;
        auto V_val = 1.0/(1.0+config_.ed_data().cached_exp(i)) * Ve;
        if (DEBUG_LEVEL > 1) FKDEBUG("T_val : " << T_val << " V_val : " << V_val);
        T+=T_val;
        V+=V_val;
        }
    double stiffness = V-T;
    if (DEBUG_LEVEL > 0) FKDEBUG("T = " << T << "V = " << V << " stiffness = " << stiffness);
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

