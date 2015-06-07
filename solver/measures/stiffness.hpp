#pragma once 

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>
#include <Eigen/Dense>

namespace fk {

template <typename T> bool almost_equal(T x, T y, double tol = std::numeric_limits<double>::epsilon()) { return std::abs(x - y) < tol; }

// approximate delta(x)
inline constexpr double lorentzian(double w, double eta) { return eta / M_PI / (w*w + eta*eta); } 

struct resonant_term { 
public:
    resonant_term(double energy, double coeff):energy(energy), coeff(coeff){};
    resonant_term() = default;
    bool operator<(resonant_term const& r) const {return this->energy < r.energy; } 
    double operator()(double w, double eta) const {return lorentzian(w + energy, eta) * coeff; }
    resonant_term& operator+=(resonant_term const&rhs){coeff += rhs.coeff; return *this; }
    resonant_term& operator/=(double x){coeff /= x; return *this; }
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & energy; ar & coeff; }
    friend std::ostream& operator<<(std::ostream& in, resonant_term const& x){in << x.coeff << " \u0394_(w - " << x.energy << ")"; return in; }

    double energy;
    double coeff; 
};



/** A stiffness measurement. 
 * In order to obtain the stiffness for a given configuration, we need it's eigenvectors. 
 * The stiffness is obtained by obtaining kinetic and current operators.
 * We first define the matrices of operators and then obtain averages by matrix-matrix multiplication
 */
template <typename lattice_t>
struct measure_stiffness {
    typedef typename configuration_t::real_array_t  real_array_t;
    typedef Eigen::MatrixXd matrix_type;
    typedef Eigen::SparseMatrix<double> sparse_m;

    measure_stiffness(configuration_t& in, const lattice_t& lattice, std::vector<double>& stiffness_vals, 
                      std::vector<std::vector<double>>& cond_history,
                      std::vector<double> wgrid, double offset); 
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

protected:
    // lattice
    const lattice_t& lattice_;
    // configuration
    configuration_t& config_;

    std::vector<double>& stiffness_vals_;
    std::vector<std::vector<double>>& cond_history_;
    double t_ = 1.; // fix for nontrivial lattices
    static constexpr int DEBUG_LEVEL = 2;

    std::vector<double> wgrid_;
    double offset_;
    sparse_m Tm_, Jm_;
    //matrix_type Tm_, Jm_;
};


template <typename lattice_t>
measure_stiffness<lattice_t>::measure_stiffness(configuration_t& in, const lattice_t& lattice, 
        std::vector<double>& stiffness_vals,
        std::vector<std::vector<double>>& cond_history, 
        std::vector<double> wgrid, double offset
        ): 
    config_(in), 
    lattice_(lattice),
    stiffness_vals_(stiffness_vals),
    cond_history_(cond_history),
    wgrid_(wgrid),
    offset_(offset),
    Tm_(lattice_.get_msize(), lattice_.get_msize()), //matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize())),
    Jm_(lattice_.get_msize(), lattice_.get_msize()) //matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize()))
    //Tm_(matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize())),
    //Jm_(matrix_type::Zero(lattice_.get_msize(), lattice_.get_msize()))
{
    cond_history_.resize(wgrid.size());
    auto dims = lattice_.dims;
    size_t Lx = dims[0];
    auto pos0(dims),pos_left(dims), pos_right(dims);
    std::vector<Eigen::Triplet<double>> T_triplets, J_triplets;
    T_triplets.reserve(2*lattice_.get_msize());
    J_triplets.reserve(2*lattice_.get_msize());
    for (size_t i=0; i<lattice_.get_msize();i++) { 
        pos0 = lattice_.index_to_pos(i);
        pos_left = pos0; pos_right = pos0;
        pos_left[0] = ((pos0[0] - 1)+Lx)%Lx; // pos_left = |r-x>
        pos_right[0] = ((pos0[0] + 1)+Lx)%Lx; // pos_right = |r+x>
        auto index0 = lattice_.pos_to_index(pos0);
        assert(lattice_.pos_to_index(pos0) == i);
        auto index_left = lattice_.pos_to_index(pos_left);
        auto index_right = lattice_.pos_to_index(pos_right);

        Eigen::Triplet<double> t_left(i, index_left, -t_);
        Eigen::Triplet<double> t_right(i, index_right, -t_);
        Eigen::Triplet<double> j_left(i, index_left, -1); // in fact it is -I, but j is only used as |J|^2
        Eigen::Triplet<double> j_right(i, index_right, 1);

        T_triplets.push_back(t_left);
        T_triplets.push_back(t_right);
        J_triplets.push_back(j_left);
        J_triplets.push_back(j_right);

        //FKDEBUG(pos_left << " <- " << pos0 << " -> " << pos_right); 
        //FKDEBUG(index_left << " <- " << index0 << " -> " << index_right);
    //    Tm_(i, index_left) = -t_;
    //    Tm_(i, index_right) = -t_;

    //    Jm_(i, index_left) = -I;
    //    Jm_(i, index_right) = I;
    }

    Tm_.setFromTriplets(T_triplets.begin(), T_triplets.end());
    Jm_.setFromTriplets(J_triplets.begin(), J_triplets.end());

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
    double coeff_tolerance_ = 1e-13;
    double reduce_tolerance_ = 1e-4;
    
    double T = 0, V = 0;

    matrix_type mJ = evecs.transpose() * Jm_ * evecs;
    T = -M_PI*(Eigen::RowVectorXd((evecs.transpose() * Tm_ * evecs).diagonal()) * fermi.matrix())(0,0);
    std::list<resonant_term> cond_terms;

    for (size_t i=0; i<evals.size(); i++) { 
        //double T_val = (-M_PI) * (evecs.col(i).transpose()*Tm_*evecs.col(i) * fermi(i)).real()(0,0);
        //T+=T_val;
        for (size_t j=0; j<i; j++) { 
            bool neq = std::abs(evals[i] - evals[j]) > 1e-12;
            if (neq) { 
                double sigma_v = neq ? M_PI * (fermi(j) - fermi(i)) * mJ(j,i) * mJ(i,j) : 0;
                double V_val = neq ? 2.*sigma_v / (evals(i) - evals(j)) : 0.0;
                if (std::abs(sigma_v) > coeff_tolerance_) { 
                    cond_terms.emplace_back(evals(j) - evals(i), sigma_v);
                    cond_terms.emplace_back(evals(i) - evals(j), -sigma_v);
                    V+=V_val;
                    }
                }
            }
        }

    double stiffness = (V + T)/Volume;


    FKDEBUG("T = " << T << "; V = " << V << "; stiffness = " << stiffness, DEBUG_LEVEL, 2);
    stiffness_vals_.push_back(stiffness);

    for (int wn = 0; wn < wgrid_.size(); wn++) { 
        cond_history_[wn].push_back(std::accumulate(cond_terms.begin(), cond_terms.end(), 0.0, 
                                     [&](double s, resonant_term const& d)->double { 
                                        //FKDEBUG(wgrid_[wn] << " : " << s << " + " << d << " " << offset_ << " : " << d(wgrid_[wn], offset_)); 
                                        return s+d(wgrid_[wn], offset_);
                                        }
                                    )); 
        }
}

template <typename lattice_t>
void measure_stiffness<lattice_t>::collect_results(boost::mpi::communicator const &comm)
{
    using namespace boost::mpi;
    int nmeasures = stiffness_vals_.size()*comm.size();
    std::vector<double> stiffness_vals(nmeasures);
    gather(comm, stiffness_vals_.data(), stiffness_vals_.size(), stiffness_vals, 0);
    stiffness_vals_.swap(stiffness_vals);

    std::vector<double> current_cond;
    assert(cond_history_.size() == wgrid_.size());
    for (size_t i=0; i<cond_history_.size(); ++i) {
        current_cond.resize(cond_history_[i].size()*comm.size());
        boost::mpi::gather(comm, cond_history_[i].data(), cond_history_[i].size(), current_cond, 0);
        cond_history_[i].swap(current_cond);
        };
}

} // end of namespace fk

/*
    // Reduction of terms 
    std::cout << "was : " << cond_terms_.size() << std::endl;
    for(std::list<resonant_term>::iterator it1 = cond_terms_.begin(); it1 != cond_terms_.end();++it1){
        std::cout << *it1 << " " << std::flush; } 
    for(std::list<resonant_term>::iterator it1 = cond_terms_.begin(); it1 != cond_terms_.end();){
        std::list<resonant_term>::iterator it2 = it1;
        for(it2++; it2 != cond_terms_.end();){
            if(std::abs(it1->energy - it2->energy) < reduce_tolerance_){
                *it1 += *it2;
                it2 = cond_terms_.erase(it2);
            }else
                it2++;
        }
        if(std::abs(it1->coeff) < coeff_tolerance_)
            it1 = cond_terms_.erase(it1);
        else
            it1++;
    }
    std::cout << "now : " << cond_terms_.size() << " / " << m_size*m_size << std::endl;
   */ 

