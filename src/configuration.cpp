#include "fk_mc/configuration.hpp"

#include "../eigen/ArpackSupport"

namespace fk {

bool config_params::operator== ( const config_params& rhs) const
{
    double tol = std::numeric_limits<double>::epsilon(); 
    return (std::abs(beta - rhs.beta) < tol && std::abs(U - rhs.U) < tol && std::abs(mu_c - rhs.mu_c) < tol && std::abs(mu_f - rhs.mu_f) < tol);
}

configuration_t::configuration_t(
    const abstract_lattice &lattice, double beta, double U, double mu_c, double mu_f, std::vector<double> W):
    lattice_(lattice),
    f_config_(lattice_.msize()),
    params_(config_params({beta, U, mu_c, mu_f, W})),
    hamilt_(lattice_.hopping_m().rows(), lattice_.hopping_m().cols())
{ 
    f_config_.setZero(); 
}

/*
void configuration_t::swap(configuration_t &rhs)
{
    f_config_.swap(rhs.f_config_); 
    cached_spectrum.swap(rhs.cached_spectrum);
    cached_exp.swap(rhs.cached_exp);
    hamilt_.swap(rhs.hamilt_);
}*/

configuration_t& configuration_t::operator=(const configuration_t& rhs) 
{
    f_config_ = rhs.f_config_; 
    hamilt_ = rhs.hamilt_;
    ed_data_ = rhs.ed_data_;
    cheb_data_ = rhs.cheb_data_;
    if (!(params_ == rhs.params_)) throw (std::logic_error("Mismatched parameters in config assignment"));
    return *this;
};

size_t configuration_t::get_nf() const
{
    return std::accumulate(f_config_.data(), f_config_.data()+ lattice_.msize(), 0);
}

void configuration_t::randomize_f(random_generator &rnd, size_t nf){
    std::uniform_int_distribution<> distr(0, lattice_.msize() - 1);
    if (!nf) nf = distr(rnd);//(lattice_.msize());
    f_config_.setZero();
    for (size_t i=0; i<nf; ++i) {  
        //size_t ind = rnd(lattice_.msize());
        size_t ind = distr(rnd); //rnd(lattice_.msize());
        //while (f_config_(ind)==1) ind = rnd(lattice_.msize());
        while (f_config_(ind)==1) ind = distr(rnd);//(lattice_.msize());
        f_config_(ind) = 1; 
    };
}


double configuration_t::calc_ff_energy() const 
{
    // 1D - easy to add f-f interactions
    if (this->lattice_.ndim() != 1) return 0;
    double e = 0;
    size_t V = lattice_.msize();
    for (int i = 0; i < V; ++i) {  
        if (!f_config_(i)) continue;
        for (int l = 0; l < params_.W.size(); ++l) { 
            int left = (i - l + V)%V;
            int right = (i + l)%V;
            double el = params_.W[l] * f_config_(left);
            double er = params_.W[l] * f_config_(right); 
            e+=el;
            e+=er * (l>0);
            }
        }
    return e;
}

const typename configuration_t::sparse_m& configuration_t::calc_hamiltonian()
{
    reset_cache();
    hamilt_.reserve(lattice_.hopping_m().nonZeros() + lattice_.msize());
    hamilt_ = lattice_.hopping_m();
    for (size_t i=0; i< lattice_.msize(); ++i) hamilt_.coeffRef(i,i)+= -params_.mu_c + params_.U*f_config_(i); // unoptimized
    return hamilt_;
}

void configuration_t::calc_chebyshev( const chebyshev::chebyshev_eval& cheb)
{
    if (int(cheb_data_.status) >= int(chebyshev_cache::logz)) return;
    sparse_m x = hamilt_; 
    size_t msize = lattice_.msize();
    double e_min = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(hamilt_,1,"SA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    double e_max = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(hamilt_,1,"LA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    double a = (e_max - e_min)/2.;
    double b = (e_max + e_min)/2.; 
    double beta = params_.beta;
    cheb_data_.e_min = e_min;
    cheb_data_.e_max = e_max;
    cheb_data_.a = a;
    cheb_data_.b = b;

    for (size_t i=0; i<msize; ++i) x.coeffRef(i,i)+= -b; // unoptimized
    x/=a;


    size_t cheb_size = cheb.cheb_size();
    assert(cheb_size%2 == 0);

    cheb_data_.moments.resize(cheb_size);
    sparse_m cm0(msize,msize);
    cm0.reserve(msize*1.5);
    for (size_t i=0; i<msize; ++i) cm0.insert(i,i)= 1.0; // unoptimized


//    cm0.makeCompressed();
    std::vector<bool> is_set(cheb_size,false);

    cheb_data_.moments[0] = 1.0;
    is_set[0] = true;
    sparse_m cm1 = (x*cm0).pruned(1.0);
    cheb_data_.moments[1] = cm1.diagonal().sum()/msize;
    is_set[1] = true;
    FKDEBUG(cm1.nonZeros() << " nonzero elems [" << msize*msize << "] = " << (double(cm1.nonZeros())/msize/msize));

    sparse_m cm_tmp;

    int m=2;
    bool still_sparse = true;
    // first try to get as much from the sparse matrices as we can - then move on to dense matrices
    for (; m<=cheb_size/2 && still_sparse; ++m) {
            cm_tmp = (x*2.*cm1).pruned(1.0) - cm0; cm0.swap(cm1); cm1.swap(cm_tmp);
            if (!is_set[m]) { 
                cheb_data_.moments[m] = cm1.diagonal().sum()/msize;
                is_set[m] = true;
                FKDEBUG("moment [" << m << "] = " << cheb_data_.moments[m]);
                };
            //FKDEBUG(cm1.nonZeros() << " nonzero elems [" << msize*msize << "] = " << (double(cm1.nonZeros())/msize/msize));
//            std::cout << cm1.nonZeros() << " nonzero elems [" << msize*msize << "] = " << (double(cm1.nonZeros())/msize/msize) << std::endl;
            // add moments for 2*m using relations for Chebyshev polynomials
            int k = 2*(m)-1;
            if (k < cheb_size && k>=cheb_size/2) { 
                double moment_k = (sparse_m(cm0 * cm1).diagonal().sum()*2. - x.diagonal().sum())/msize;
                is_set[k] = true;
                cheb_data_.moments[k] = moment_k;
                FKDEBUG(m << " + moment [" << k << "] = " << cheb_data_.moments[k]);

                if (k!=cheb_size-1) { 
                    ++k;
                    moment_k = (sparse_m(cm1 * cm1).diagonal().sum()/msize*2. - 1.0);
                    cheb_data_.moments[k] = moment_k;
                    is_set[k] = true;
                    FKDEBUG(m << " + moment [" << k << "] = " << cheb_data_.moments[k]);
                };
            }
            still_sparse = (double(cm1.nonZeros())/msize/msize < 0.5);
        }
    // now go on with dense
    if (m <=cheb_size/2) {
        FKDEBUG("Evaluated " << m << " moments with sparse matrices");
        dense_m dm0 = cm0;
        dense_m dm1 = cm1;
        dense_m dm_tmp;
        for (; m<=cheb_size/2; m++) {
                dm_tmp = (x*2*dm1) - dm0; dm0.swap(dm1); dm1.swap(dm_tmp); 
                if (!is_set[m]) { 
                    cheb_data_.moments[m] = dm1.diagonal().sum()/msize;
                    is_set[m] = true;
                    FKDEBUG("moment [" << m << "] = " << cheb_data_.moments[m]);
                    };

                int k = 2*(m)-1;
                if (k < cheb_size && k>=cheb_size/2) { 
                    double moment_k = ((dm0 * dm1).diagonal().sum()*2. - x.diagonal().sum())/msize;
                    is_set[k] = true;
                    cheb_data_.moments[k] = moment_k;
                    FKDEBUG(m << " + moment [" << k << "] = " << cheb_data_.moments[k]);

                    if (k!=cheb_size-1) { 
                        ++k;
                        moment_k = ((dm1 * dm1).diagonal().sum()/msize*2. - 1.0);
                        cheb_data_.moments[k] = moment_k;
                        is_set[k] = true;
                        FKDEBUG(m << " + moment [" << k << "] = " << cheb_data_.moments[k]);
                    };
                }
            }
        }

    assert(m==cheb_size/2+1);

    std::function<double(double)> logz_f = [a,b,beta,msize](double w){return msize*log(1. + exp(-beta*(a*w+b)));}; 
    double s = cheb.moment_f(logz_f, 0);
    for (m=1; m<cheb_size; m++) s+=2.*cheb.moment_f(logz_f, m)*cheb_data_.moments[m];

    cheb_data_.logZ = s;
    cheb_data_.x.swap(x);
    cheb_data_.status = chebyshev_cache::logz;
}

void configuration_t::calc_ed(bool calc_evecs)
{
    if ( (ed_data_.status == ed_cache::spectrum && !calc_evecs) || (ed_data_.status == ed_cache::full && calc_evecs)) return;

    dense_m h(hamilt_);
    Eigen::SelfAdjointEigenSolver<dense_m> s(h,(calc_evecs?Eigen::ComputeEigenvectors:Eigen::EigenvaluesOnly));
    ed_data_.cached_spectrum = s.eigenvalues();
    ed_data_.status = ed_cache::spectrum;
    if (calc_evecs) {
        ed_data_.cached_evecs = s.eigenvectors();
        ed_data_.status = ed_cache::full;
        };
    //auto s2 = cached_spectrum;
    //std::sort (cached_spectrum.data(), cached_spectrum.data()+cached_spectrum.size());  
    //FKDEBUG((Eigen::VectorXd(cached_spectrum - s2)).squaredNorm());

    const auto& cached_spectrum = ed_data_.cached_spectrum;

    double beta = params_.beta;

    double e0 = cached_spectrum[0];
    double logw0 = beta * e0;
    double weight0 = exp(logw0);

    ed_data_.cached_exp.resize(cached_spectrum.size());
    ed_data_.cached_fermi.resize(cached_spectrum.size());
    double logz = 0.0;
    for (size_t i=0; i<cached_spectrum.size(); ++i) { 
        double e = cached_spectrum[i]; 
        double w = exp(-beta*(e-e0)); 
        ed_data_.cached_exp[i] = exp(beta * e);
        ed_data_.cached_fermi[i] = 1.0 / (1.0 + ed_data_.cached_exp[i]);
        assert(!std::isnan(ed_data_.cached_exp[i]));
        logz += std::log(weight0 + w) - logw0;
        };

    ed_data_.logZ = logz; 

}



} // end of namespace fk
