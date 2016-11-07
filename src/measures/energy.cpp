#include <boost/mpi/collectives.hpp>
#include "fk_mc/measures/energy.hpp"

namespace fk {

void measure_energy::accumulate(double sign) 
{
    config.calc_ed(false);
    const auto& spectrum = config.ed_data_.cached_spectrum;
    const auto& exp_e = config.ed_data().cached_exp;
    _Z++;

    real_array_t e_nf(spectrum.size()), d2e_nf(spectrum.size());
    for (size_t i=0; i<e_nf.size(); ++i) {
        e_nf(i) = spectrum(i) / (1.0+exp_e[i]);
        d2e_nf(i) = spectrum(i)*spectrum(i) / (1.0+0.5*(exp_e[i] + 1./exp_e[i]));
        };
    double e_val_c = e_nf.sum();
    double e_val = e_val_c - double(config.params_.mu_f)*config.get_nf();
    double d2e_val = d2e_nf.sum()/2.0;
    _average_energy += e_val;
    _average_d2energy += d2e_val;
    _energies.push_back(e_val);
    _d2energies.push_back(d2e_val);
    _c_energies.push_back(e_val_c);
}

void measure_energy::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    double sum_E, sum_d2E;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    boost::mpi::reduce(c, _average_energy, sum_E, std::plus<double>(), 0);
    boost::mpi::reduce(c, _average_d2energy, sum_d2E, std::plus<double>(), 0);

    std::vector<double> energies(_energies.size()*c.size());
    boost::mpi::gather(c, _energies.data(), _energies.size(), energies, 0);
    _energies.swap(energies);

    c.barrier();
    std::vector<double> d2energies(_d2energies.size()*c.size());
    boost::mpi::gather(c, _d2energies.data(), _d2energies.size(), d2energies, 0);
    _d2energies.swap(d2energies);

    c.barrier();
    std::vector<double> c_energies(_c_energies.size()*c.size());
    boost::mpi::gather(c, _c_energies.data(), _c_energies.size(), c_energies, 0);
    _c_energies.swap(c_energies);

    if (c.rank() == 0) {
    INFO("Total energy: " << sum_E / sum_Z);
    INFO("Total d2energy: " << sum_d2E / sum_Z);
    //INFO("Total c_energy: " << std::accumulate(_c_energies.begin(), _c_energies.end(), 0.0, std::plus<double>()) / sum_Z);
    }
}



} // end of namespace FK
