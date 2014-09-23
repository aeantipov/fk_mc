#pragma once
#include <triqs/h5.hpp>

#include "fk_mc.hpp"

namespace fk { 

// save data from solver to hdf5 file 
template <typename MC>
void save_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false);

// construct a solver from given hdf5 file. The observables are populated with existing data
template <typename MC>
MC load_data(std::string output_file, triqs::utility::parameters pnew)
{
    boost::mpi::communicator world;
    bool success = true; 

    // first we construct mc from parameters (mostly the ones saved in hdf5 file, except of number of cycles, random seed, etc
    H5::H5File input(output_file.c_str(),H5F_ACC_RDONLY);
    triqs::h5::group top(input);

    triqs::utility::parameters pold, p;
    h5_read(top, "parameters", pold);

    //if (!bool(p["measure_history"])) TRIQS_RUNTIME_ERROR << "No data to load. Need measure_history = 1 for that" << std::end

    success = 
        (!pold.has_key("t") || std::abs(double(pnew["t"]) - double(pold["t"])) < 1e-4) && 
        (!pold.has_key("L") || int(pold["L"]) == int(pnew["L"])) && 
        std::abs(double(pnew["U"]) - double(pold["U"])) < 1e-4 && 
        std::abs(double(pnew["beta"]) - double(pold["beta"])) < 1e-4 && 
        bool(pnew["measure_history"]) == bool(pold["measure_history"]) && 
        bool(pnew["measure_ipr"]) == bool(pold["measure_ipr"]);

    if (!success) { 
            if (!world.rank()) std::cout << "Parameters mismatch" << std::endl << "old: " << pold << std::endl << "new: " << pnew << std::endl; 
            TRIQS_RUNTIME_ERROR << "Parameters mismatch";
            };

    p.update(pnew);
    p.update(pold);
    // update ncycles and max_time but keep old cycle length. No warmups
    p["n_cycles"] = std::max(int(pnew["n_cycles"]) - int(pold["n_cycles"]),0); 
    p["max_time"] = pnew["max_time"];
    p["n_warmup_cycles"] = int(0);
    p["random_seed"] = pnew["random_seed"];
    if (!world.rank()) std::cout << "params : " << p << std::endl;

    int L = p["L"];
    typedef typename MC::lattice_type lattice_t;
    lattice_t lattice(L);

    fk_mc<lattice_t> mc(lattice,p); 
    // now load actual data

    auto h5_mc_data = top.open_group("mc_data");
    if (!world.rank()) std::cout << "Loading mc observables... " << std::flush;
    h5_read(h5_mc_data,"energies", mc.observables.energies);
    h5_read(h5_mc_data,"d2energies", mc.observables.d2energies);
    h5_read(h5_mc_data,"nf0", mc.observables.nf0);
    h5_read(h5_mc_data,"nfpi", mc.observables.nfpi);
    std::vector<double> spectrum;
    h5_read(h5_mc_data,"spectrum", spectrum);
    mc.observables.spectrum.resize(spectrum.size());
    std::copy(spectrum.begin(), spectrum.end(), mc.observables.spectrum.data());

    if (p["measure_history"]) { 
        triqs::arrays::array<double, 2> t_spectrum_history; 
        h5_read(h5_mc_data,"spectrum_history", t_spectrum_history);
        mc.observables.spectrum_history.resize(t_spectrum_history.shape()[0]);
        for (int i=0; i<mc.observables.spectrum_history.size(); i++) { 
            mc.observables.spectrum_history[i].resize(t_spectrum_history.shape()[1]);
            for (int j=0; j< mc.observables.spectrum_history[0].size(); j++)
                mc.observables.spectrum_history[i][j] = t_spectrum_history(i,j);
            }
        
        triqs::arrays::array<double, 2> focc_history;
        h5_read(h5_mc_data,"focc_history", focc_history);
        mc.observables.focc_history.resize(focc_history.shape()[0]);
        for (int i=0; i<mc.observables.focc_history.size(); i++) { 
            mc.observables.focc_history[i].resize(focc_history.shape()[1]);
            for (int j=0; j< mc.observables.focc_history[0].size(); j++)
                focc_history(i,j) =  mc.observables.focc_history[i][j];
            }
        };

    // Inverse participation ratio
    if (p["measure_ipr"] && p["measure_history"]) {
        std::cout << "loading ipr..." << std::flush;
        triqs::arrays::array<double, 2> t_ipr_history;
        h5_read(h5_mc_data,"ipr_history", t_ipr_history);
        mc.observables.ipr_history.resize(t_ipr_history.shape()[0]);
        for (int i=0; i<mc.observables.ipr_history.size(); i++) { 
            mc.observables.ipr_history[i].resize(t_ipr_history.shape()[1]);
            for (int j=0; j< mc.observables.ipr_history[0].size(); j++)
                mc.observables.ipr_history[i][j] = t_ipr_history(i,j); 
            };
        };
    if (!world.rank()) std::cout << "done." << std::endl;

    exit(0);
}



} // end of namespace fk

#include "data_save.hxx"


