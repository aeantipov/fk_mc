#pragma once
#include <triqs/h5.hpp>

#include "fk_mc.hpp"

namespace fk { 


// construct a solver from given hdf5 file. The observables are populated with existing data
triqs::utility::parameters load_parameters(std::string output_file, triqs::utility::parameters pnew)
{
    boost::mpi::communicator world;
    bool success = true; 

    // first we construct mc from parameters (mostly the ones saved in hdf5 file, except of number of cycles, random seed, etc
    H5::H5File input(output_file.c_str(),H5F_ACC_RDONLY);
    triqs::h5::group top(input);

    triqs::utility::parameters pold, p;
    h5_read(top, "parameters", pold);

    // Check that the most important parameters are the same in new and old runs
    success = 
        (!pold.has_key("t") || std::abs(double(pnew["t"]) - double(pold["t"])) < 1e-4) && 
        (!pold.has_key("L") || int(pold["L"]) == int(pnew["L"])) && 
        std::abs(double(pnew["U"]) - double(pold["U"])) < 1e-4 && 
        std::abs(double(pnew["beta"]) - double(pold["beta"])) < 1e-4 && 
        bool(pnew["measure_history"]) == bool(pold["measure_history"]) && 
        bool(pnew["measure_ipr"]) == bool(pold["measure_ipr"]) && 
        bool(pnew["cheb_moves"]) == bool(pold["cheb_moves"]) &&
        double(pnew["cheb_prefactor"]) == double(pold["cheb_prefactor"]);

    if (!success) { 
            if (!world.rank()) std::cerr << "Parameters mismatch" << std::endl << "old: " << pold << std::endl << "new: " << pnew << std::endl; 
            TRIQS_RUNTIME_ERROR << "Parameters mismatch";
            };

    // first put old parameters and then override them with new ones
    p.update(pold);
    p.update(pnew);
    // update ncycles and max_time but keep old cycle length
    p["n_cycles"] = std::max(int(pnew["n_cycles"]) + int(pold["n_cycles"]),0); 
    p["length_cycle"] = pold["length_cycle"];
    if (!world.rank()) std::cout << "params : " << p << std::endl;

    return p;
}

observables_t load_observables(std::string output_file, triqs::utility::parameters p)
{
    // now load actual data
    H5::H5File input(output_file.c_str(),H5F_ACC_RDONLY);
    triqs::h5::group top(input);
    observables_t obs;

    auto h5_mc_data = top.open_group("mc_data");
    std::cout << "Loading observables... " << std::flush;
    if (h5_mc_data.has_key("energies")) h5_read(h5_mc_data,"energies", obs.energies);
    if (h5_mc_data.has_key("d2energies")) h5_read(h5_mc_data,"d2energies", obs.d2energies);
    h5_read(h5_mc_data,"nf0", obs.nf0);
    h5_read(h5_mc_data,"nfpi", obs.nfpi);
    if (h5_mc_data.has_key("specteum")) h5_read(h5_mc_data,"spectrum", obs.spectrum);

    if (p["measure_history"]) { 
        std::cout << "spectrum_history... " << std::flush;
        triqs::arrays::array<double, 2> t_spectrum_history; 
        h5_read(h5_mc_data,"spectrum_history", t_spectrum_history);
        obs.spectrum_history.resize(t_spectrum_history.shape()[0]);
        for (int i=0; i<obs.spectrum_history.size(); i++) { 
            obs.spectrum_history[i].resize(t_spectrum_history.shape()[1]);
            for (int j=0; j< obs.spectrum_history[0].size(); j++)
                obs.spectrum_history[i][j] = t_spectrum_history(i,j);
            }
        
        std::cout << "focc_history... " << std::flush;
        triqs::arrays::array<double, 2> focc_history;
        h5_read(h5_mc_data,"focc_history", focc_history);
        obs.focc_history.resize(focc_history.shape()[0]);
        for (int i=0; i<obs.focc_history.size(); i++) { 
            obs.focc_history[i].resize(focc_history.shape()[1]);
            for (int j=0; j< obs.focc_history[0].size(); j++)
                focc_history(i,j) =  obs.focc_history[i][j];
            }
        };

    // Inverse participation ratio
    if (p["measure_ipr"] && p["measure_history"]) {
        std::cout << "ipr..." << std::flush;
        triqs::arrays::array<double, 2> t_ipr_history;
        h5_read(h5_mc_data,"ipr_history", t_ipr_history);
        obs.ipr_history.resize(t_ipr_history.shape()[0]);
        for (int i=0; i<obs.ipr_history.size(); i++) { 
            obs.ipr_history[i].resize(t_ipr_history.shape()[1]);
            for (int j=0; j< obs.ipr_history[0].size(); j++)
                obs.ipr_history[i][j] = t_ipr_history(i,j); 
            };
        };
    std::cout << "done." << std::endl;
    return std::move(obs);
}



} // end of namespace fk


