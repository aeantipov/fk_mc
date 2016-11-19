#include <boost/mpi/environment.hpp>
#include <chrono>
#include <random>

#include "fk_mc.hpp"
#include "data_save.hpp"
//#include "data_load.hpp"
#include "measures/polarization.hpp"
#include <alps/mc/mpiadapter.hpp>
#include <alps/mc/stop_callback.hpp>


using namespace fk;

#ifdef LATTICE_triangular
    #include "lattice/triangular.hpp"
    typedef triangular_lattice lattice_t;
#elif LATTICE_cubic1d
    #include "lattice/hypercubic.hpp"
    typedef hypercubic_lattice<1> lattice_t;
#elif LATTICE_cubic2d
    #include "lattice/hypercubic.hpp"
    typedef hypercubic_lattice<2> lattice_t;
#elif LATTICE_cubic3d
    #include "lattice/hypercubic.hpp"
    typedef hypercubic_lattice<3> lattice_t;
#elif LATTICE_chain
    #include "lattice/chain.hpp"
    typedef chain_lattice lattice_t;
#elif LATTICE_honeycomb
    #include "lattice/honeycomb.hpp"
    typedef honeycomb_lattice lattice_t;
#endif


size_t _myrank;

#define MINFO(MSG)            if (_myrank==0) std::cout << std::boolalpha << MSG << std::endl;
#define MINFO2(MSG)            if (_myrank==0) std::cout << "    " << std::boolalpha << MSG << std::endl;
#define mpi_cout if (!comm.rank()) std::cout
#define mpi_cerr if (!comm.rank()) std::cerr
#define mpi_exit() MPI_Finalize()

using namespace std::chrono;

typedef alps::mcmpiadapter<fk_mc<lattice_t>> qmc_t;

// params from command line
alps::params cmdline_params(int argc, char *argv[]);

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    _myrank = comm.rank();
        
    //alps::params p;
    alps::params p(argc, (const char **) argv);

    // Parameters from command line -> auto convert to alps::params
    try { p = cmdline_params(argc, argv); }
    catch (std::exception &e) {
        mpi_exit();
        std::cout << e.what() << std::endl;
        exit(1);
    };
    if (p.help_requested(std::cerr)) {
        mpi_exit();
        exit(1);
    };

    mpi_cout << "All parameters : " << std::endl << p << std::endl;

    MINFO("Falicov-Kimball Monte Carlo");

    double beta = p["beta"]; 
    double U = p["U"]; 
    bool resume = false; 
    bool dry_run = p["exit"];
    //p["random_seed"] = (random_seed_switch.getValue()?std::random_device()():(rnd_seed_arg.getValue()+comm.rank()));

    p["measure_history"] = (p["measure_history"].as<bool>() && !bool(p["cheb_moves"])) || bool(p["measure_ipr"]) || bool(p["measure_eigenfunctions"]);
    
    int nsweeps_new = std::max(int(p["nsweeps"]), 0);
    int nsweeps_old = 0; 

// Now we construct and run mc 
try{
    // try resuming calc
    observables_t obs_old;
/*
    if (resume) {
        p = load_parameters(h5fname, p);
        nsweeps_old = p["nsweeps"]; // this is old + new number of sweeps
        if (!comm.rank()) {
             std::cout << "Resuming calculation from " << nsweeps_old << " sweeps. "<< std::endl;
             obs_old = load_observables(h5fname, p);
             std::cout << "parameters for run : " << p << std::endl;
        };
        comm.barrier();
       // boost::mpi::broadcast(comm, p, 0); 
        //std::cout << "proc " << comm.rank() << " parameters : " << p << std::endl;
        //exit(0);
        }

*/
    int nsweeps_total = std::max(nsweeps_new, nsweeps_old);
    p["nsweeps"] = std::max(nsweeps_total - nsweeps_old, 0); // do a calc with only new number of sweeps
    if (!comm.rank()) std::cout << nsweeps_old << " -> " << nsweeps_total << " sweeps" << std::endl; 
    size_t L = p["L"];
    double U = p["U"]; 
    double t = p["t"]; 
    double beta = p["beta"];
    double mu_c = p["mu_c"];
    double mu_f = p["mu_f"];

    if (!_myrank) print_section("Model parameters:");
    MINFO2("System size                  : " << L);
    MINFO2("U                            : " << U);
    MINFO2("t                            : " << t);
    MINFO2("Temperature                  : " << 1.0/beta);

    MINFO2("Chemical potential (mu_c)    : " << mu_c);
    MINFO2("e_f                          : " << mu_f - mu_c);
    MINFO2("mu_f (mu_c + e_f)            : " << mu_f);
    MINFO2("beta                         : " << beta);

    if (!_myrank) print_section("Monte Carlo parameters:");
    MINFO2("Total number of sweeps       : " << p["nsweeps"]); 
    MINFO2("MC steps in a sweep          : " << p["sweep_len"]); 
    MINFO2("Warmup sweeps                : " << p["ntherm_sweeps"]);
    MINFO2("Random seed                  : " << p["seed"]);
    MINFO2("MC flip moves weight         : " << p["mc_flip"]); 
    MINFO2("MC add/remove moves weight   : " << p["mc_add_remove"]);
    MINFO2("MC reshuffle moves weight    : " << p["mc_reshuffle"]);

    if (dry_run) exit(0);

    lattice_t lattice(p["L"].as<size_t>());
    #ifdef LATTICE_triangular
        MINFO2("tp                           : " << p["tp"]);
        lattice.fill(double(p["t"]),double(p["tp"]));
        p["Nf_start"] = int(L*L/2);
    #elif LATTICE_chain
        lattice.fill(t,eta_arg.getValue(),delta_arg.getValue());
        p["Nf_start"] = int(L/2);
    #elif LATTICE_cubic1d 
        lattice.fill(p["t"]);
        p["Nf_start"] = int(L/2);
    #elif LATTICE_cubic2d 
        lattice.fill(p["t"]);
        p["Nf_start"] = int(L*L/2);
    #elif LATTICE_cubic3d 
        lattice.fill(p["t"]);
        p["Nf_start"] = int(L*L*L/2);
    #elif LATTICE_honeycomb 
        lattice.fill(p["t"]);
        p["Nf_start"] = int(L*L/2);
    #endif

    if (!comm.rank()) std::cout << "All parameters: " << p << std::endl;
    
    // create a log grid for conductivity
    int nw_size = p["cond_npoints"];
    double wmax = std::max(8.0, 2*double(p["U"]));
    double base = M_E;
    int min_power = -15;
    Eigen::VectorXd wgrid1 (2*nw_size + 1);
    { 
        Eigen::VectorXd wgrid2 = Eigen::VectorXd::LinSpaced(nw_size, min_power, 0);
        for (int i=0; i<wgrid2.size(); i++) { wgrid2[i] = wmax * std::pow(base, wgrid2[i]); }
        wgrid1 << -wgrid2.reverse(), Eigen::VectorXd::Zero(1), wgrid2;
    }
        
    std::vector<double> wgrid_conductivity ({wgrid1.data(), wgrid1.data() + wgrid1.size()});

    fk_mc<lattice_t> mc(p, _myrank); 
    mc.initialize(lattice, true, wgrid_conductivity);

    //#ifdef LATTICE_chain
    //mc.add_measure(measure_polarization<lattice_t>(mc.config,mc.lattice),"polarization");
    //#endif
        
    steady_clock::time_point start, end;
    start = steady_clock::now();
    mc.run(alps::stop_callback(p["max_time"].as<size_t>())); // this runs monte-carlo
    end = steady_clock::now();

    comm.barrier();
    if (comm.rank() == 0) {
        mc.observables.merge(obs_old);
        p["nsweeps"] = nsweeps_total;
        std::cout << "Calculation lasted : " 
            << duration_cast<hours>(end-start).count() << "h " 
            << duration_cast<minutes>(end-start).count()%60 << "m " 
            << duration_cast<seconds>(end-start).count()%60 << "s " 
            << duration_cast<milliseconds>(end-start).count()%1000 << "ms " 
            << std::endl;

        start = steady_clock::now();
        save_all_data(mc,p,wgrid_conductivity);
        end = steady_clock::now();
        std::cout << "Saving lasted : " 
            << duration_cast<hours>(end-start).count() << "h " 
            << duration_cast<minutes>(end-start).count()%60 << "m " 
            << duration_cast<seconds>(end-start).count()%60 << "s " 
            << duration_cast<milliseconds>(end-start).count()%1000 << "ms " 
            << std::endl;
        }
    }
    // Any exceptions related with command line parsing.
    catch (std::exception const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
return 0;
}


alps::params cmdline_params(int argc, char *argv[]) {
    alps::params p(argc, (const char **) argv);

    p.description("Falicov-Kimball Monte Carlo - parameters from command line");

    qmc_t::define_parameters(p);

    p.define<size_t> ("L", 4, "System linear size");
    p.define<double> ("t", 1.0, "Hopping");

    #ifdef LATTICE_triangular
        p.define<double> ("tp", 1.0, "Triangular lattice : NNN Hopping");
    #elif LATTICE_chain
        p.define<double> ("delta", 0.0, "chain : delta");
        p.define<double> ("eta", 0.0, "chain : eta");
    #endif

    p.define<std::string>("output", "output.h5", "archive to read/write data to");
    p.define<bool>("plaintext", false, "plaintext output level");
    p.define<int>("maxtime", 24*30, "max evaluation time");
    // DOS args
    p.define<double>("dos_width", 6.0, "Width of DOS");
    p.define<int>("dos_npts", 240, "Number of points for DOS sampling");
    p.define<double>("dos_offset", 0.05, "DOS offset from real axis");
    // stiffness args
    p.define<int>("cond_npoints", 150, "number of points to sample conductivity");
    // eigenfunctions storage
    p.define<bool>("save_eigenfunctions", false, "Store eigenfunctions?");

    p.define<bool>("exit", false, "dry_run");
    p.define<size_t>("max_time",size_t(14*24*3600), "Maximum running time in seconds");

    return p;
}


