#include <boost/mpi/environment.hpp>
#include <chrono>
#include <tclap/CmdLine.h>

//#include <triqs/arrays/indexmaps/cuboid/domain.hpp>
#include <random>

#include "fk_mc.hpp"
#include "data_saveload.hpp"
#include "measures/polarization.hpp"


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
#elif LATTICE_chain
    #include "lattice/chain.hpp"
    typedef chain_lattice lattice_t;
#endif


size_t _myrank;

#define MINFO(MSG)            if (_myrank==0) std::cout << std::boolalpha << MSG << std::endl;
#define MINFO2(MSG)            if (_myrank==0) std::cout << "    " << std::boolalpha << MSG << std::endl;

using namespace std::chrono;

int main(int argc, char* argv[])
{
boost::mpi::environment env(argc, argv);
boost::mpi::communicator world;
_myrank = world.rank();

try {
    TCLAP::CmdLine cmd("Falicov-Kimball Monte Carlo - parameters from command line", ' ', "");
    //TCLAP::ValueArg<std::string> nameArg("n","name","Name to print",true,"homer","string");
    // Required model flags.
    TCLAP::ValueArg<double> U_arg("U","U","value of U",false,1.0,"double",cmd);
    TCLAP::ValueArg<double> T_arg("T","T","Temperature",false,0.1,"double",cmd);
    TCLAP::ValueArg<size_t> L_arg("L","L","system size",false,4,"int",cmd);
    TCLAP::ValueArg<double> t_arg("t","t","hopping",false,1.0,"double",cmd);
    #ifdef LATTICE_triangular
        TCLAP::ValueArg<double> tp_arg("","tp","next nearest hopping",false,0.0,"double",cmd);
    #elif LATTICE_chain
        TCLAP::ValueArg<double> delta_arg("","delta","delta",false,0.0,"double",cmd);
        TCLAP::ValueArg<double> eta_arg("","eta","eta",false,0.0,"double",cmd);
    #endif

    TCLAP::ValueArg<std::string> h5file_arg("o","output","archive to read/write data to",false,"output.h5","string",cmd);
    TCLAP::SwitchArg             resume_switch("r","resume","Attepmpt to resume a calculation", cmd, false);

    // Optional flags.
    TCLAP::ValueArg<double> mu_arg("","mu","chemical potential",false,0.5,"double", cmd);
    //TCLAP::ValueArg<int> nc_arg("","nc","total number of c-electrons",false,4,"int", cmd);
    TCLAP::ValueArg<double> eps_f_arg("","ef","chemical potential",false,0.0,"double", cmd);
    //TCLAP::ValueArg<int> nf_arg("","nf","total number of f-electrons",false,4,"int", cmd);
    

    TCLAP::ValueArg<int>      ncycles_arg("","ncycles","total number of cycles",false,100,"int",cmd);
    TCLAP::ValueArg<int>      nwarmup_arg("","nwarmup","Number of warmup cycles (no measure)",false,1,"int",cmd);
    TCLAP::ValueArg<int>      cycle_len_arg("l","cyclelen","Number of steps in one cycle",false,100,"int",cmd);
    TCLAP::SwitchArg          random_seed_switch("s","seed","Make a random or fixed seed?", cmd, false);

    TCLAP::ValueArg<double>   move_flips_switch("","flip","Make flip (conserving)", false, 0.0, "double", cmd);
    TCLAP::ValueArg<double>   move_add_remove_switch("","addremove","Make add/remove step (non conserving)", false, 1.0, "double", cmd);
    TCLAP::ValueArg<double>   move_reshuffle_switch("","reshuffle","Make reshuffle step (non conserving)", false, 0.0, "double", cmd);

    TCLAP::SwitchArg          plaintext_switch("p","plaintext","Save data to plaintext format?", cmd, false);
    TCLAP::ValueArg<int>      maxtime_arg("","maxtime","Maximum evaluation time (in hours)", false, 24*30, "int", cmd);

    TCLAP::ValueArg<bool>     calc_history_switch("","calc_history","Calculate data history (for errorbars)", false, false, "bool", cmd);
    TCLAP::ValueArg<bool>     calc_ipr_switch("","calc_ipr","Calculate inverse participation ratio", false, false, "bool", cmd);
    // dos-related args
    TCLAP::ValueArg<double>   dos_width_arg("","dos_width","width of dos", false, 6.0, "double", cmd);
    TCLAP::ValueArg<int>      dos_npts_arg("","dos_npts","npts dos", false, 1000, "int", cmd);
    TCLAP::ValueArg<double>   dos_offset_arg("","dos_offset","offset of dos from real axis", false, 0.05, "double", cmd);

    // chebyshev flags
    TCLAP::SwitchArg          chebyshev_switch("c","chebyshev","Make chebyshev moves?", cmd, false);
    TCLAP::ValueArg<double>   chebyshev_prefactor("","cheb_prefactor","Prefactor of log(N) chebyshev polynomials", false, 2.2, "double", cmd);

    TCLAP::SwitchArg exit_switch("","exit","Dry run", cmd, false);
    cmd.parse( argc, argv );

    MINFO("Falicov-Kimball Monte Carlo");
    std::unique_ptr<fk_mc<lattice_t>> mc_ptr;
    triqs::utility::parameters p;

    double T = T_arg.getValue();    
    double beta = 1.0/T;
    double U = U_arg.getValue();

    // convert input to triqs::parameters
    p["L"] = L_arg.getValue();
    p["U"] = U_arg.getValue();
    p["t"] = t_arg.getValue();
    p["beta"] = beta;
    p["t"] = t_arg.getValue();
    double mu_c = (U_arg.isSet() && (!mu_arg.isSet())?U/2.:mu_arg.getValue()); 
    p["mu_c"] = mu_c; 
    double mu_f = (U_arg.isSet() && (!eps_f_arg.isSet())?mu_c:mu_c + eps_f_arg.getValue()); 
    p["mu_f"] = mu_f;
    //p["random_name"] = ""; 
    p["random_seed"] = (random_seed_switch.getValue()?std::random_device()():(32167+world.rank()));
    p["verbosity"] = (!world.rank()?1:0);
    p["length_cycle"] = cycle_len_arg.getValue(); 
    p["n_warmup_cycles"] = nwarmup_arg.getValue();
    p["n_cycles"] = ncycles_arg.getValue();
    p["max_time"]=maxtime_arg.getValue()*3600;

    p["measure_ipr"] = calc_ipr_switch.getValue();
    p["measure_history"] = calc_history_switch.getValue() || bool(p["measure_ipr"]);
    p["dos_width"] = dos_width_arg.getValue();
    p["dos_npts"] = dos_npts_arg.getValue();
    p["dos_offset"] = dos_offset_arg.getValue();
    
    p["mc_flip"] = move_flips_switch.getValue();
    p["mc_add_remove"] = move_add_remove_switch.getValue();
    p["mc_reshuffle"] = move_reshuffle_switch.getValue();

    p["cheb_moves"] = chebyshev_switch.getValue(); 
    p["cheb_prefactor"] = chebyshev_prefactor.getValue();
    //if (!world.rank()) std::cout << p << std::endl;

    lattice_t lattice(p["L"]);
    if (resume_switch.getValue()) { 
        if (!world.rank()) std::cout << "Resuming calculation" << std::endl;
        mc_ptr.reset( new fk_mc<lattice_t>(load_data<fk_mc<lattice_t>>(h5file_arg.getValue(),p)) );
        exit(0);
        }
    else { 
        size_t L = p["L"];
        double U = p["U"]; 
        double t = p["t"]; 

        if (!_myrank) print_section("Model parameters:");
        MINFO2("System size                  : " << L);
        MINFO2("U                            : " << U);
        MINFO2("t                            : " << t);
        MINFO2("Temperature                  : " << T);

        MINFO2("Chemical potential (mu_c)    : " << mu_c);
        MINFO2("e_f                          : " << eps_f_arg.getValue());
        MINFO2("mu_f (mu_c + e_f)            : " << mu_f);
        MINFO2("beta                         : " << beta);

        if (!_myrank) print_section("Monte Carlo parameters:");
        MINFO2("Total number of cycles       : " << ncycles_arg.getValue()); 
        MINFO2("MC steps in a cycle          : " << cycle_len_arg.getValue()); 
        MINFO2("Warmup cycles                : " << nwarmup_arg.getValue()); 
        MINFO2("Pure random (entropic) seed  : " << random_seed_switch.getValue()); 
        MINFO2("MC flip moves weight         : " << move_flips_switch.getValue());
        MINFO2("MC add/remove moves weight   : " << move_add_remove_switch.getValue());
        MINFO2("MC reshuffle moves weight    : " << move_reshuffle_switch.getValue());
        if (exit_switch.getValue()) exit(0);


        #ifdef LATTICE_triangular
            double tp = tp_arg.getValue();   MINFO2("tp                           : " << tp);
            lattice.fill(t,tp);
            p["Nf_start"] = L*L/2;
            p["tp"] = tp;
        #elif LATTICE_chain
            lattice.fill(t,eta_arg.getValue(),delta_arg.getValue());
            p["Nf_start"] = L/2;
        #elif LATTICE_cubic1d 
            lattice.fill(t);
            p["Nf_start"] = L/2;
        #elif LATTICE_cubic2d 
            lattice.fill(t);
            p["Nf_start"] = L*L/2;
        #endif

        if (!world.rank()) std::cout << "All parameters: " << p << std::endl;
        mc_ptr.reset(new fk_mc<lattice_t> (lattice,p));
        }; // end - noresume

    fk_mc<lattice_t>& mc = *mc_ptr;

    #ifdef LATTICE_chain
    mc.add_measure(measure_polarization<lattice_t>(mc.config,lattice),"polarization");
    #endif
        
    steady_clock::time_point start, end;
    start = steady_clock::now();
    mc.solve();
    end = steady_clock::now();

    world.barrier();
    if (world.rank() == 0) {
        save_data(mc,p,h5file_arg.getValue(),plaintext_switch.getValue());
        std::cout << "Calculation lasted : " 
            << duration_cast<hours>(end-start).count() << "h " 
            << duration_cast<minutes>(end-start).count()%60 << "m " 
            << duration_cast<seconds>(end-start).count()%60 << "s " 
            << duration_cast<milliseconds>(end-start).count()%1000 << "ms " 
            << std::endl;
        }
    mc_ptr.release();
    }
    // Any exceptions related with command line parsing.
    catch (TCLAP::ArgException &e) {std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; } 
    catch (triqs::runtime_error const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
return 0;
}



