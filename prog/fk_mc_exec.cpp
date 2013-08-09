#include "fk_mc.hpp"
#include <boost/mpi/environment.hpp>
#include <chrono>
#include <tclap/CmdLine.h>
#include <fstream>
#include <triqs/h5.hpp>
#include <triqs/arrays/indexmaps/cuboid/domain.hpp>

#include "binning.hpp"

#ifndef FK_MC_DEBUG
#undef DEBUG
#endif

using namespace fk;


typedef fk_mc<triangular_lattice_traits> mc_t;

void print_section (const std::string& str); // fancy screen output
void save_data(const mc_t& mc, std::string output_file);

int main(int argc, char* argv[])
{
boost::mpi::environment env(argc, argv);
boost::mpi::communicator world;

try {
    TCLAP::CmdLine cmd("Falicov-Kimball Monte Carlo - parameters from command line", ' ', "");
    //TCLAP::ValueArg<std::string> nameArg("n","name","Name to print",true,"homer","string");
    /* Required model flags. */
    TCLAP::ValueArg<double> U_arg("U","U","value of U",false,4.0,"double",cmd);
    TCLAP::ValueArg<double> T_arg("T","T","Temperature",false,0.15,"double",cmd);
    TCLAP::ValueArg<size_t> L_arg("L","L","system size",false,4,"int",cmd);
    TCLAP::ValueArg<double> t_arg("t","t","hopping",false,-1.0,"double",cmd);
    TCLAP::ValueArg<double> tp_arg("","tp","next nearest hopping",false,-0.0,"double",cmd);

    /* Optional flags. */
    TCLAP::ValueArg<double> mu_arg("","mu","chemical potential",false,2.0,"double", cmd);
    //TCLAP::ValueArg<int> nc_arg("","nc","total number of c-electrons",false,4,"int", cmd);
    
    TCLAP::ValueArg<double> eps_f_arg("","ef","chemical potential",false,0.0,"double", cmd);
    //TCLAP::ValueArg<int> nf_arg("","nf","total number of f-electrons",false,4,"int", cmd);
    

    

    TCLAP::ValueArg<int> ncycles_arg("","ncycles","total number of cycles",false,50000,"int",cmd);
    TCLAP::ValueArg<int> nwarmup_arg("","nwarmup","Number of warmup cycles (no measure)",false,10000,"int",cmd);
    TCLAP::ValueArg<int> cycle_len_arg("l","cyclelen","Number of steps in one cycle",false,1,"int",cmd);
    TCLAP::SwitchArg     time_seed_switch("s","seed","Make a time seed flag", cmd, false);

    TCLAP::ValueArg<double>     move_flips_switch("","flip","Make flip (conserving)", false, 0.0, "double", cmd);
    TCLAP::ValueArg<double>     move_add_remove_switch("","addremove","Make add/remove step (non conserving)", false, 1.0, "double", cmd);
    TCLAP::ValueArg<double>     move_reshuffle_switch("","reshuffle","Make reshuffle step (non conserving)", false, 0.0, "double", cmd);

    TCLAP::SwitchArg exit_switch("","exit","Dry run", cmd, false);
    cmd.parse( argc, argv );

    INFO("Falicov-Kimball Monte Carlo");
    print_section("Model parameters:");
    size_t L = L_arg.getValue();     INFO2("System size                  : " << L);
    double U = U_arg.getValue();     INFO2("U                            : " << U);
    double t = t_arg.getValue();     INFO2("t                            : " << t);
    double tp = tp_arg.getValue();   INFO2("tp                           : " << tp);
    double T = T_arg.getValue();     INFO2("Temperature                  : " << T);

    double mu_c = (U_arg.isSet() && (!mu_arg.isSet())?U/2.:mu_arg.getValue()); 
                                     INFO2("Chemical potential (mu_c)    : " << mu_c);
    double mu_f = (U_arg.isSet() && (!eps_f_arg.isSet())?mu_c:mu_c + eps_f_arg.getValue()); 
                                     INFO2("e_f                          : " << eps_f_arg.getValue());
                                     INFO2("mu_f (mu_c + e_f)            : " << mu_f);
    double beta = 1.0/T;             INFO2("beta                         : " << beta);

    print_section("Monte Carlo parameters:");
    INFO2("Total number of cycles       : " << ncycles_arg.getValue()); 
    INFO2("MC steps in a cycle          : " << cycle_len_arg.getValue()); 
    INFO2("Warmup cycles                : " << nwarmup_arg.getValue()); 
    INFO2("Time-based seed              : " << time_seed_switch.getValue()); 
    INFO2("MC flip moves weight         : " << move_flips_switch.getValue());
    INFO2("MC add/remove moves weight   : " << move_add_remove_switch.getValue());
    INFO2("MC reshuffle moves weight    : " << move_reshuffle_switch.getValue());
    if (exit_switch.getValue()) exit(0);
    triangular_lattice_traits lattice(L);
    lattice.fill(t,tp);

    mc_t mc(lattice);

    triqs::utility::parameters p;
    p["U"] = U;
    p["mu_c"] = mu_c; p["mu_f"] = mu_f;
    p["beta"] = beta;
    p["Nf_start"] = L*L/2;
    p["Random_Generator_Name"] = ""; 
    p["Random_Seed"] = (34788+std::chrono::system_clock::now().time_since_epoch().count()*time_seed_switch.getValue());
    p["Verbosity"] = 3;
    p["Length_Cycle"] = cycle_len_arg.getValue(); 
    p["N_Warmup_Cycles"] = nwarmup_arg.getValue();
    p["N_Cycles"] = ncycles_arg.getValue();
    p["max_time"]=3600*5;
    
    p["mc_flip"] = move_flips_switch.getValue();
    p["mc_add_remove"] = move_add_remove_switch.getValue();
    p["mc_reshuffle"] = move_reshuffle_switch.getValue();

    mc.solve(p);

    world.barrier();
    if (world.rank() == 0) {
        save_data(mc,"output.h5");
        }
    }
    // Any exceptions related with command line parsing.
    catch (TCLAP::ArgException &e) {std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; } 
    catch (triqs::runtime_error const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
return 0;
}

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
  std::cout << std::string(str.size(),'=') << std::endl;
}

void save_binning(const binning::bin_data_t& binning_data, triqs::h5::group& h5_group, std::string name, bool save_plaintext = false)
{
    auto cor_lens = binning::calc_cor_length(binning_data);
    tqa::array<double, 2> data_arr(binning_data.size(),5);
    std::ofstream out; out.setf(std::ios::scientific); //out << std::setw(9);
    if (save_plaintext) out.open(name+"_binning.dat",std::ios::out);
    for (size_t i=0; i<binning_data.size(); i++){
        auto e = binning_data[i]; double c = cor_lens[i];
        std::array<double, 5> t ({double(std::get<binning::_SIZE>(e)), std::get<binning::_MEAN>(e), std::get<binning::_DISP>(e), std::get<binning::_SQERROR>(e), c });
        std::copy(t.begin(),t.end(),data_arr(i,tqa::range()).begin());
        if (save_plaintext) out << i << " " << std::get<binning::_SIZE>(e) << "  " << std::get<binning::_MEAN>(e) 
            << "  " << std::get<binning::_DISP>(e) << "  " << std::get<binning::_SQERROR>(e) << "  " << c << "\n";
       };
    h5_write(h5_group,name, data_arr);
    out.close();
}

void save_data(const mc_t& mc, std::string output_file)
{
    print_section("Statistics");
    H5::H5File output(output_file,H5F_ACC_TRUNC);
    triqs::h5::group top(output);
    top.create_group("mc_data");
    top.create_group("stats");

    auto h5_mc_data = top.open_group("mc_data");
    h5_write(h5_mc_data,"energies", mc.observables.energies);

    auto h5_stats = top.open_group("stats");
    INFO("Energy binning");
    size_t size = mc.observables.energies.size();
    int maxbin = std::min(15,std::max(int(std::log(size/16)/std::log(2.)-1),1));
    INFO("\tBinning " << maxbin <<" times.");
    auto energy_binning_data = binning::accumulate_binning(mc.observables.energies.rbegin(), mc.observables.energies.rend(), maxbin); 
    save_binning(energy_binning_data,h5_stats,"energies",true);
    auto d2energy_binning_data = binning::accumulate_binning(mc.observables.d2energies.rbegin(),mc.observables.d2energies.rend(), maxbin); 
    save_binning(d2energy_binning_data,h5_stats,"d2energies",true);
}

