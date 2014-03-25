#include <boost/mpi/environment.hpp>
#include <chrono>
#include <tclap/CmdLine.h>
#include <fstream>
#include <triqs/h5.hpp>
#include <triqs/arrays/indexmaps/cuboid/domain.hpp>
#include <random>

#include "fk_mc.hpp"
#include "measures/fsusc0pi.hpp"

#include "binning.hpp"
#include "jackknife.hpp"

using namespace fk;

#ifdef LATTICE_triangular
    typedef triangular_lattice lattice_t;
#endif

size_t _myrank;

#define MINFO(MSG)            if (_myrank==0) std::cout << std::boolalpha << MSG << std::endl;
#define MINFO2(MSG)            if (_myrank==0) std::cout << "    " << std::boolalpha << MSG << std::endl;

void print_section (const std::string& str); // fancy screen output
void savetxt (std::string fname, const triqs::arrays::array<double,1>& in);
void savetxt (std::string fname, const triqs::arrays::array<double,2>& in);
void save_data(const fk_mc& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false);

int main(int argc, char* argv[])
{
boost::mpi::environment env(argc, argv);
boost::mpi::communicator world;
_myrank = world.rank();

try {
    TCLAP::CmdLine cmd("Falicov-Kimball Monte Carlo - parameters from command line", ' ', "");
    //TCLAP::ValueArg<std::string> nameArg("n","name","Name to print",true,"homer","string");
    /* Required model flags. */
    TCLAP::ValueArg<double> U_arg("U","U","value of U",false,1.0,"double",cmd);
    TCLAP::ValueArg<double> T_arg("T","T","Temperature",false,0.1,"double",cmd);
    TCLAP::ValueArg<size_t> L_arg("L","L","system size",false,4,"int",cmd);
    TCLAP::ValueArg<double> t_arg("t","t","hopping",false,1.0,"double",cmd);
    TCLAP::ValueArg<double> tp_arg("","tp","next nearest hopping",false,0.0,"double",cmd);

    /* Optional flags. */
    TCLAP::ValueArg<double> mu_arg("","mu","chemical potential",false,0.5,"double", cmd);
    //TCLAP::ValueArg<int> nc_arg("","nc","total number of c-electrons",false,4,"int", cmd);
    
    TCLAP::ValueArg<double> eps_f_arg("","ef","chemical potential",false,0.0,"double", cmd);
    //TCLAP::ValueArg<int> nf_arg("","nf","total number of f-electrons",false,4,"int", cmd);
    

    TCLAP::ValueArg<int> ncycles_arg("","ncycles","total number of cycles",false,50000,"int",cmd);
    TCLAP::ValueArg<int> nwarmup_arg("","nwarmup","Number of warmup cycles (no measure)",false,0,"int",cmd);
    TCLAP::ValueArg<int> cycle_len_arg("l","cyclelen","Number of steps in one cycle",false,1,"int",cmd);
    TCLAP::SwitchArg     random_seed_switch("s","seed","Make a random or fixed seed?", cmd, false);

    TCLAP::ValueArg<double>     move_flips_switch("","flip","Make flip (conserving)", false, 0.0, "double", cmd);
    TCLAP::ValueArg<double>     move_add_remove_switch("","addremove","Make add/remove step (non conserving)", false, 1.0, "double", cmd);
    TCLAP::ValueArg<double>     move_reshuffle_switch("","reshuffle","Make reshuffle step (non conserving)", false, 0.0, "double", cmd);

    //TCLAP::ValueArg<double>     eval_tolerance_switch("","evaltol","Tolerance for eigenvalue weights", false, std::numeric_limits<double>::epsilon(), "double", cmd);
    TCLAP::SwitchArg     plaintext_switch("p","plaintext","Save data to plaintext format?", cmd, false);

    TCLAP::ValueArg<bool>     calc_history_switch("","calc_history","Calculate data history (for errorbars)", false, true, "bool", cmd);
    // dos-related args
    TCLAP::ValueArg<double>     dos_width_arg("","dos_width","width of dos", false, 6.0, "double", cmd);
    TCLAP::ValueArg<int>        dos_npts_arg("","dos_npts","npts dos", false, 1000, "int", cmd);
    TCLAP::ValueArg<double>     dos_offset_arg("","dos_offset","offset of dos from real axis", false, 0.05, "double", cmd);

    TCLAP::SwitchArg exit_switch("","exit","Dry run", cmd, false);
    cmd.parse( argc, argv );

    MINFO("Falicov-Kimball Monte Carlo");
    print_section("Model parameters:");
    size_t L = L_arg.getValue();     MINFO2("System size                  : " << L);
    double U = U_arg.getValue();     MINFO2("U                            : " << U);
    double t = t_arg.getValue();     MINFO2("t                            : " << t);
    double tp = tp_arg.getValue();   MINFO2("tp                           : " << tp);
    double T = T_arg.getValue();     MINFO2("Temperature                  : " << T);

    double mu_c = (U_arg.isSet() && (!mu_arg.isSet())?U/2.:mu_arg.getValue()); 
                                     MINFO2("Chemical potential (mu_c)    : " << mu_c);
    double mu_f = (U_arg.isSet() && (!eps_f_arg.isSet())?mu_c:mu_c + eps_f_arg.getValue()); 
                                     MINFO2("e_f                          : " << eps_f_arg.getValue());
                                     MINFO2("mu_f (mu_c + e_f)            : " << mu_f);
    double beta = 1.0/T;             MINFO2("beta                         : " << beta);

    print_section("Monte Carlo parameters:");
    MINFO2("Total number of cycles       : " << ncycles_arg.getValue()); 
    MINFO2("MC steps in a cycle          : " << cycle_len_arg.getValue()); 
    MINFO2("Warmup cycles                : " << nwarmup_arg.getValue()); 
    MINFO2("Pure random (entropic) seed  : " << random_seed_switch.getValue()); 
    MINFO2("MC flip moves weight         : " << move_flips_switch.getValue());
    MINFO2("MC add/remove moves weight   : " << move_add_remove_switch.getValue());
    MINFO2("MC reshuffle moves weight    : " << move_reshuffle_switch.getValue());
    //MINFO2("Eigenvalue Boltzmann weight cutoff    : " << eval_tolerance_switch.getValue());
    if (exit_switch.getValue()) exit(0);
    lattice_t lattice(L); // create a lattice
    lattice.fill(t,tp);


    triqs::utility::parameters p;
    p["U"] = U;
    p["mu_c"] = mu_c; p["mu_f"] = mu_f;
    p["beta"] = beta;
    p["Nf_start"] = L*L/2;
    p["random_name"] = ""; 
    //p["eval_tol"] = eval_tolerance_switch.getValue(); 
    p["random_seed"] = (random_seed_switch.getValue()?std::random_device()():(34788+world.rank()));
    p["verbosity"] = (!world.rank()?3:0);
    p["length_cycle"] = cycle_len_arg.getValue(); 
    p["n_warmup_cycles"] = nwarmup_arg.getValue();
    p["n_cycles"] = ncycles_arg.getValue();
    p["max_time"]=3600*5;

    p["measure_history"] = calc_history_switch.getValue();
    p["dos_width"] = dos_width_arg.getValue();
    p["dos_npts"] = dos_npts_arg.getValue();
    p["dos_offset"] = dos_offset_arg.getValue();
    
    p["mc_flip"] = move_flips_switch.getValue();
    p["mc_add_remove"] = move_add_remove_switch.getValue();
    p["mc_reshuffle"] = move_reshuffle_switch.getValue();

    INFO(p);

    fk_mc mc(lattice,p);
    mc.add_measure(measure_nf0pi<lattice_t>(mc.config, lattice, mc.observables.nf0, mc.observables.nfpi), "nf0pi");
    mc.solve();

    world.barrier();
    if (world.rank() == 0) {
        save_data(mc,p,"output.h5",plaintext_switch.getValue());
        }
    }
    // Any exceptions related with command line parsing.
    catch (TCLAP::ArgException &e) {std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; } 
    catch (triqs::runtime_error const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
return 0;
}

void print_section (const std::string& str)
{
  MINFO(std::string(str.size(),'='));
  MINFO(str)
  MINFO(std::string(str.size(),'='));
}

void save_bin_data(const binning::bin_stats_t& data, triqs::h5::group& h5_group, std::string name, bool save_plaintext = false)
{
    std::array<double, 4> tmp (
        {double(std::get<binning::_SIZE>(data)), std::get<binning::_MEAN>(data), std::get<binning::_DISP>(data), std::get<binning::_SQERROR>(data) });
    tqa::array<double, 1> data_arr(4);
    std::copy(tmp.begin(),tmp.end(),data_arr(tqa::range()).begin());

    std::cout << name << " = " << std::get<binning::_MEAN>(data) << " +/- " << std::get<binning::_SQERROR>(data) 
              << " (" << std::get<binning::_SIZE>(data) << ") samples." << std::endl ;
    h5_write(h5_group, name, data_arr);
    if (save_plaintext) savetxt(name+"_error.dat",data_arr);
}
void save_binning(const binning::bin_data_t& binning_data, triqs::h5::group& h5_group, std::string name, bool save_plaintext = false)
{
    auto cor_lens = binning::calc_cor_length(binning_data);
    tqa::array<double, 2> data_arr(binning_data.size(),5);
    //std::ofstream out; out.setf(std::ios::scientific); //out << std::setw(9);
    //if (save_plaintext) out.open(name+"_binning.dat",std::ios::out);
    for (size_t i=0; i<binning_data.size(); i++){
        auto e = binning_data[i]; double c = cor_lens[i];
        std::array<double, 5> t ({double(std::get<binning::_SIZE>(e)), std::get<binning::_MEAN>(e), std::get<binning::_DISP>(e), std::get<binning::_SQERROR>(e), c });
        std::copy(t.begin(),t.end(),data_arr(i,tqa::range()).begin());
        //if (save_plaintext) out << i << " " << std::get<binning::_SIZE>(e) << "  " << std::get<binning::_MEAN>(e) 
            //<< "  " << std::get<binning::_DISP>(e) << "  " << std::get<binning::_SQERROR>(e) << "  " << c << "\n";
       };
    if (save_plaintext) savetxt(name+"_binning.dat",data_arr);
    h5_write(h5_group,name, data_arr);
}

/** Estimate the bin, at which the error bar is saturated. */
size_t estimate_bin(const fk::binning::bin_data_t& data)
{
    std::vector<double> errors(data.size());
    for (size_t i=0; i<errors.size(); i++) errors[i] = std::get<binning::bin_m::_SQERROR>(data[i]);
    double rel_error = 1.;
    bool f = true; 
    size_t ind = errors.size()-1;
    while (f && ind>0){
        double rel_error_current = std::abs(errors[ind-1]/errors[ind]-1.);
        f = (rel_error_current < 0.05 && rel_error_current < rel_error);
        rel_error = f?rel_error_current:rel_error;
        if (f) ind--;
    }
    return ind;
}

void save_data(const fk_mc& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext)
{

    double beta = p["beta"];
    double Volume = mc.lattice.get_msize();

    print_section("Statistics");
    H5::H5File output(output_file.c_str(),H5F_ACC_TRUNC);
    triqs::h5::group top(output);

    //===== save direct measures ===== //
    top.create_group("mc_data");
    auto h5_mc_data = top.open_group("mc_data");
    h5_write(h5_mc_data,"energies", mc.observables.energies);
    h5_write(h5_mc_data,"d2energies", mc.observables.d2energies);
    h5_write(h5_mc_data,"nf0", mc.observables.nf0);
    h5_write(h5_mc_data,"nfpi", mc.observables.nfpi);

    std::vector<double> spectrum(mc.observables.spectrum.size());
    std::copy(mc.observables.spectrum.data(), mc.observables.spectrum.data()+spectrum.size(), spectrum.begin());
    h5_write(h5_mc_data,"spectrum", spectrum);


    if (p["measure_history"]) { 
        triqs::arrays::array<double, 2> t_spectrum_history(mc.observables.spectrum_history.size(), mc.observables.spectrum_history[0].size());
        for (int i=0; i<mc.observables.spectrum_history.size(); i++)
            for (int j=0; j< mc.observables.spectrum_history[0].size(); j++)
                t_spectrum_history(i,j) =  mc.observables.spectrum_history[i][j];
        h5_write(h5_mc_data,"spectrum_history", t_spectrum_history);

        triqs::arrays::array<double, 2> focc_history(mc.observables.focc_history.size(), mc.observables.focc_history[0].size());
        for (int i=0; i<mc.observables.focc_history.size(); i++)
            for (int j=0; j< mc.observables.focc_history[0].size(); j++)
                focc_history(i,j) =  mc.observables.focc_history[i][j];
        h5_write(h5_mc_data,"focc_history", focc_history);
        };



    //===== save statistics ===== //
    top.create_group("stats");
    auto h5_stats = top.open_group("stats");
    top.create_group("binning");
    auto h5_binning = top.open_group("binning");

    const std::vector<double>& energies = mc.observables.energies;
    const std::vector<double>& d2energies = mc.observables.d2energies;
    int maxbin = std::min(15,std::max(int(std::log(energies.size()/16)/std::log(2.)-1),1));
    size_t energy_bin = 0;

    { /* Energy binning */
        INFO("Energy binning");
        size_t size = energies.size();
        INFO("Binning " << maxbin <<" times.");
        auto energy_binning_data = binning::accumulate_binning(energies.rbegin(), energies.rend(), maxbin); 
        save_binning(energy_binning_data,h5_binning,"energies",save_plaintext);
        auto d2energy_binning_data = binning::accumulate_binning(d2energies.rbegin(),d2energies.rend(), maxbin); 
        save_binning(d2energy_binning_data,h5_binning,"d2energies",save_plaintext);

        energy_bin = estimate_bin(energy_binning_data);
        save_bin_data(energy_binning_data[energy_bin],h5_stats,"energy",save_plaintext);
        save_bin_data(d2energy_binning_data[energy_bin],h5_stats,"d2energy",save_plaintext);
    };
    
    std::vector<double> energies_square(energies.size());
    std::transform(energies.begin(), energies.end(), energies_square.begin(), [](double x){ return x*x; });
    { /* Specific heat */
        typedef std::function<double(double, double, double)> cf_t;

        cf_t cv_function = [beta, Volume](double e, double e2, double de2){return beta*beta*(e2 - de2 - e*e)/Volume;}; 

        typedef decltype(energies.rbegin()) it_t;
        std::vector<std::pair<it_t,it_t>> c_data = {
            std::make_pair(energies.rbegin(), energies.rend()), 
            std::make_pair(energies_square.rbegin(), energies_square.rend()), 
            std::make_pair(d2energies.rbegin(), d2energies.rend())
        }; 
        auto cv_stats = jackknife::accumulate_jackknife(cv_function,c_data,maxbin);
        save_binning(cv_stats,h5_binning,"cv",save_plaintext);
        save_bin_data(cv_stats[energy_bin],h5_stats,"cv",save_plaintext);
    };

    {   /* Save nf(q=0), nf(q=pi) */
        const std::vector<double>& nf0 = mc.observables.nf0;
        const std::vector<double>& nfpi = mc.observables.nfpi;
        std::vector<double> nn_0(nf0.size()), nn_pi(nfpi.size());
        std::transform(nf0.begin(), nf0.end(), nn_0.begin(), [](double x){return x*x;});
        std::transform(nfpi.begin(), nfpi.end(), nn_pi.begin(), [](double x){return x*x;});
        std::function<double(double,double)> disp_f = [](double x, double x2){return x2 - x*x;};
        
        auto nf0_stats = binning::accumulate_binning(nf0.rbegin(), nf0.rend(), maxbin);
        save_binning( nf0_stats, h5_binning,"nf_0",save_plaintext);
        auto nfpi_stats = binning::accumulate_binning(nfpi.rbegin(), nfpi.rend(), maxbin);
        save_binning( nfpi_stats, h5_binning,"nf_pi",save_plaintext);
        auto fsusc0_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nf0,nn_0}),maxbin);
        save_binning(fsusc0_stats,h5_binning,"fsusc_0",save_plaintext);
        auto fsuscpi_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nfpi,nn_pi}),maxbin);
        save_binning(fsuscpi_stats,h5_binning,"fsusc_pi",save_plaintext);

        auto nf_bin = estimate_bin(fsusc0_stats);
        save_bin_data(nf0_stats[nf_bin],h5_stats,"nf_0",save_plaintext);
        save_bin_data(fsusc0_stats[nf_bin],h5_stats,"fsusc_0",save_plaintext);
        save_bin_data(nfpi_stats[nf_bin],h5_stats,"nf_pi",save_plaintext);
        save_bin_data(fsuscpi_stats[nf_bin],h5_stats,"fsusc_pi",save_plaintext);
    
        /* Binder cumulant (q=0, pi). */
        std::vector<double> nnnn_0(nf0.size()), nnnn_pi(nfpi.size());
        std::transform(nf0.begin(), nf0.end(), nnnn_0.begin(), [](double x){return x*x*x*x;});
        std::transform(nfpi.begin(), nfpi.end(), nnnn_pi.begin(), [](double x){return x*x*x*x;});
        std::function<double(double,double)> binder_f = [](double x2, double x4){return 1. - x4/3./x2/x2;};
        auto binder_0 = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_0, nnnn_0}), nf_bin);
        auto binder_pi = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_pi, nnnn_pi}), nf_bin);
        save_bin_data(binder_0,h5_stats,"binder_0",save_plaintext);
        save_bin_data(binder_pi,h5_stats,"binder_pi",save_plaintext);
    }


    // Local green's functions
    auto gf_im_f = [&](const std::vector<double>& spec, std::complex<double> z, double offset, int norm)->double {
        std::complex<double> d = 0.0;
        for (size_t i=0; i<spec.size(); i++) d+=1./(z - spec[i] + I*offset); 
            return imag(d)/norm;
        };

    auto gf_re_f = [&](const std::vector<double>& spec, std::complex<double> z, double offset, int norm)->double {
            std::complex<double> d = 0.0;
            for (size_t i=0; i<spec.size(); i++) d+=1./(z - spec[i] + I*offset); 
            return real(d)/norm;
            };
    auto dos0_f = [&](const std::vector<double>& spec, std::complex<double> z, double offset, int norm)->double { 
        return -gf_im_f(spec, z, offset, norm)/PI; };

    size_t dos_npts = p["dos_npts"];
    double dos_width = p["dos_width"];
    std::vector<double> grid_real(dos_npts); for (size_t i=0; i<dos_npts; i++) grid_real[i] = -dos_width+2.*dos_width*i/(1.*dos_npts);
    std::vector<double> grid_imag(std::max(int(beta)*10,1024)); for (size_t i=0; i<grid_imag.size(); i++) grid_imag[i] = PI/beta*(2.*i + 1);

    { // gf_matsubara - no errorbars
        triqs::arrays::array<double, 2> gf_im_v(grid_imag.size(),3);
        for (size_t i=0; i<grid_imag.size(); i++) { 
            std::complex<double> z = I*grid_imag[i]; 
            gf_im_v(i,0) = std::imag(z); 
            gf_im_v(i,1) = std::bind(gf_re_f, std::placeholders::_1, z, 0.0, mc.lattice.get_msize())(spectrum); 
            gf_im_v(i,2) = std::bind(gf_im_f, std::placeholders::_1, z, 0.0, mc.lattice.get_msize())(spectrum); 
            };
        h5_write(h5_stats,"gw_imfreq",gf_im_v);
        if (save_plaintext) savetxt("gw_imfreq.dat",gf_im_v);
    };

    { // gf_refreq - no errorbars
        triqs::arrays::array<double, 2> dos_v(grid_real.size(),2), gf_re_v(grid_real.size(),3);
        for (size_t i=0; i<grid_real.size(); i++) { 
            std::complex<double> z = grid_real[i]; 
            dos_v(i,0) = std::real(z); gf_re_v(i,0) = std::real(z); 
            gf_re_v(i,1) = std::bind(gf_re_f, std::placeholders::_1, z, p["dos_offset"], mc.lattice.get_msize())(spectrum); 
            gf_re_v(i,2) = std::bind(gf_re_f, std::placeholders::_1, z, p["dos_offset"], mc.lattice.get_msize())(spectrum); 
            dos_v(i,1) = std::bind(dos0_f, std::placeholders::_1, z, p["dos_offset"], mc.lattice.get_msize())(spectrum);
            };
        h5_write(h5_stats,"dos",dos_v);
        h5_write(h5_stats,"gf_re",gf_re_v);
        if (save_plaintext) { savetxt("dos.dat",dos_v); savetxt("gf_refreq.dat", gf_re_v); };
    };

    if (p["measure_history"])
    {
        const auto &spectrum_history = mc.observables.spectrum_history;
        
        auto dos0_stats = jackknife::accumulate_jackknife(

            std::function<double(std::vector<double>)> 
            (std::bind(dos0_f, std::placeholders::_1, 0.0, p["dos_offset"], mc.lattice.get_msize()))
            ,spectrum_history,maxbin);
        save_binning(dos0_stats,h5_stats,"dos0",save_plaintext);

        size_t dos_bin = estimate_bin(dos0_stats);
        INFO("Using data from bin = " << dos_bin);
        save_bin_data(dos0_stats[dos_bin],h5_stats,"dos0",save_plaintext);

        {
            INFO("\tLocal DOS w errorbars");
            triqs::arrays::array<double, 2> dos_ev(grid_real.size(),3);
            for (size_t i=0; i<grid_real.size(); i++) {
                std::complex<double> z = grid_real[i];
                auto dosz_data = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(dos0_f, std::placeholders::_1, z, p["dos_offset"], mc.lattice.get_msize()))
                    ,spectrum_history,dos_bin);
                dos_ev(i,0) = std::real(z); 
                dos_ev(i,1) = std::get<binning::bin_m::_MEAN>(dosz_data);
                dos_ev(i,2) = std::get<binning::bin_m::_SQERROR>(dosz_data); 
                }
            h5_write(h5_stats,"dos_err",dos_ev);
            if (save_plaintext) savetxt("dos_err.dat",dos_ev);
            }
    }
}

void savetxt (std::string fname, const triqs::arrays::array<double,1>& in)
{
    std::ofstream out(fname);
    out.setf(std::ios::scientific); //out << std::setw(9);
    for (size_t i=0; i<in.shape()[0]; i++){
            out << in(i) << " ";
        };
    out.close();
}

void savetxt (std::string fname, const triqs::arrays::array<double,2>& in)
{
    std::ofstream out(fname);
    out.setf(std::ios::scientific); //out << std::setw(9);
    for (size_t i=0; i<in.shape()[0]; i++){
        for (size_t j=0; j<in.shape()[1]; j++){
            out << in(i,j) << " ";
        };
        out << std::endl;
    };
    out.close();
}

