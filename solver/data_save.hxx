#pragma once

#include <iostream>
#include <fstream>

#include "data_saveload.hpp"

#include "binning.hpp"
#include "jackknife.hpp"

namespace fk {

triqs::utility::parameter_defaults save_defaults() {
  triqs::utility::parameter_defaults pdef;
  pdef.optional
   ("dos_npts", int(100), "Number of points for dos")
   ("dos_width", double(6), "Energy window to save dos")
   ("dos_offset", double(0.05), "dos offset from the real axis")
   ("measure_ipr", bool(false), "Measure inverse participation ratio")
   ;
  return pdef;
 }

void savetxt (std::string fname, const triqs::arrays::array<double,1>& in);
void savetxt (std::string fname, const triqs::arrays::array<double,2>& in);

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

inline void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
}

template <typename MC>
void save_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext)
{
    p.update(save_defaults());
    double beta = p["beta"];
    double Volume = mc.lattice.get_msize();

    print_section("Statistics");
    H5::H5File output(output_file.c_str(),H5F_ACC_TRUNC);
    triqs::h5::group top(output);
    //===== save parameters ===== //
    h5_write(top, "parameters", p);

    //===== save direct measures ===== //
    top.create_group("mc_data");
    auto h5_mc_data = top.open_group("mc_data");
    if (mc.observables.energies.size()) h5_write(h5_mc_data,"energies", mc.observables.energies);
    if (mc.observables.d2energies.size()) h5_write(h5_mc_data,"d2energies", mc.observables.d2energies);
    if (mc.observables.nf0.size()) h5_write(h5_mc_data,"nf0", mc.observables.nf0);
    if (mc.observables.nfpi.size()) h5_write(h5_mc_data,"nfpi", mc.observables.nfpi);

    std::vector<double> spectrum(mc.observables.spectrum.size());
    if (spectrum.size()) { 
        std::copy(mc.observables.spectrum.data(), mc.observables.spectrum.data()+spectrum.size(), spectrum.begin());
        h5_write(h5_mc_data,"spectrum", spectrum);
        };

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

    // Inverse participation ratio
    if (p["measure_ipr"] && p["measure_history"]) {
        std::cout << "Inverse participation ratio" << std::endl;
        triqs::arrays::array<double, 2> t_ipr_history(mc.observables.ipr_history.size(), mc.observables.ipr_history[0].size());
        for (int i=0; i<mc.observables.ipr_history.size(); i++)
            for (int j=0; j< mc.observables.ipr_history[0].size(); j++)
                t_ipr_history(i,j) =  mc.observables.ipr_history[i][j];
        h5_write(h5_mc_data,"ipr_history", t_ipr_history);
        };


    //===== save statistics ===== //
    top.create_group("stats");
    auto h5_stats = top.open_group("stats");
    top.create_group("binning");
    auto h5_binning = top.open_group("binning");

    int maxbin = std::min(15,std::max(int(std::log(double(p["n_cycles"])/16)/std::log(2.)-1),1));

    if (mc.observables.energies.size()) { 
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
    };

    if (mc.observables.nf0.size() && mc.observables.nfpi.size()) { 
          /* Save nf(q=0), nf(q=pi) */
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

    if (mc.observables.energies.size()){ // gf_matsubara - no errorbars
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

    if (mc.observables.energies.size()) { // gf_refreq - no errorbars
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
    { // dos(w=0)
        const auto &spectrum_history = mc.observables.spectrum_history;
        
        auto dos0_stats = jackknife::accumulate_jackknife(

            std::function<double(std::vector<double>)> 
            (std::bind(dos0_f, std::placeholders::_1, 0.0, p["dos_offset"], mc.lattice.get_msize()))
            ,spectrum_history,maxbin);
        save_binning(dos0_stats,h5_stats,"dos0",save_plaintext);

        size_t dos_bin = estimate_bin(dos0_stats);
        INFO("Using data from bin = " << dos_bin);
        save_bin_data(dos0_stats[dos_bin],h5_stats,"dos0",save_plaintext);
        // dos(w)
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
    
    
    // Inverse participation ratio
    if (p["measure_ipr"]) {
            auto ipr_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
                std::complex<double> nom = 0.0, denom = 0.0;
                for (size_t i=0; i<Volume; i++) {
                    denom+=1./(z - ipr_spec[i] + I*offset); 
                    nom+=1./(z - ipr_spec[i] + I*offset)*ipr_spec[i+Volume]; 
                    };
                return imag(nom)/imag(denom);
                };

            const auto& ipr_vals = mc.observables.ipr_history;

            typedef std::vector<double>::const_iterator iter_t;
            assert(Volume == ipr_vals.size());
            std::vector<std::pair<iter_t,iter_t>> ipr_and_spectrum(2*Volume); // create a vector of pair of 2*Volume size
            for (size_t i=0; i<Volume; ++i) { 
                ipr_and_spectrum[i]=std::make_pair(spectrum_history[i].begin(),spectrum_history[i].end());
                ipr_and_spectrum[i+Volume]=std::make_pair(ipr_vals[i].begin(),ipr_vals[i].end());
            }

            { // save ipr at w=0
                auto ipr0_stats = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(ipr_f, std::placeholders::_1, 0.0, p["dos_offset"]))
                    ,ipr_and_spectrum,dos_bin);

                //save_binning(ipr0_stats,h5_stats,"ipr0",save_plaintext);
                save_bin_data(ipr0_stats,h5_stats,"ipr0",save_plaintext);
            }

            triqs::arrays::array<double, 2> ipr_ev(grid_real.size(),3);
            for (size_t i=0; i<grid_real.size(); i++) {
                std::complex<double> z = grid_real[i];
                auto ipr_data = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(ipr_f, std::placeholders::_1, z, p["dos_offset"]))
                    ,ipr_and_spectrum,dos_bin);
                ipr_ev(i,0) = std::real(z); 
                ipr_ev(i,1) = std::get<binning::bin_m::_MEAN>(ipr_data);
                ipr_ev(i,2) = std::get<binning::bin_m::_SQERROR>(ipr_data); 

                };

            h5_write(h5_stats,"ipr_err",ipr_ev);
            if (save_plaintext) savetxt("ipr_err.dat",ipr_ev);

        } // end measure_ipr
    } // end measure_history
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

} // end of namespace fk
