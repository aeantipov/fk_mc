#pragma once

#include <iostream>
#include <fstream>


#include "binning.hpp"
#include "jackknife.hpp"

namespace fk {

// save data from solver to hdf5 file 
template <typename MC>
void save_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false, std::vector<double> wgrid_cond = {0.0});

// save g(w,r) and g(w,k)
template <typename LATTICE>
void save_gwr(observables_t const& obs, LATTICE const& lattice, triqs::utility::parameters const& p, std::vector<std::complex<double>> wgrid);

triqs::utility::parameter_defaults save_defaults() {
  triqs::utility::parameter_defaults pdef;
  pdef
   .optional("dos_npts", int(100), "Number of points for dos")
   .optional("dos_width", double(6), "Energy window to save dos")
   .optional("measure_ipr", bool(false), "Measure inverse participation ratio")
   .optional("dos_offset", double(0.05), "dos offset from the real axis")
   ;
  return pdef;
 }

void savetxt (std::string fname, const triqs::arrays::array<double,1>& in);
void savetxt (std::string fname, const triqs::arrays::array<double,2>& in);

/// Estimate the bin, at which the error bar is saturated.
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



void save_bin_data(const binning::bin_stats_t& data, triqs::h5::group& h5_group, std::string name, int bin, bool save_plaintext = false)
{
    std::array<double, 4> tmp (
        {{double(std::get<binning::_SIZE>(data)), std::get<binning::_MEAN>(data), std::get<binning::_DISP>(data), std::get<binning::_SQERROR>(data) }});
    tqa::array<double, 1> data_arr(4);
    std::copy(tmp.begin(),tmp.end(),data_arr(tqa::range()).begin());

    std::cout << name << " = " << std::get<binning::_MEAN>(data) << " +/- " << std::get<binning::_SQERROR>(data) 
              << " (" << std::get<binning::_SIZE>(data) << ") samples [" << bin << "]." << std::endl;
    h5_write(h5_group, name, data_arr);
    if (save_plaintext) savetxt(name+"_error.dat",data_arr);
}
void save_binning(const binning::bin_data_t& binning_data, triqs::h5::group& h5_group, triqs::h5::group& h5_stats, std::string name, bool save_plaintext = false)
{
    auto cor_lens = binning::calc_cor_length(binning_data);
    tqa::array<double, 2> data_arr(binning_data.size(),5);
    //std::ofstream out; out.setf(std::ios::scientific); //out << std::setw(9);
    //if (save_plaintext) out.open(name+"_binning.dat",std::ios::out);
    for (size_t i=0; i<binning_data.size(); i++){
        auto e = binning_data[i]; double c = cor_lens[i];
        std::array<double, 5> t ({{double(std::get<binning::_SIZE>(e)), std::get<binning::_MEAN>(e), std::get<binning::_DISP>(e), std::get<binning::_SQERROR>(e), c }});
        std::copy(t.begin(),t.end(),data_arr(i,tqa::range()).begin());
        //if (save_plaintext) out << i << " " << std::get<binning::_SIZE>(e) << "  " << std::get<binning::_MEAN>(e) 
            //<< "  " << std::get<binning::_DISP>(e) << "  " << std::get<binning::_SQERROR>(e) << "  " << c << "\n";
       };
    if (save_plaintext) savetxt(name+"_binning.dat",data_arr);
    h5_write(h5_group,name, data_arr);
            
    int data_bin = estimate_bin(binning_data);
    save_bin_data(binning_data[data_bin],h5_stats,name,data_bin,save_plaintext);
}

inline void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
}

template <typename MC>
void save_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext, std::vector<double> wgrid_cond)
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
    if (mc.observables.stiffness.size()) h5_write(h5_mc_data,"stiffness", mc.observables.stiffness);

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

    // Conductivity
    if (p["measure_stiffness"]) {
        std::cout << "Conductivity" << std::endl;
        triqs::arrays::array<double, 2> t_cond_history(mc.observables.cond_history.size(), mc.observables.cond_history[0].size());
        for (int i=0; i<mc.observables.cond_history.size(); i++)
            for (int j=0; j< mc.observables.cond_history[0].size(); j++)
                t_cond_history(i,j) =  mc.observables.cond_history[i][j];
        h5_write(h5_mc_data,"cond_history", t_cond_history);
        };

    if (p["measure_eigenfunctions"]) {
        std::cout << "Eigenfunctions" << std::endl;
        const auto& eig_hist = mc.observables.eigenfunctions_history;
        triqs::arrays::array<double, 3> t_eig_history(eig_hist.size(), eig_hist[0].rows(), eig_hist[0].cols() );
        for (int i=0; i<eig_hist.size(); i++)
            for (int j=0; j< eig_hist[0].rows(); j++) 
                for (int k=0; k< eig_hist[0].cols(); k++) 
                    t_eig_history(i,j,k) =  eig_hist[i](j,k);
        h5_write(h5_mc_data,"eig_history", t_eig_history);
        }


    //===== save statistics ===== //
    top.create_group("stats");
    auto h5_stats = top.open_group("stats");
    top.create_group("binning");
    auto h5_binning = top.open_group("binning");

    int ncycles = mc.observables.nfpi.size();
    int maxbin = std::min(15,std::max(int(std::log(double(ncycles)/16)/std::log(2.)-1),1));

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
            save_binning(energy_binning_data,h5_binning,h5_stats,"energy",save_plaintext);
            auto d2energy_binning_data = binning::accumulate_binning(d2energies.rbegin(),d2energies.rend(), maxbin); 
            save_binning(d2energy_binning_data,h5_binning,h5_stats,"d2energy",save_plaintext);
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
            save_binning(cv_stats,h5_binning,h5_stats,"cv",save_plaintext);
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
        save_binning( nf0_stats,h5_binning,h5_stats,"nf_0",save_plaintext);
        auto nfpi_stats = binning::accumulate_binning(nfpi.rbegin(), nfpi.rend(), maxbin);
        save_binning( nfpi_stats,h5_binning,h5_stats,"nf_pi",save_plaintext);
        auto fsusc0_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nf0,nn_0}),maxbin);
        save_binning(fsusc0_stats,h5_binning,h5_stats,"fsusc_0",save_plaintext);
        auto fsuscpi_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nfpi,nn_pi}),maxbin);
        save_binning(fsuscpi_stats,h5_binning,h5_stats,"fsusc_pi",save_plaintext);

    
        /* Binder cumulant (q=0, pi). */
        std::vector<double> nnnn_0(nf0.size()), nnnn_pi(nfpi.size());
        std::transform(nf0.begin(), nf0.end(), nnnn_0.begin(), [](double x){return x*x*x*x;});
        std::transform(nfpi.begin(), nfpi.end(), nnnn_pi.begin(), [](double x){return x*x*x*x;});
        std::function<double(double,double)> binder_f = [](double x2, double x4){return 1. - x4/3./x2/x2;};
        auto nf_bin = estimate_bin(fsusc0_stats);
        auto binder_0 = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_0, nnnn_0}), nf_bin);
        auto binder_pi = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_pi, nnnn_pi}), nf_bin);
        save_bin_data(binder_0,h5_stats,"binder_0",save_plaintext);
        save_bin_data(binder_pi,h5_stats,"binder_pi",save_plaintext);
    }

    if (p["measure_stiffness"]) { // Stiffness
            INFO("Stiffness");
            const auto& stiffness = mc.observables.stiffness;
            size_t size = stiffness.size();
            auto stiffness_data = binning::accumulate_binning(stiffness.rbegin(), stiffness.rend(), maxbin); 
            save_binning(stiffness_data,h5_binning,h5_stats,"stiffness",save_plaintext);

            // conductivity
            std::vector<std::vector<double>> const& cond_history = mc.observables.cond_history;
            std::vector<double> const& wgrid = wgrid_cond; 

            // cond(w)
            auto cond_stats0 = binning::accumulate_binning(cond_history[wgrid.size()/2].begin(), cond_history[wgrid.size()/2].end(), maxbin);
            int cond_bin = estimate_bin(cond_stats0); 
            {
                INFO("Saving w*conductivity (w)");
                triqs::arrays::array<double, 2> wcond_ev(wgrid.size(),3);
                for (size_t i=0; i<wgrid.size(); i++) {
                    std::complex<double> z = wgrid[i];
                    auto cond_binning_data = binning::bin(cond_history[i].begin(), cond_history[i].end(), cond_bin);
                    wcond_ev(i,0) = std::real(z); 
                    wcond_ev(i,1) = std::get<binning::bin_m::_MEAN>(cond_binning_data);
                    wcond_ev(i,2) = std::get<binning::bin_m::_SQERROR>(cond_binning_data); 
                    }
                h5_write(h5_stats,"wcond_err",wcond_ev);
                if (save_plaintext) savetxt("wcond_err.dat",wcond_ev);

                int index_zero = (wgrid.size() - 1)/2;
                std::cout << "wgrid center = " << wgrid[index_zero] << std::endl;
                //if (std::abs(wgrid[index_zero]) > std::numeric_limits<double>::epsilon()) index_zero = std::distance(wgrid.begin(), wgrid
                double cond_left = -wcond_ev(index_zero - 1,1)/wgrid[index_zero - 1];
                double cond_right = -wcond_ev(index_zero + 1,1)/wgrid[index_zero + 1];
                double cond0 = (cond_left + cond_right)/2;
                double cond0_err = wcond_ev(index_zero + 1, 2)/wgrid[index_zero + 1];

                // save dc conductivity
                triqs::arrays::array<double, 1> t_cond0(4); 
                t_cond0(binning::bin_m::_SIZE) = std::get<binning::bin_m::_SIZE>(cond_stats0[cond_bin]); 
                t_cond0(binning::bin_m::_MEAN) = cond0; 
                t_cond0(binning::bin_m::_DISP) = t_cond0(binning::bin_m::_SIZE)*cond0_err*cond0_err; 
                t_cond0(binning::bin_m::_SQERROR)=cond0_err;
                h5_write(h5_stats,"cond0",t_cond0);
                if (save_plaintext) savetxt("cond0.dat",t_cond0);

                triqs::arrays::array<double, 2> cond_ev(wcond_ev), cond_ev_subtract(wcond_ev);
                for (size_t i=0; i<wgrid.size(); i++) { 
                    cond_ev(i,1) = -wcond_ev(i,1) / wgrid[i];
                    cond_ev(i,2) = -wcond_ev(i,2) / wgrid[i];
                    cond_ev_subtract(i,1) = cond_ev(i,1) - cond0;
                    cond_ev_subtract(i,2) = cond_ev(i,2);
                    }
                
                h5_write(h5_stats,"cond_err",cond_ev);
                if (save_plaintext) savetxt("cond_err.dat",cond_ev);
                h5_write(h5_stats,"cond_dynamic",cond_ev_subtract);
                if (save_plaintext) savetxt("cond_dynamic.dat",cond_ev_subtract);

                }
        };
        
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


    // G(w,k)
    if (p["measure_eigenfunctions"]) { save_gwr(mc.observables, mc.lattice, p, {std::complex<double>(0,0)}); } 
    
    if (p["measure_history"])
    { // dos(w=0)
        const auto &spectrum_history = mc.observables.spectrum_history;
        
        auto dos0_stats = jackknife::accumulate_jackknife(

            std::function<double(std::vector<double>)> 
            (std::bind(dos0_f, std::placeholders::_1, 0.0, p["dos_offset"], mc.lattice.get_msize()))
            ,spectrum_history,maxbin);
        save_binning(dos0_stats,h5_binning,h5_stats,"dos0",save_plaintext);
        size_t dos_bin = estimate_bin(dos0_stats);
        // dos(w)
        {
            INFO("Saving local DOS w errorbars");
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
                    // as norm4 is measured -> take the power of 4 to extract ipr
                    nom+=1./(z - ipr_spec[i] + I*offset)*std::pow(ipr_spec[i+Volume],4); 
                    };
                return imag(nom)/imag(denom);
                };

            auto ipr_thermal_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
                double out = 0.0;
                double ipr_state, state_weight;
                for (size_t i=0; i<Volume; i++) {
                    ipr_state = std::pow(ipr_spec[i+Volume],4);
                    state_weight = beta / (1. + std::exp(beta*ipr_spec[i])) / (1. + std::exp(-beta*ipr_spec[i]));
                    out += ipr_state * state_weight / Volume;  
                    };
                return out; 
                };
            auto dos_thermal_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
                double out = 0.0;
                double energy_state, state_weight;
                for (size_t i=0; i<Volume; i++) {
                    energy_state = ipr_spec[i];
                    state_weight = beta / (1. + std::exp(beta*ipr_spec[i])) / (1. + std::exp(-beta*ipr_spec[i]));
                    out += state_weight / Volume;  
                    };
                return out; 
                };


            const auto& ipr_vals = mc.observables.ipr_history;

            typedef std::vector<double>::const_iterator iter_t;
            assert(Volume == ipr_vals.size());
            std::vector<std::pair<iter_t,iter_t>> ipr_and_spectrum(2*Volume); // create a vector of pair of 2*Volume size
            for (size_t i=0; i<Volume; ++i) { 
                ipr_and_spectrum[i]=std::make_pair(spectrum_history[i].begin(),spectrum_history[i].end());
                ipr_and_spectrum[i+Volume]=std::make_pair(ipr_vals[i].begin(),ipr_vals[i].end());
            }

            typename binning::bin_data_t ipr0_binning(maxbin);
            for (int i=0; i<maxbin; i++) 
                { // save ipr at w=0
                    auto ipr0_stats = jackknife::jack(
                        std::function<double(std::vector<double>)>( 
                        std::bind(ipr_f, std::placeholders::_1, 0.0, p["dos_offset"]))
                        ,ipr_and_spectrum,i);

                    ipr0_binning[i] = ipr0_stats;
                }
            save_binning(ipr0_binning,h5_binning,h5_stats,"ipr0",save_plaintext);

            // ipr - thermal
            typename binning::bin_data_t ipr_thermal_binning(maxbin);
            for (int i=0; i<maxbin; i++) 
                { // save ipr at w=0
                    auto ipr_th_stats = jackknife::jack(
                        std::function<double(std::vector<double>)>( 
                        std::bind(ipr_thermal_f, std::placeholders::_1, 0.0, p["dos_offset"]))
                        ,ipr_and_spectrum,i);

                    ipr_thermal_binning[i] = ipr_th_stats;
                }
            save_binning(ipr_thermal_binning,h5_binning,h5_stats,"ipr_thermal",save_plaintext);

            // dos - thermal
            typename binning::bin_data_t dos_thermal_binning(maxbin);
            for (int i=0; i<maxbin; i++) 
                { // save ipr at w=0
                    auto dos_th_stats = jackknife::jack(
                        std::function<double(std::vector<double>)>( 
                        std::bind(dos_thermal_f, std::placeholders::_1, 0.0, p["dos_offset"]))
                        ,ipr_and_spectrum,i);

                    dos_thermal_binning[i] = dos_th_stats;
                }
            save_binning(dos_thermal_binning,h5_binning,h5_stats,"dos_thermal",save_plaintext);


            auto ipr0_bin=estimate_bin(ipr0_binning);
            triqs::arrays::array<double, 2> ipr_ev(grid_real.size(),3);
            for (size_t i=0; i<grid_real.size(); i++) {
                std::complex<double> z = grid_real[i];
                auto ipr_data = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(ipr_f, std::placeholders::_1, z, p["dos_offset"]))
                    ,ipr_and_spectrum,ipr0_bin);
                ipr_ev(i,0) = std::real(z); 
                ipr_ev(i,1) = std::get<binning::bin_m::_MEAN>(ipr_data);
                ipr_ev(i,2) = std::get<binning::bin_m::_SQERROR>(ipr_data); 

                };

            h5_write(h5_stats,"ipr_err",ipr_ev);
            if (save_plaintext) savetxt("ipr_err.dat",ipr_ev);

        } // end measure_ipr

    // f-electron correlation functions
    // assuming x <-> y symmetry
    const auto& fhistory = mc.observables.focc_history;
    const auto dims = mc.lattice.dims;

    int nf_bin = estimate_bin(binning::accumulate_binning(fhistory[0].rbegin(), fhistory[0].rend(), maxbin));

    std::vector<double> nf_mean(Volume, 0);
    for (int i = 0; i < Volume; i++) { 
        auto nf_stats = binning::bin(fhistory[i].rbegin(), fhistory[i].rend(), nf_bin); 
        nf_mean[i] = std::get<binning::_MEAN>(nf_stats);
        }
        

    auto fcorrel_f = [&](const std::vector<double> focc_history, int l)->double { 
        double out = 0.0;
        for (size_t i=0; i<Volume; i++) {
            double nf_i = focc_history[i];
            double nf_i_mean = nf_mean[i];
            auto current_pos = mc.lattice.index_to_pos(i); 
            for (int d = 0; d < mc.lattice.Ndim; d++) { 
                auto pos_l(current_pos), pos_r(current_pos);
                pos_l[d]=(current_pos[d] - l + dims[d]) % dims[d];
                pos_r[d]=(current_pos[d] + l + dims[d]) % dims[d];
                size_t index_l = mc.lattice.pos_to_index(pos_l);
                size_t index_r = mc.lattice.pos_to_index(pos_r);

                double nf_l = focc_history[index_l]; 
                double nf_r = focc_history[index_r]; 
                double nf_l_mean = nf_mean[index_l];
                double nf_r_mean = nf_mean[index_r];
                out += (nf_i - nf_i_mean) * (nf_l - nf_l_mean);
                out += (nf_i - nf_i_mean) * (nf_r - nf_r_mean);
                };
            }
            return out / Volume / (2.0 * mc.lattice.Ndim); 
        };

        auto fcorrel0_stats = jackknife::accumulate_jackknife(std::function<double(std::vector<double>)>(std::bind(fcorrel_f, std::placeholders::_1, 0)),fhistory,maxbin);
        nf_bin = estimate_bin(fcorrel0_stats);
        double fcorrel0_mean = std::get<binning::_MEAN>(fcorrel0_stats[nf_bin]);
        double fcorrel0_error = std::get<binning::_SQERROR>(fcorrel0_stats[nf_bin]);

        triqs::arrays::array<double, 2> fcorrel_out(mc.lattice.dims[0] / 2, 5);
        for (int l = 0; l < mc.lattice.dims[0] / 2; l++) { 
            auto fcorrel_stats = jackknife::jack(std::function<double(std::vector<double>)>(std::bind(fcorrel_f, std::placeholders::_1, l)),fhistory,nf_bin);
            save_bin_data(fcorrel_stats,h5_stats,"fcorrel_" + std::to_string(l),save_plaintext);
            double fcorrel_mean = std::get<binning::_MEAN>(fcorrel_stats);
            double fcorrel_error = std::get<binning::_SQERROR>(fcorrel_stats);
            fcorrel_out(l,0) = l;
            fcorrel_out(l,1) = fcorrel_mean;
            fcorrel_out(l,2) = fcorrel_error;
            fcorrel_out(l,3) = fcorrel_mean / fcorrel0_mean;
            fcorrel_out(l,4) = std::sqrt(std::pow(fcorrel_error / fcorrel0_mean, 2) + std::pow(fcorrel_mean / (fcorrel0_mean * fcorrel0_mean) * fcorrel0_error, 2));
            }

        h5_write(h5_stats,"fcorrel",fcorrel_out);
        if (save_plaintext) savetxt("fcorrel.dat",fcorrel_out);

    } // end measure_history
}


// Save G(w,r) and G(w,k) to plaintext files
// Warning: only works for 2d
template <typename LATTICE>
void save_gwr(observables_t const& obs, LATTICE const& lattice, triqs::utility::parameters const& p, std::vector<std::complex<double>> wgrid) 
{ 
    typedef observables_t::dense_m dense_m;
    const std::vector<dense_m>& eigs = obs.eigenfunctions_history;
    std::cout << "Saving G(w,k)" << std::endl;

    // create a grid for gw
    std::vector<double> grid_real2; 
    grid_real2.push_back(0);

    dense_m gwr_re(eigs[0].rows(), eigs[0].cols()); 
    dense_m gwr_im(eigs[0].rows(), eigs[0].cols()); 
    //Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> evals(eigs[0].rows());
    //Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ident_v(eigs[0].rows());
    dense_m evals(eigs[0].rows(), eigs[0].cols()); evals.setZero();
    dense_m ident_v(eigs[0].rows(), eigs[0].cols());
    ident_v.setIdentity();

    double xi = p["dos_offset"];
    dense_m xi_m = ident_v * xi;

    for (std::complex<double> w : wgrid) { 
        std::string wstring = std::to_string(float(w.real())) + "_" + std::to_string(float(w.imag()));
        gwr_re.setZero();
        gwr_im.setZero();
        for (int m = 0; m < eigs.size(); ++m) { 
            for (int i = 0; i < eigs[0].rows(); ++i) { evals.diagonal()(i) = obs.spectrum_history[i][m]; }
            auto wminuseps = ident_v * w.real() - evals;
            std::cout << (wminuseps * wminuseps + (xi_m + w.imag() * ident_v)*(xi_m + w.imag() * ident_v)).inverse() << std::endl;// wminuseps*wminuseps + xi_m*xi_m << std::endl;
            // eigenvalues are in columns
            std::cout << "norm = " << eigs[0].row(0) * eigs[0].row(0).transpose() << std::endl;
            std::cout << "norm2 = " << eigs[0].col(0).transpose() * eigs[0].col(0) << std::endl;
            gwr_im -= eigs[m].transpose() * xi * (wminuseps * wminuseps + xi_m*xi_m).inverse() * eigs[m] / eigs.size();
            gwr_re += eigs[m].transpose() * wminuseps * (wminuseps * wminuseps + xi_m*xi_m).inverse() * eigs[m] / eigs.size();
            }

        std::ofstream gwr_re_str("gr_full_w"+wstring+"_re.dat");
        std::ofstream gwr_im_str("gr_full_w"+wstring+"_im.dat");
        gwr_re_str << gwr_re << std::endl;
        gwr_im_str << gwr_im << std::endl;
        gwr_re_str.close();
        gwr_im_str.close();

        // now gwr_re, gwr_im contain Gw(r1,r2)
        // let's now convert it to Gw(r1 - r2)
        auto dims = lattice.dims;
        dense_m gwr_im2 = dense_m::Zero(dims[0], dims[1]);
        dense_m gwr_re2 = dense_m::Zero(dims[0], dims[1]);
        for (int i=0; i<gwr_im.rows(); ++i) { 
            auto pos_i = lattice.index_to_pos(i);
            for (int j =0; j < gwr_im.cols(); ++j) { 
                auto pos_j = lattice.index_to_pos(j);
                //std::cout << " i = " << pos_i << " j = " << pos_j << std::endl;
        
                int p0 = (dims[0] + pos_j[0] - pos_i[0]) % dims[0];
                int p1 = (dims[1] + pos_j[1] - pos_i[1]) % dims[1];

                gwr_im2(p0,p1) += gwr_im(i, j) / gwr_im.cols();
                gwr_re2(p0,p1) += gwr_re(i, j) / gwr_im.cols();
                }
            }

        gwr_re_str.open("gr_w"+wstring+"_re.dat");
        gwr_im_str.open("gr_w"+wstring+"_im.dat");
        gwr_re_str << gwr_re2 << std::endl;
        gwr_im_str << gwr_im2 << std::endl;
        gwr_re_str.close();
        gwr_im_str.close();
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

} // end of namespace fk
