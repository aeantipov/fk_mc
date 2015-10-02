#pragma once

#include "data_save.hpp"

namespace fk {

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

template <typename MC>
void data_saver<MC>::save_all(std::vector<double> wgrid_cond) 
{ 
    double beta = p_["beta"];
    std::string output_file = p_["output_file"];
    bool save_plaintext = p_["save_plaintext"];

    print_section("Statistics");
    H5::H5File output(output_file.c_str(),H5F_ACC_TRUNC);
    top_ = triqs::h5::group(output);
    // save parameters to "/parameters" dataset
    h5_write(top_, "parameters", p_);

    // save measurement to "/mc_data" group
    top_.create_group("mc_data");
    h5_mc_data_ = top_.open_group("mc_data");
    save_observables();
    
    // save final statistics to "/stats" group
    top_.create_group("stats");
    h5_stats_ = top_.open_group("stats");
    // save binned data to "/binning" group
    top_.create_group("binning");
    h5_binning_ = top_.open_group("binning");

    std::vector<double> const& spectrum = observables_.spectrum;
    // save energy and the specific heat
    if (observables_.energies.size()) { this->save_energy(); }
    // save f-electron susc and binder cumulant
    if (observables_.nf0.size() && observables_.nfpi.size()) { this->save_fstats(); }
    // f-electron correlation functions
    if (p_["measure_history"]) { this->save_fcorrel(); }

    size_t dos_npts = p_["dos_npts"];
    double dos_width = p_["dos_width"];
    std::vector<double> grid_real(dos_npts); for (size_t i=0; i<dos_npts; i++) grid_real[i] = -dos_width+2.*dos_width*i/(1.*dos_npts);

    // Save glocal
    this->save_glocal(grid_real);
    // Inverse participation ratio
    if (p_["measure_ipr"]) { this->save_ipr(grid_real); }
    // Conductivity and Drude weight
    if (p_["measure_stiffness"]) { this->save_conductivity(wgrid_cond); };
    // G(w,k)
    if (p_["measure_eigenfunctions"]) { save_gwr({std::complex<double>(0,0)}); } 
    
} // end measure_history

template <typename MC>
void data_saver<MC>::save_fstats()
{
    bool save_plaintext = p_["save_plaintext"];
    if (observables_.nf0.size() && observables_.nfpi.size()) { 
          /* Save nf(q=0), nf(q=pi) */
        const std::vector<double>& nf0 = observables_.nf0;
        const std::vector<double>& nfpi = observables_.nfpi;
        std::vector<double> nn_0(nf0.size()), nn_pi(nfpi.size());
        std::transform(nf0.begin(), nf0.end(), nn_0.begin(), [](double x){return x*x;});
        std::transform(nfpi.begin(), nfpi.end(), nn_pi.begin(), [](double x){return x*x;});
        std::function<double(double,double)> disp_f = [](double x, double x2){return x2 - x*x;};
        
        auto nf0_stats = binning::accumulate_binning(nf0.rbegin(), nf0.rend(), max_bin_);
        save_binning( nf0_stats,h5_binning_,h5_stats_,"nf_0",save_plaintext);
        auto nfpi_stats = binning::accumulate_binning(nfpi.rbegin(), nfpi.rend(), max_bin_);
        save_binning( nfpi_stats,h5_binning_,h5_stats_,"nf_pi",save_plaintext);
        auto fsusc0_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nf0,nn_0}),max_bin_);
        save_binning(fsusc0_stats,h5_binning_,h5_stats_,"fsusc_0",save_plaintext);
        auto fsuscpi_stats = jackknife::accumulate_jackknife(disp_f,std::vector<std::vector<double>>({nfpi,nn_pi}),max_bin_);
        save_binning(fsuscpi_stats,h5_binning_,h5_stats_,"fsusc_pi",save_plaintext);

    
        /* Binder cumulant (q=0, pi). */
        std::vector<double> nnnn_0(nf0.size()), nnnn_pi(nfpi.size());
        std::transform(nf0.begin(), nf0.end(), nnnn_0.begin(), [](double x){return x*x*x*x;});
        std::transform(nfpi.begin(), nfpi.end(), nnnn_pi.begin(), [](double x){return x*x*x*x;});
        std::function<double(double,double)> binder_f = [](double x2, double x4){return 1. - x4/3./x2/x2;};
        auto nf_bin = estimate_bin(fsusc0_stats);
        auto binder_0 = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_0, nnnn_0}), nf_bin);
        auto binder_pi = jackknife::jack(binder_f, std::vector<std::vector<double>>({nn_pi, nnnn_pi}), nf_bin);
        save_bin_data(binder_0,h5_stats_,"binder_0",save_plaintext);
        save_bin_data(binder_pi,h5_stats_,"binder_pi",save_plaintext);
    }
}

template <typename MC>
void data_saver<MC>::save_glocal(std::vector<double> grid_real)
{
    bool save_plaintext = p_["save_plaintext"];
    double beta = mc_.config.params().beta;
    std::vector<double> const& spectrum = observables_.spectrum;
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

    std::vector<double> grid_imag(std::max(int(beta)*10,1024)); for (size_t i=0; i<grid_imag.size(); i++) grid_imag[i] = PI/beta*(2.*i + 1);

    if (observables_.energies.size()){ // gf_matsubara - no errorbars
        triqs::arrays::array<double, 2> gf_im_v(grid_imag.size(),3);
        for (size_t i=0; i<grid_imag.size(); i++) { 
            std::complex<double> z = I*grid_imag[i]; 
            gf_im_v(i,0) = std::imag(z); 
            gf_im_v(i,1) = std::bind(gf_re_f, std::placeholders::_1, z, 0.0, lattice_.get_msize())(spectrum); 
            gf_im_v(i,2) = std::bind(gf_im_f, std::placeholders::_1, z, 0.0, lattice_.get_msize())(spectrum); 
            };
        h5_write(h5_stats_,"gw_imfreq",gf_im_v);
        if (save_plaintext) savetxt("gw_imfreq.dat",gf_im_v);
    };

    if (observables_.energies.size()) { // gf_refreq - no errorbars
        triqs::arrays::array<double, 2> dos_v(grid_real.size(),2), gf_re_v(grid_real.size(),3);
        for (size_t i=0; i<grid_real.size(); i++) { 
            std::complex<double> z = grid_real[i]; 
            dos_v(i,0) = std::real(z); gf_re_v(i,0) = std::real(z); 
            gf_re_v(i,1) = std::bind(gf_re_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum); 
            gf_re_v(i,2) = std::bind(gf_re_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum); 
            dos_v(i,1) = std::bind(dos0_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum);
            };
        h5_write(h5_stats_,"dos",dos_v);
        h5_write(h5_stats_,"gf_re",gf_re_v);
        if (save_plaintext) { savetxt("dos.dat",dos_v); savetxt("gf_refreq.dat", gf_re_v); };
    };

    /// Save glocal with error-bars
    if (p_["measure_history"])
    { // dos(w=0)
        const auto &spectrum_history = observables_.spectrum_history;
        
        auto dos0_stats = jackknife::accumulate_jackknife(

            std::function<double(std::vector<double>)> 
            (std::bind(dos0_f, std::placeholders::_1, 0.0, p_["dos_offset"], lattice_.get_msize()))
            ,spectrum_history,max_bin_);
        save_binning(dos0_stats,h5_binning_,h5_stats_,"dos0",save_plaintext);
        size_t dos_bin = estimate_bin(dos0_stats);
        // dos(w)
        {
            INFO("Saving local DOS w errorbars");
            triqs::arrays::array<double, 2> dos_ev(grid_real.size(),3);
            for (size_t i=0; i<grid_real.size(); i++) {
                std::complex<double> z = grid_real[i];
                auto dosz_data = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(dos0_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize()))
                    ,spectrum_history,dos_bin);
                dos_ev(i,0) = std::real(z); 
                dos_ev(i,1) = std::get<binning::bin_m::_MEAN>(dosz_data);
                dos_ev(i,2) = std::get<binning::bin_m::_SQERROR>(dosz_data); 
                }
            h5_write(h5_stats_,"dos_err",dos_ev);
            if (save_plaintext) savetxt("dos_err.dat",dos_ev);
            }
      }
}

template <typename MC>
void data_saver<MC>::save_fcorrel()
{
    bool save_plaintext = p_["save_plaintext"];
    // f-electron correlation functions
    // assuming x <-> y symmetry
    const auto& fhistory = observables_.focc_history;
    if (fhistory.size()==0) return;
    const auto dims = lattice_.dims;

    int nf_bin = estimate_bin(binning::accumulate_binning(fhistory[0].rbegin(), fhistory[0].rend(), max_bin_));

    std::vector<double> nf_mean(volume_, 0);
    for (int i = 0; i < volume_; i++) { 
        auto nf_stats = binning::bin(fhistory[i].rbegin(), fhistory[i].rend(), nf_bin); 
        nf_mean[i] = std::get<binning::_MEAN>(nf_stats);
        }
        

    auto fcorrel_f = [&](const std::vector<double> focc_history, int l)->double { 
        double out = 0.0;
        for (size_t i=0; i<volume_; i++) {
            double nf_i = focc_history[i];
            double nf_i_mean = nf_mean[i];
            auto current_pos = lattice_.index_to_pos(i); 
            for (int d = 0; d < lattice_.Ndim; d++) { 
                auto pos_l(current_pos), pos_r(current_pos);
                pos_l[d]=(current_pos[d] - l + dims[d]) % dims[d];
                pos_r[d]=(current_pos[d] + l + dims[d]) % dims[d];
                size_t index_l = lattice_.pos_to_index(pos_l);
                size_t index_r = lattice_.pos_to_index(pos_r);

                double nf_l = focc_history[index_l]; 
                double nf_r = focc_history[index_r]; 
                double nf_l_mean = nf_mean[index_l];
                double nf_r_mean = nf_mean[index_r];
                out += (nf_i - nf_i_mean) * (nf_l - nf_l_mean);
                out += (nf_i - nf_i_mean) * (nf_r - nf_r_mean);
                };
            }
            return out / volume_ / (2.0 * lattice_.Ndim); 
        };

        auto fcorrel0_stats = jackknife::accumulate_jackknife(std::function<double(std::vector<double>)>(std::bind(fcorrel_f, std::placeholders::_1, 0)),fhistory,max_bin_);
        nf_bin = estimate_bin(fcorrel0_stats);
        double fcorrel0_mean = std::get<binning::_MEAN>(fcorrel0_stats[nf_bin]);
        double fcorrel0_error = std::get<binning::_SQERROR>(fcorrel0_stats[nf_bin]);

        triqs::arrays::array<double, 2> fcorrel_out(lattice_.dims[0] / 2, 5);
        for (int l = 0; l < lattice_.dims[0] / 2; l++) { 
            auto fcorrel_stats = jackknife::jack(std::function<double(std::vector<double>)>(std::bind(fcorrel_f, std::placeholders::_1, l)),fhistory,nf_bin);
            save_bin_data(fcorrel_stats,h5_stats_,"fcorrel_" + std::to_string(l),save_plaintext);
            double fcorrel_mean = std::get<binning::_MEAN>(fcorrel_stats);
            double fcorrel_error = std::get<binning::_SQERROR>(fcorrel_stats);
            fcorrel_out(l,0) = l;
            fcorrel_out(l,1) = fcorrel_mean;
            fcorrel_out(l,2) = fcorrel_error;
            fcorrel_out(l,3) = fcorrel_mean / fcorrel0_mean;
            fcorrel_out(l,4) = std::sqrt(std::pow(fcorrel_error / fcorrel0_mean, 2) + std::pow(fcorrel_mean / (fcorrel0_mean * fcorrel0_mean) * fcorrel0_error, 2));
            }

        h5_write(h5_stats_,"fcorrel",fcorrel_out);
        if (save_plaintext) savetxt("fcorrel.dat",fcorrel_out);
}


template <typename MC>
void data_saver<MC>::save_conductivity(std::vector<double> wgrid_cond)
{
    bool save_plaintext = p_["save_plaintext"];
    INFO("Stiffness");
    const auto& stiffness = observables_.stiffness;
    size_t size = stiffness.size();
    auto stiffness_data = binning::accumulate_binning(stiffness.rbegin(), stiffness.rend(), max_bin_); 
    save_binning(stiffness_data,h5_binning_,h5_stats_,"stiffness",save_plaintext);

    // conductivity
    std::vector<std::vector<double>> const& cond_history = observables_.cond_history;
    std::vector<double> const& wgrid = wgrid_cond; 

    // cond(w)
    auto cond_stats0 = binning::accumulate_binning(cond_history[wgrid.size()/2].begin(), cond_history[wgrid.size()/2].end(), max_bin_);
    int cond_bin = estimate_bin(cond_stats0); 

    INFO("Saving w*conductivity (w)");
    triqs::arrays::array<double, 2> wcond_ev(wgrid.size(),3);
    for (size_t i=0; i<wgrid.size(); i++) {
        std::complex<double> z = wgrid[i];
        auto cond_binning_data = binning::bin(cond_history[i].begin(), cond_history[i].end(), cond_bin);
        wcond_ev(i,0) = std::real(z); 
        wcond_ev(i,1) = std::get<binning::bin_m::_MEAN>(cond_binning_data);
        wcond_ev(i,2) = std::get<binning::bin_m::_SQERROR>(cond_binning_data); 
        }
    h5_write(h5_stats_,"wcond_err",wcond_ev);
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
    h5_write(h5_stats_,"cond0",t_cond0);
    if (save_plaintext) savetxt("cond0.dat",t_cond0);

    triqs::arrays::array<double, 2> cond_ev(wcond_ev), cond_ev_subtract(wcond_ev);
    for (size_t i=0; i<wgrid.size(); i++) { 
        cond_ev(i,1) = -wcond_ev(i,1) / wgrid[i];
        cond_ev(i,2) = -wcond_ev(i,2) / wgrid[i];
        cond_ev_subtract(i,1) = cond_ev(i,1) - cond0;
        cond_ev_subtract(i,2) = cond_ev(i,2);
        }
    
    h5_write(h5_stats_,"cond_err",cond_ev);
    if (save_plaintext) savetxt("cond_err.dat",cond_ev);
    h5_write(h5_stats_,"cond_dynamic",cond_ev_subtract);
    if (save_plaintext) savetxt("cond_dynamic.dat",cond_ev_subtract);
}

template <typename MC>
void data_saver<MC>::save_energy()
{
    double beta = mc_.config.params().beta;

    if (observables_.energies.size()) { 
        const std::vector<double>& energies = observables_.energies;
        const std::vector<double>& d2energies = observables_.d2energies;
        size_t energy_bin = 0;

        { /* Energy binning */
            INFO("Energy binning");
            size_t size = energies.size();
            INFO("Binning " << max_bin_ <<" times.");
            auto energy_binning_data = binning::accumulate_binning(energies.rbegin(), energies.rend(), max_bin_); 
            save_binning(energy_binning_data,h5_binning_,h5_stats_,"energy",p_["save_plaintext"]);
            auto d2energy_binning_data = binning::accumulate_binning(d2energies.rbegin(),d2energies.rend(), max_bin_); 
            save_binning(d2energy_binning_data,h5_binning_,h5_stats_,"d2energy",p_["save_plaintext"]);
        };
        
        std::vector<double> energies_square(energies.size());
        std::transform(energies.begin(), energies.end(), energies_square.begin(), [](double x){ return x*x; });
        { /* Specific heat */
            typedef std::function<double(double, double, double)> cf_t;

            cf_t cv_function = [beta, this](double e, double e2, double de2){return beta*beta*(e2 - de2 - e*e)/volume_;}; 

            typedef decltype(energies.rbegin()) it_t;
            std::vector<std::pair<it_t,it_t>> c_data = {
                std::make_pair(energies.rbegin(), energies.rend()), 
                std::make_pair(energies_square.rbegin(), energies_square.rend()), 
                std::make_pair(d2energies.rbegin(), d2energies.rend())
            }; 
            auto cv_stats = jackknife::accumulate_jackknife(cv_function,c_data,max_bin_);
            save_binning(cv_stats,h5_binning_,h5_stats_,"cv",p_["save_plaintext"]);
        };
    };
}

template <typename MC>
void data_saver<MC>::save_observables()
{
    //===== save direct measures ===== //
    if (observables_.energies.size()) h5_write(h5_mc_data_,"energies", observables_.energies);
    if (observables_.d2energies.size()) h5_write(h5_mc_data_,"d2energies", observables_.d2energies);
    if (observables_.nf0.size()) h5_write(h5_mc_data_,"nf0", observables_.nf0);
    if (observables_.nfpi.size()) h5_write(h5_mc_data_,"nfpi", observables_.nfpi);
    if (observables_.stiffness.size()) h5_write(h5_mc_data_,"stiffness", observables_.stiffness);

    std::vector<double> const& spectrum = observables_.spectrum;
    if (spectrum.size()) { 
        h5_write(h5_mc_data_,"spectrum", spectrum);
        };

    if (p_["measure_history"]) { 
        triqs::arrays::array<double, 2> t_spectrum_history(observables_.spectrum_history.size(), observables_.spectrum_history[0].size());
        for (int i=0; i<observables_.spectrum_history.size(); i++)
            for (int j=0; j< observables_.spectrum_history[0].size(); j++)
                t_spectrum_history(i,j) =  observables_.spectrum_history[i][j];
        h5_write(h5_mc_data_,"spectrum_history", t_spectrum_history);

        triqs::arrays::array<double, 2> focc_history(observables_.focc_history.size(), observables_.focc_history[0].size());
        for (int i=0; i<observables_.focc_history.size(); i++)
            for (int j=0; j< observables_.focc_history[0].size(); j++)
                focc_history(i,j) =  observables_.focc_history[i][j];
        h5_write(h5_mc_data_,"focc_history", focc_history);
        };

    // Inverse participation ratio
    if (p_["measure_ipr"] && p_["measure_history"]) {
        std::cout << "Inverse participation ratio" << std::endl;
        triqs::arrays::array<double, 2> t_ipr_history(observables_.ipr_history.size(), observables_.ipr_history[0].size());
        for (int i=0; i<observables_.ipr_history.size(); i++)
            for (int j=0; j< observables_.ipr_history[0].size(); j++)
                t_ipr_history(i,j) =  observables_.ipr_history[i][j];
        h5_write(h5_mc_data_,"ipr_history", t_ipr_history);
        };

    // Conductivity
    if (p_["measure_stiffness"]) {
        std::cout << "Conductivity" << std::endl;
        triqs::arrays::array<double, 2> t_cond_history(observables_.cond_history.size(), observables_.cond_history[0].size());
        for (int i=0; i<observables_.cond_history.size(); i++)
            for (int j=0; j< observables_.cond_history[0].size(); j++)
                t_cond_history(i,j) =  observables_.cond_history[i][j];
        h5_write(h5_mc_data_,"cond_history", t_cond_history);
        };

    if (p_["measure_eigenfunctions"]) {
        std::cout << "Eigenfunctions" << std::endl;
        const auto& eig_hist = observables_.eigenfunctions_history;
        triqs::arrays::array<double, 3> t_eig_history(eig_hist.size(), eig_hist[0].rows(), eig_hist[0].cols() );
        for (int i=0; i<eig_hist.size(); i++)
            for (int j=0; j< eig_hist[0].rows(); j++) 
                for (int k=0; k< eig_hist[0].cols(); k++) 
                    t_eig_history(i,j,k) =  eig_hist[i](j,k);
        h5_write(h5_mc_data_,"eig_history", t_eig_history);
        }
}

template <typename MC>
void data_saver<MC>::save_ipr(std::vector<double> grid_real) 
{
    double beta = mc_.config.params().beta;
    const auto &spectrum_history = observables_.spectrum_history;
    auto ipr_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
            std::complex<double> nom = 0.0, denom = 0.0;
            for (size_t i=0; i<volume_; i++) {
                // as norm4 is measured -> take the power of 4 to extract ipr
                double ipr_state = std::pow(ipr_spec[i+volume_],4);
                denom+=1./(z - ipr_spec[i] + I*offset); 
                nom+=1./(z - ipr_spec[i] + I*offset)*ipr_state; 
                };
            return imag(nom)/imag(denom);
            };

        auto ipr_thermal_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
            double out = 0.0;
            double ipr_state, state_weight;
            for (size_t i=0; i<volume_; i++) {
                ipr_state = std::pow(ipr_spec[i+volume_],4);
                state_weight = beta / (1. + std::exp(beta*ipr_spec[i])) / (1. + std::exp(-beta*ipr_spec[i]));
                out += ipr_state * state_weight / volume_;  
                };
            return out; 
            };
        auto dos_thermal_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
            double out = 0.0;
            double energy_state, state_weight;
            for (size_t i=0; i<volume_; i++) {
                energy_state = ipr_spec[i];
                state_weight = beta / (1. + std::exp(beta*ipr_spec[i])) / (1. + std::exp(-beta*ipr_spec[i]));
                out += state_weight / volume_;  
                };
            return out; 
            };

        auto ipr2_f = [&](const std::vector<double> ipr_spec, std::complex<double> z, double offset)->double { 
            std::complex<double> nom = 0.0, denom = 0.0;
            for (size_t i=0; i<volume_; i++) {
                double e = ipr_spec[i];
                // as norm4 is measured -> take the power of 4 to extract ipr
                double ipr = std::pow(ipr_spec[i+volume_],4);
                denom+=1./(z - ipr_spec[i] + I*offset); 
                nom+=1./(z - ipr_spec[i] + I*offset)*ipr;
                };
            return boost::math::pow<2>(imag(nom)/imag(denom));
            };


        const auto& ipr_vals = mc_.observables.ipr_history;

        typedef std::vector<double>::const_iterator iter_t;
        assert(volume_ == ipr_vals.size());
        std::vector<std::pair<iter_t,iter_t>> ipr_and_spectrum(2*volume_); // create a vector of pair of 2*volume_ size
        for (size_t i=0; i<volume_; ++i) { 
            ipr_and_spectrum[i]=std::make_pair(spectrum_history[i].begin(),spectrum_history[i].end());
            ipr_and_spectrum[i+volume_]=std::make_pair(ipr_vals[i].begin(),ipr_vals[i].end());
        }

        typename binning::bin_data_t ipr0_binning(max_bin_);
        for (int i=0; i<max_bin_; i++) 
            { // save ipr at w=0
                auto ipr0_stats = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(ipr_f, std::placeholders::_1, 0.0, p_["dos_offset"]))
                    ,ipr_and_spectrum,i);

                ipr0_binning[i] = ipr0_stats;
            }
        save_binning(ipr0_binning,h5_binning_,h5_stats_,"ipr0",p_["save_plaintext"]);

        // ipr - thermal
        typename binning::bin_data_t ipr_thermal_binning(max_bin_);
        for (int i=0; i<max_bin_; i++) 
            { // save ipr at w=0
                auto ipr_th_stats = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(ipr_thermal_f, std::placeholders::_1, 0.0, p_["dos_offset"]))
                    ,ipr_and_spectrum,i);

                ipr_thermal_binning[i] = ipr_th_stats;
            }
        save_binning(ipr_thermal_binning,h5_binning_,h5_stats_,"ipr_thermal",p_["save_plaintext"]);

        // dos - thermal
        typename binning::bin_data_t dos_thermal_binning(max_bin_);
        for (int i=0; i<max_bin_; i++) 
            { // save ipr at w=0
                auto dos_th_stats = jackknife::jack(
                    std::function<double(std::vector<double>)>( 
                    std::bind(dos_thermal_f, std::placeholders::_1, 0.0, p_["dos_offset"]))
                    ,ipr_and_spectrum,i);

                dos_thermal_binning[i] = dos_th_stats;
            }
        save_binning(dos_thermal_binning,h5_binning_,h5_stats_,"dos_thermal",p_["save_plaintext"]);


        auto ipr0_bin=estimate_bin(ipr0_binning);
        triqs::arrays::array<double, 2> ipr_ev(grid_real.size(),3);
        for (size_t i=0; i<grid_real.size(); i++) {
            std::complex<double> z = grid_real[i];
            auto ipr_data = jackknife::jack(
                std::function<double(std::vector<double>)>( 
                std::bind(ipr_f, std::placeholders::_1, z, p_["dos_offset"]))
                ,ipr_and_spectrum,ipr0_bin);
            ipr_ev(i,0) = std::real(z); 
                ipr_ev(i,1) = std::get<binning::bin_m::_MEAN>(ipr_data);
                ipr_ev(i,2) = std::get<binning::bin_m::_SQERROR>(ipr_data); 

                };

            h5_write(h5_stats_,"ipr_err",ipr_ev);
            if (p_["save_plaintext"]) savetxt("ipr_err.dat",ipr_ev);



}


template <typename MC>
void data_saver<MC>::save_gwr(std::vector<std::complex<double>> wgrid) 
{
    typedef observables_t::dense_m dense_m;
    const std::vector<dense_m>& eigs = observables_.eigenfunctions_history;
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

    double xi = p_["dos_offset"];
    dense_m xi_m = ident_v * xi;

    for (std::complex<double> w : wgrid) { 
        std::string wstring = std::to_string(float(w.real())) + "_" + std::to_string(float(w.imag()));
        gwr_re.setZero();
        gwr_im.setZero();
        for (int m = 0; m < eigs.size(); ++m) { 
            for (int i = 0; i < eigs[0].rows(); ++i) { evals.diagonal()(i) = observables_.spectrum_history[i][m]; }
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
        auto dims = lattice_.dims;
        dense_m gwr_im2 = dense_m::Zero(dims[0], dims[1]);
        dense_m gwr_re2 = dense_m::Zero(dims[0], dims[1]);
        for (int i=0; i<gwr_im.rows(); ++i) { 
            auto pos_i = lattice_.index_to_pos(i);
            for (int j =0; j < gwr_im.cols(); ++j) { 
                auto pos_j = lattice_.index_to_pos(j);
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



} // end namespace fk