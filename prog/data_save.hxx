#pragma once

#include <fftw3.h>
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
    save_measurements();
    
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
    std::vector<std::complex<double>> grid_real2(dos_npts); for (size_t i=0; i<dos_npts; i++) grid_real2[i] = grid_real[i];

    // G(w,k)
    if (p_["measure_eigenfunctions"]) { 
        // Save gw at zero
        save_gwr({std::complex<double>(0,0)}, p_["dos_offset"], false); 
        // Save dos and typical dos
        save_gwr(grid_real2, p_["dos_offset"], true); 
    } 
    // Save glocal
    this->save_glocal(grid_real);
    // Inverse participation ratio
    if (p_["measure_ipr"]) { this->save_ipr(grid_real); }
    // Conductivity and Drude weight
    if (p_["measure_stiffness"]) { this->save_conductivity(wgrid_cond); };
}

template <typename MC>
void data_saver<MC>::save_measurements()
{
    //===== save direct measures ===== //
    if (observables_.energies.size()) h5_write(h5_mc_data_,"energies", observables_.energies);
    if (observables_.d2energies.size()) h5_write(h5_mc_data_,"d2energies", observables_.d2energies);
    if (observables_.c_energies.size()) h5_write(h5_mc_data_,"c_energies", observables_.c_energies);
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

    if (bool(p_["measure_eigenfunctions"]) && bool(p_["save_eigenfunctions"])) {
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

//
// save processed data (not required for loading)
//

template <typename MC>
void data_saver<MC>::save_energy()
{
    double beta = mc_.config.params().beta;

    if (observables_.energies.size()) { 
        const std::vector<double>& energies = observables_.energies;
        const std::vector<double>& d2energies = observables_.d2energies;
        const std::vector<double>& c_energies = observables_.c_energies;
        size_t energy_bin = 0;

        { /* Energy binning */
            INFO("Energy binning");
            size_t size = energies.size();
            INFO("Binning " << max_bin_ <<" times.");
            auto energy_binning_data = binning::accumulate_binning(energies.rbegin(), energies.rend(), max_bin_); 
            save_binning(energy_binning_data,h5_binning_,h5_stats_,"energy",p_["save_plaintext"]);
            auto d2energy_binning_data = binning::accumulate_binning(d2energies.rbegin(),d2energies.rend(), max_bin_); 
            save_binning(d2energy_binning_data,h5_binning_,h5_stats_,"d2energy",p_["save_plaintext"]);
            if (observables_.c_energies.size()) { 
                auto c_energy_binning_data = binning::accumulate_binning(c_energies.rbegin(),c_energies.rend(), max_bin_); 
                save_binning(c_energy_binning_data,h5_binning_,h5_stats_,"c_energy",p_["save_plaintext"]);
            }
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
template <int N>
double data_saver<MC>::ipr_moment_f(std::vector<double> const& ipr_spec, std::complex<double> z, double offset, int volume, double ipr_mean)
{
    std::complex<double> nom = 0.0, denom = 0.0;
    for (size_t i=0; i<volume; i++) {
        // as Eigen::Matrix::norm(4) was used to measure ipr -> take the power of 4 to extract ipr
        double ipr_state = boost::math::pow<4>(ipr_spec[i+volume]);
        denom+=1./(z - ipr_spec[i] + I*offset); 
        nom+=1./(z - ipr_spec[i] + I*offset)*ipr_state; 
        };
    return boost::math::pow<N>(imag(nom)/imag(denom) - ipr_mean);
}

template <typename MC>
template <int N>
double data_saver<MC>::dos_moment_f(std::vector<double> const& ipr_spec, std::complex<double> z, double offset, int volume, double dos_mean)
{
    std::complex<double> nom = 0.0, denom = 0.0;
    for (size_t i=0; i<volume; i++) {
        nom+=1./(z - ipr_spec[i] + I*offset); 
        };
    return boost::math::pow<N>(-imag(nom)/M_PI - dos_mean);
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

    std::cout << "DOS imag offset = " << p_["dos_offset"] << std::endl; 

    /// Save glocal with error-bars
    if (p_["measure_history"])
    { 
        std::vector<double> grid_imag(std::max(int(beta)*10,1024)); for (size_t i=0; i<grid_imag.size(); i++) grid_imag[i] = PI/beta*(2.*i + 1);
        // dos(w=0)

        std::vector<std::vector<double>> spec_hist_transposed(nmeasures_, std::vector<double>(volume_, 0.0));
        for (size_t m=0; m<nmeasures_; ++m) for (int i=0; i<volume_; ++i) spec_hist_transposed[m][i] = observables_.spectrum_history[i][m];
        const auto &spectrum_history = observables_.spectrum_history;
        //const auto &spectrum_history = spec_hist_transposed; //observables_.spectrum_history;
        size_t norm1 = spectrum_history.size();

        std::vector<double> dos_data(nmeasures_,0.0);
        for (size_t m =0; m<nmeasures_; ++m) { 
            dos_data[m] = dos0_f(spec_hist_transposed[m], 0.0, p_["dos_offset"], volume_);
        }
        auto dos0_stats = binning::accumulate_binning(dos_data.rbegin(), dos_data.rend(), max_bin_);
        save_binning(dos0_stats,h5_binning_,h5_stats_,"dos0",save_plaintext);
        size_t dos_bin = estimate_bin(dos0_stats);
        // dos(w)
        {
            INFO("Saving local DOS w errorbars");
            triqs::arrays::array<double, 2> dos_ev(grid_real.size(),3);
            for (size_t i=0; i<grid_real.size(); i++) {
                std::complex<double> z = grid_real[i];
                for (size_t m=0; m<nmeasures_; ++m) { 
                    dos_data[m] = dos0_f(spec_hist_transposed[m], z, p_["dos_offset"], volume_);
                }
                auto dosz_data = binning::bin(dos_data.rbegin(), dos_data.rend(), dos_bin);
                dos_ev(i,0) = std::real(z); 
                dos_ev(i,1) = std::get<binning::bin_m::_MEAN>(dosz_data);
                dos_ev(i,2) = std::get<binning::bin_m::_SQERROR>(dosz_data); 
                }
            h5_write(h5_stats_,"dos_err",dos_ev);
            if (save_plaintext) savetxt("dos_err.dat",dos_ev);
            }
        //double gw0 = 0.0;
        //for (auto spec : spectrum_history) for (auto e : spec) gw0+=std::imag(1.0/(I*double(p_["dos_offset"]) - e)) / double(volume_) / double(nmeasures_); 
        //std::cout << "Gw0 = " << gw0 << std::endl;
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
void data_saver<MC>::save_ipr(std::vector<double> grid_real) 
{
    double beta = mc_.config.params().beta;

    auto ipr_f = [&](std::vector<double> const& ipr_spec, double z, double offset, int volume)->double  {
        double nom = 0.0, denom = 0.0;
        for (size_t i=0; i<volume; i++) {
            // as Eigen::Matrix::norm(4) was used to measure ipr -> take the power of 4 to extract ipr
            double ipr_state = boost::math::pow<4>(ipr_spec[i+volume]);
            double wminuseps = z - ipr_spec[i]; 
            double lorentz = offset / ( wminuseps * wminuseps + offset * offset);
            denom+=lorentz;
            nom+=lorentz*ipr_state; 
            };
        return nom/denom;
        };

    std::vector<std::vector<double>> ipr_spec_tr(nmeasures_, std::vector<double>(2*volume_, 0.0));
    for (size_t m=0; m<nmeasures_; ++m) for (int i=0; i<volume_; ++i) { 
        ipr_spec_tr[m][i+volume_] = observables_.ipr_history[i][m];
        ipr_spec_tr[m][i] = observables_.spectrum_history[i][m];
        }

    std::vector<double> ipr_data(nmeasures_);
    for (size_t m=0; m<nmeasures_; ++m) for (int i=0; i<volume_; ++i) { 
        ipr_data[m] = ipr_f(ipr_spec_tr[m], 0.0, p_["dos_offset"], this->lattice_.get_msize());
        }

    auto ipr0_binning = binning::accumulate_binning(ipr_data.rbegin(), ipr_data.rend(), max_bin_);
    save_binning(ipr0_binning,h5_binning_,h5_stats_,"ipr0",p_["save_plaintext"]);
    auto ipr0_bin=estimate_bin(ipr0_binning);
    triqs::arrays::array<double, 2> ipr_ev(grid_real.size(),3);
    for (size_t i=0; i<grid_real.size(); i++) {
        double z = grid_real[i];
        for (size_t m=0; m<nmeasures_; ++m) for (int i=0; i<volume_; ++i) { 
            ipr_data[m] = ipr_f(ipr_spec_tr[m], z, p_["dos_offset"], this->lattice_.get_msize());
            }
        auto ipr_bin_data = binning::bin(ipr_data.rbegin(), ipr_data.rend(), ipr0_bin);
        ipr_ev(i,0) = z;
        ipr_ev(i,1) = std::get<binning::bin_m::_MEAN>(ipr_bin_data);
        ipr_ev(i,2) = std::get<binning::bin_m::_SQERROR>(ipr_bin_data); 
        };

    h5_write(h5_stats_,"ipr_err",ipr_ev);
    if (p_["save_plaintext"]) savetxt("ipr_err.dat",ipr_ev);
}

template <typename MC>
void data_saver<MC>::save_gwr(std::vector<std::complex<double>> wgrid, double imag_offset, bool save_only_dos) 
{
    typedef observables_t::dense_m dense_m;
    const std::vector<dense_m>& eigs = observables_.eigenfunctions_history;
    std::cout << "Saving G(w,r)" << std::endl;

    // create a grid for gw
    std::vector<double> grid_real2; 
    grid_real2.push_back(0);

    dense_m gwr_re(volume_, volume_); 
    dense_m gwr_im(volume_, volume_); 
    typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> diag_m;
    diag_m wminuseps(volume_);
    diag_m denom(volume_);
    diag_m evals(volume_);
    Eigen::ArrayXi fconf(volume_); 
    fconf.setZero(); 

    bool save_plaintext = p_["save_plaintext"];

    double xi = imag_offset;
    std::cout << "measurements : " << nmeasures_ << std::endl;
    std::cout << "Gwr imag offset = " << xi << std::endl;

    triqs::arrays::array<double, 2> tdos_vals(wgrid.size(), 8);

    for (int windex = 0; windex < wgrid.size(); ++windex) { 
        std::complex<double> w = wgrid[windex];  
        std::cout << "w = " << w << std::endl;
        w+=I*xi;
        std::string wstr_re = std::to_string(float(w.real()));
        std::string wstr_im = std::to_string(float(w.imag()));
        std::string wstring = 
            wstr_re.erase ( wstr_re.find_last_not_of("0") + 1, std::string::npos ) + "_" + 
            wstr_im.erase ( wstr_im.find_last_not_of("0") + 1, std::string::npos ); 
            

        gwr_re.setZero();
        gwr_im.setZero();
        double gw_test = 0.0;
        Eigen::ArrayXd evals1(volume_);
        for (int m = 0; m < nmeasures_; ++m) { 
            for (int i = 0; i < volume_; ++i) { 
                double eps = observables_.spectrum_history[i][m];
                double wmeps = w.real() - eps;
                evals.diagonal()(i) = eps;
                denom.diagonal()(i) = 1.0 / ( (w.real() - eps)*(w.real() - eps) + w.imag()*w.imag()); 
                wminuseps.diagonal()(i) = w.real() - eps;
                fconf[i] = observables_.focc_history[i][m];
                gw_test += std::imag(1.0 / (w - observables_.spectrum_history[i][m])) / double(volume_) / double(nmeasures_);
            }

            #ifndef NDEBUG
            // below is a check that (w - H)^{-1} is summed correctly
            configuration_t config1(lattice_,p_["beta"],p_["U"],p_["mu_c"],p_["mu_f"]);
            config1.f_config_=fconf;
            config1.calc_hamiltonian();
            config1.calc_ed(true);
            assert(sqrt(Eigen::VectorXd(config1.ed_data().cached_spectrum - Eigen::ArrayXd(evals.diagonal())).squaredNorm()) < 1e-5);
            // H = Evecs * Evals * Evecs^T
            assert(sqrt(Eigen::MatrixXd(config1.hamilt_) - eigs[m] * evals * eigs[m].transpose().squaredNorm()) < 1e-5);

            std::cout << "H - U.E.U^+ diff = " << sqrt((Eigen::MatrixXd(config1.hamilt_) - eigs[m] * (Eigen::MatrixXd::Identity(volume_, volume_)*w.real() - dense_m(wminuseps))  * eigs[m].transpose()).squaredNorm()) << std::endl;
            std::cout << "U^+.H.U - E diff = " << sqrt(( eigs[m].transpose() * Eigen::MatrixXd(config1.hamilt_) * eigs[m] - (Eigen::MatrixXd::Identity(volume_, volume_)*w.real() - dense_m(wminuseps))).squaredNorm()) << std::endl;

            std::cout << "[w-H] - U.(w - E).U^+ diff = " << sqrt(( 
                (Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(config1.hamilt_.cast<std::complex<double>>())) - 
                eigs[m].cast<std::complex<double>>() * ( Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(Eigen::MatrixXd(evals).cast<std::complex<double>>()) ) * eigs[m].transpose().cast<std::complex<double>>()
                ).squaredNorm()) << std::endl;

            std::cout << "[w-H]^{-1} - [U.(w - E).U^+]^{-1} diff = " << sqrt(( 
                (Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(config1.hamilt_.cast<std::complex<double>>())).inverse() - 
                (eigs[m].cast<std::complex<double>>() * ( Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(Eigen::MatrixXd(evals).cast<std::complex<double>>()) ) * eigs[m].transpose().cast<std::complex<double>>()).inverse() ).squaredNorm()) << std::endl;

            std::cout << "[w-H]^{-1} - U.(w - E)^{-1}.U^+ diff = " << sqrt(( 
                (Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(config1.hamilt_.cast<std::complex<double>>())).inverse() - 
                eigs[m].cast<std::complex<double>>() * ( Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(Eigen::MatrixXd(evals).cast<std::complex<double>>()) ).inverse() * eigs[m].transpose().cast<std::complex<double>>() ).squaredNorm()) << std::endl;

            Eigen::MatrixXd gw_add = (Eigen::MatrixXcd::Identity(volume_, volume_) * (w) -
                                      Eigen::MatrixXcd(config1.hamilt_.cast<std::complex<double>>())).inverse().imag().template cast<double>();
            std::cout << "Im[[w-H]^{-1}] - Im[U.(w - E)^{-1}.U^+] diff = " << sqrt(( 
                gw_add - 
                (eigs[m].cast<std::complex<double>>() * ( Eigen::MatrixXcd::Identity(volume_, volume_) * (w) - Eigen::MatrixXcd(Eigen::MatrixXd(evals).cast<std::complex<double>>()) ).inverse() * eigs[m].transpose().cast<std::complex<double>>()).imag().template cast<double>() ).squaredNorm()) << std::endl;
            std::cout << "Im[(w-H)^{-1}] - U.(-Im[w]).(w - E)^{-1}.U^+ diff = " << sqrt((gw_add - eigs[m] * (-w.imag()) * denom * eigs[m].transpose()).squaredNorm()) << std::endl;
            #endif
            gwr_im += eigs[m] * (-w.imag()) * denom * eigs[m].transpose() / nmeasures_;
            gwr_re += eigs[m] * wminuseps * denom * eigs[m].transpose() / nmeasures_;
            }

        // get dos'es
        Eigen::VectorXd gw_diag = gwr_im.diagonal() / (-PI);
        double dos_val = gw_diag.sum() / volume_;
        std::cout << "dos0 = " << dos_val << std::endl;
        double dos_err = 0.;
        for (int i = 0; i<volume_; ++i) { gw_diag[i] = std::log(gw_diag[i]); }
        double dos_geom =std::exp(1.0 / volume_ *  gw_diag.sum());
        std::cout << "dos_geom0 = " << dos_geom << std::endl;
        double dos_geom_err = 0;

        tdos_vals(windex, 0) = std::real(w);
        tdos_vals(windex, 1) = std::imag(w);
        tdos_vals(windex, 2) = dos_geom;
        tdos_vals(windex, 3) = dos_geom_err;
        tdos_vals(windex, 4) = dos_val;
        tdos_vals(windex, 5) = dos_err;
        tdos_vals(windex, 6) = dos_geom / dos_val;
        tdos_vals(windex, 7) = std::sqrt(boost::math::pow<2>(dos_geom_err / dos_val) + boost::math::pow<2>(dos_err * dos_geom / dos_val / dos_val));

        // if we don't need gf -> jump to the next w
        if (save_only_dos) continue;

        // if we need all the rest -> not too many points, save them
        triqs::arrays::array<double, 1> tdos_out(8);
        for (int jj=0; jj<8; ++jj) tdos_out(jj) = tdos_vals(windex, jj);
        h5_write(h5_stats_,"tdos"+wstring ,tdos_out);
        

        auto save_hdf5 = [&](Eigen::MatrixXd const& m, std::string name) {
            std::cout << "Saving " << name << " ...";
            triqs::arrays::array<double, 2> out(m.rows(),m.cols()); 
            for (int i = 0; i < m.rows(); ++i) 
                for (int j = 0; j < m.cols(); ++j) 
                    out(i,j) = m(i,j);
            h5_write(h5_stats_,name,out);
            if (save_plaintext) { 
                std::ofstream out_str(name + ".dat");
                out_str << m << std::endl;
                out_str.close();
                }
            std::cout << "done." << std::endl;
        };

        save_hdf5(gwr_re, "gr_full_w"+wstring+"_re");
        save_hdf5(gwr_im, "gr_full_w"+wstring+"_im");
        //std::cout << "test: = " << gw_test << " == " << gwr_im.diagonal().sum()/double(volume_) <<  std::endl;
        //std::cout << "test2: = " << gw_test << " == " << gwr_im(0,0) <<  std::endl;
    
        
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

        save_hdf5(gwr_re2, "gr_w"+wstring+"_re");
        save_hdf5(gwr_im2, "gr_w"+wstring+"_im");

        Eigen::MatrixXcd gwr = gwr_re2.cast<std::complex<double>>() + I*gwr_im2.cast<std::complex<double>>();
            
        fftw_plan p;
        p = fftw_plan_dft_2d(gwr.rows(), gwr.cols(),  
                         reinterpret_cast<fftw_complex*>(gwr.data()),
                         reinterpret_cast<fftw_complex*>(gwr.data()),
                         FFTW_FORWARD, FFTW_ESTIMATE); 
        fftw_execute(p);

        save_hdf5(gwr.real().cast<double>(), "gk_w"+wstring+"_re");
        save_hdf5(gwr.imag().cast<double>(), "gk_w"+wstring+"_im");
        }
    h5_write(h5_stats_,"tdos_gwr",tdos_vals);
}

/* old incorrect stuff
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

    std::cout << "DOS imag offset = " << p_["dos_offset"] << std::endl; 
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
            dos_v(i,0) = std::real(z); 
            gf_re_v(i,0) = std::real(z); 
            gf_re_v(i,1) = std::bind(gf_re_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum); 
            gf_re_v(i,2) = std::bind(gf_im_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum); 
            dos_v(i,1) = std::bind(dos0_f, std::placeholders::_1, z, p_["dos_offset"], lattice_.get_msize())(spectrum);
            };
        h5_write(h5_stats_,"dos",dos_v);
        h5_write(h5_stats_,"gw_re",gf_re_v);
        if (save_plaintext) { savetxt("dos.dat",dos_v); savetxt("gw_refreq.dat", gf_re_v); };
    };

    /// Save glocal with error-bars
    if (p_["measure_history"])
    { 
        std::vector<double> grid_imag(std::max(int(beta)*10,1024)); for (size_t i=0; i<grid_imag.size(); i++) grid_imag[i] = PI/beta*(2.*i + 1);
        // dos(w=0)

        std::vector<std::vector<double>> spec_hist_transposed(nmeasures_, std::vector<double>(volume_, 0.0));
        for (size_t m=0; m<nmeasures_; ++m) for (int i=0; i<volume_; ++i) spec_hist_transposed[m][i] = observables_.spectrum_history[i][m];
        const auto &spectrum_history = observables_.spectrum_history;
        //const auto &spectrum_history = spec_hist_transposed; //observables_.spectrum_history;
        size_t norm1 = spectrum_history.size();

        std::vector<double> dos0_data(nmeasures_,0.0);
        for (size_t m =0; m<nmeasures_; ++m) { 
            dos0_data[m] = dos0_f(spec_hist_transposed[m], 0.0, p_["dos_offset"], volume_);
        }

        auto dos0_stats = jackknife::accumulate_jackknife(

            std::function<double(std::vector<double>)> 
            (std::bind(dos0_f, std::placeholders::_1, 0.0, p_["dos_offset"], volume_))
            ,spectrum_history,max_bin_);
        save_binning(dos0_stats,h5_binning_,h5_stats_,"dos0_old",save_plaintext);


        dos0_stats = binning::accumulate_binning(dos0_data.rbegin(), dos0_data.rend(), max_bin_);
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
                    std::bind(dos0_f, std::placeholders::_1, z, p_["dos_offset"], norm1))
                    ,spectrum_history,dos_bin);
                dos_ev(i,0) = std::real(z); 
                dos_ev(i,1) = std::get<binning::bin_m::_MEAN>(dosz_data);
                dos_ev(i,2) = std::get<binning::bin_m::_SQERROR>(dosz_data); 
                }
            h5_write(h5_stats_,"dos_err",dos_ev);
            if (save_plaintext) savetxt("dos_err.dat",dos_ev);
            }

        double gw0 = 0.0;
        for (auto spec : spectrum_history) for (auto e : spec) gw0+=std::imag(1.0/(I*double(p_["dos_offset"]) - e)) / double(volume_) / double(nmeasures_); 
        std::cout << "Gw0 = " << gw0 << std::endl;
      }
}
*/


} // end namespace fk

