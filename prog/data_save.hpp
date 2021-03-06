#pragma once

#include <gftools.hpp>
#include <gftools/hdf5.hpp>

#include <iostream>
#include <fstream>
#include <alps/hdf5.hpp>

#include "fk_mc.hpp"
#include "binning.hpp"
#include "jackknife.hpp"

namespace fk {

using gftools::container;

inline void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
}

/** A static class to handle data saving + data postprocessing and saving. */
template<typename MC>
class data_saver { 
    typedef MC mc_t;
    typedef typename mc_t::lattice_type lattice_t;
    typedef typename mc_t::config_t config_t; 
protected:
    mc_t const& mc_;
    lattice_t const& lattice_;
    const observables_t& observables_;
    parameters_t p_;
    size_t volume_;
    long nmeasures_;
    int max_bin_;
    /// top of the archive
    std::string top_; 
    /// h5 group with actual measurements
    std::string h5_mc_data_; 
    /// h5 group with statistics
    std::string h5_stats_; 
    /// h5 group with binned data
    std::string h5_binning_; 
    alps::hdf5::archive ar_; //(output_file, "w");

    /// Save measured observables, no postprocessing.
    void save_measurements();
    /// Save the mean energy and the specific heat.
    void save_energy();
    /// Save f-electron stats.
    void save_fstats();
    /// Save f-f correlation functions.
    void save_fcorrel();
    /// Save the c-electron conductivity.
    void save_conductivity(std::vector<double> wgrid);
    /// Save Inverse Participation Ratio.
    void save_ipr(std::vector<double> wgrid);
    /// Save the local gf.
    void save_glocal(std::vector<double> wgrid_real);
    /// Save G(w,r) and G(w,k) to plaintext files.
    // Warning: only works for 2d.
    void save_gwr(std::vector<std::complex<double>> wgrid, double imag_offset, bool save_only_dos);

    template <int N>
    static double ipr_moment_f(std::vector<double> const& ipr_spec, std::complex<double> z, double offset, int volume, double mean);
    template <int N>
    static double dos_moment_f(std::vector<double> const& ipr_spec, std::complex<double> z, double offset, int volume, double mean);

    template <typename T>
    void h5_write(std::string top, std::string name, T &&t) { alps::hdf5::save(ar_, top + "/" + name, t); } 
    
public:
    static parameters_t& save_defaults(alps::params &pd) {
      pd.define<int>("dos_npts", int(100), "Number of points for dos")
       .define<double>("dos_width", double(6), "Energy window to save dos")
       .define<bool>("measure_ipr", bool(false), "Measure inverse participation ratio")
       .define<double>("dos_offset", double(0.05), "dos offset from the real axis")
       ;
      return pd;
    }
    
    data_saver(const mc_t& mc, parameters_t &p) : mc_(mc), lattice_(mc.lattice()), observables_(mc_.observables), p_(p), ar_(p["output"].as<std::string>(), "w") 
    {
        //p_.update(save_defaults());
        //p_["output_file"] = output_file; 
        volume_ = mc.lattice().msize();
        nmeasures_ = observables_.nfpi.size();
        max_bin_ = std::min(15,std::max(int(std::log(double(nmeasures_)/16)/std::log(2.)-1),1));
    }

    void save_all(std::vector<double> wgrid_cond = {0.0});
};

// save data from solver to hdf5 file 
template <typename MC>
void save_all_data(const MC& mc, parameters_t p, std::vector<double> wgrid_cond = {0.0})
{
    data_saver<MC> saver(mc, p);
    saver.save_all(wgrid_cond);
}

void savetxt (std::string fname, const gftools::container<double,1>& in);
void savetxt (std::string fname, const gftools::container<double,2>& in);

/// Estimate the bin, at which the error bar is saturated.
inline size_t estimate_bin(const fk::binning::bin_data_t& data)
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

inline void save_bin_data(const binning::bin_stats_t& data, alps::hdf5::archive& ar, std::string h5_group, std::string name, int bin, bool save_plaintext = false)
{
    std::array<double, 4> tmp (
        {{double(std::get<binning::_SIZE>(data)), std::get<binning::_MEAN>(data), std::get<binning::_DISP>(data), std::get<binning::_SQERROR>(data) }});
    gftools::container<double, 1> data_arr(4);
    std::copy(tmp.begin(),tmp.end(),data_arr.data());

    std::cout << name << " = " << std::get<binning::_MEAN>(data) << " +/- " << std::get<binning::_SQERROR>(data) 
              << " (" << std::get<binning::_SIZE>(data) << ") samples [" << bin << "]." << std::endl;
    alps::hdf5::save(ar, h5_group+"/"+name, data_arr);
    if (save_plaintext) savetxt(name+"_error.dat", data_arr);
}

inline void save_binning(const binning::bin_data_t& binning_data, alps::hdf5::archive& ar, std::string h5_group, std::string& h5_stats, std::string name, bool save_plaintext = false)
{
    auto cor_lens = binning::calc_cor_length(binning_data);
    gftools::container<double, 2> data_arr( int(binning_data.size()),5 );
    //std::ofstream out; out.setf(std::ios::scientific); //out << std::setw(9);
    //if (save_plaintext) out.open(name+"_binning.dat",std::ios::out);
    for (size_t i=0; i<binning_data.size(); i++){
        auto e = binning_data[i]; double c = cor_lens[i];
        std::array<double, 5> t ({{double(std::get<binning::_SIZE>(e)), std::get<binning::_MEAN>(e), std::get<binning::_DISP>(e), std::get<binning::_SQERROR>(e), c }});
        std::copy(t.begin(),t.end(),&(data_arr[i][0]));
        //if (save_plaintext) out << i << " " << std::get<binning::_SIZE>(e) << "  " << std::get<binning::_MEAN>(e) 
            //<< "  " << std::get<binning::_DISP>(e) << "  " << std::get<binning::_SQERROR>(e) << "  " << c << "\n";
       };
    if (save_plaintext) savetxt(name+"_binning.dat", data_arr);
    //h5_write(h5_group,name, data_arr);
    alps::hdf5::save(ar, h5_group+"/"+name, data_arr);
            
    int data_bin = estimate_bin(binning_data);
    save_bin_data(binning_data[data_bin],ar,h5_stats,name,data_bin,save_plaintext);
}

} // end of namespace fk

#include "data_save.hxx"

