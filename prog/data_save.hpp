#pragma once

#include <iostream>
#include <fstream>

#include "fk_mc.hpp"
#include "binning.hpp"
#include "jackknife.hpp"

namespace fk {

inline void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
}

template<typename MC>
class data_saver { 
    typedef MC mc_t;
    typedef typename mc_t::lattice_type lattice_t;
    typedef typename mc_t::config_t config_t; 
protected:
    mc_t const& mc_;
    lattice_t const& lattice_;
    const observables_t& observables_;
    triqs::utility::parameters p_;
    double volume_;
    long nmeasures_;
    int max_bin_;
    /// top of the archive
    triqs::h5::group top_; 
    /// h5 group with actual measurements
    triqs::h5::group h5_mc_data_; 
    /// h5 group with statistics
    triqs::h5::group h5_stats_; 
    /// h5 group with binned data
    triqs::h5::group h5_binning_; 

    /// Save measured observables, no postprocessing
    void save_observables();
    void save_energy();
    /// Save Inverse Participation Ratio
    void save_ipr(std::vector<double> grid_real);
    /// Save G(w,r) and G(w,k) to plaintext files
    // Warning: only works for 2d
    void save_gwr(std::vector<std::complex<double>> wgrid);
    
public:
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
    
    data_saver(const mc_t& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false) : mc_(mc), lattice_(mc.lattice), observables_(mc_.observables), p_(p) { 
        p_.update(save_defaults());
        p_["output_file"] = output_file; 
        p_["save_plaintext"] = save_plaintext; 
        volume_ = mc.lattice.get_msize(); 
        nmeasures_ = observables_.nfpi.size();
        max_bin_ = std::min(15,std::max(int(std::log(double(nmeasures_)/16)/std::log(2.)-1),1));
    }

    void save_all(std::vector<double> wgrid_cond = {0.0});
};

// save data from solver to hdf5 file 
template <typename MC>
void save_all_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false, std::vector<double> wgrid_cond = {0.0})
{
    data_saver<MC> saver(mc, p, output_file, save_plaintext);
    saver.save_all(wgrid_cond);
}

void savetxt (std::string fname, const triqs::arrays::array<double,1>& in);
void savetxt (std::string fname, const triqs::arrays::array<double,2>& in);

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

inline void save_bin_data(const binning::bin_stats_t& data, triqs::h5::group& h5_group, std::string name, int bin, bool save_plaintext = false)
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

inline void save_binning(const binning::bin_data_t& binning_data, triqs::h5::group& h5_group, triqs::h5::group& h5_stats, std::string name, bool save_plaintext = false)
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

} // end of namespace fk

#include "data_save.hxx"

