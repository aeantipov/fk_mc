#pragma once

#include <functional>
#include <string>

namespace fk {
namespace chebyshev { 

inline double chebyshev_t(double x, size_t n)
{
    assert(std::abs(x)<=1);
    return cos(n*acos(x));
}

struct chebyshev_eval 
{

    int cheb_size() const { return chebt_cache.rows(); }
    int grid_size() const { return chebt_cache.cols(); }

    chebyshev_eval(int max_moment, int grid_size):
        angle_grid(Eigen::VectorXd::LinSpaced(Eigen::Sequential, grid_size, 0., 1.)),
        lobatto_grid(grid_size),
        chebt_cache(max_moment+max_moment%2, grid_size)
    {
        max_moment = max_moment+max_moment%2;
        std::cout << "Initializing Chebyshev evaluator with " << this->cheb_size() << " Chebyshev polynomials on a Lobatto grid with " 
                  <<   this->grid_size() << " points." << std::endl;
        for (int i=0; i<angle_grid.size(); i++) { 
            //angle_grid[i] = (1.*i/n);//(2.*i + 1.) / (2.*n);
            lobatto_grid[i] = -cos(M_PI*angle_grid[i]);
            for (int k=0; k<max_moment; k++) chebt_cache(k,i) = chebyshev_t(lobatto_grid[i], k); 
            }
    }

    template <typename F>
    inline auto moment(const F& op, int order) const -> 
        typename std::remove_reference<typename std::result_of<F(double)>::type>::type { // trapezoidal
            typedef typename std::remove_reference<typename std::result_of<F(double)>::type>::type value_type;
            std::vector<value_type> vals(angle_grid.size());
            for (size_t i =0; i<angle_grid.size(); ++i) vals[i] = op(lobatto_grid[i]);
            return moment(vals, order); 
            }

    template <class F>// decltype (std::declval<F>()[0])>
        auto moment(F &&in, int order) const -> 
            typename std::remove_reference<decltype (std::declval<F>()[0])>::type {
            typedef typename std::remove_reference<decltype (std::declval<F>()[0])>::type value_type;
            value_type s = 0.0;
            for (size_t i =0; i<angle_grid.size()-1; ++i) { 
                s+=(in[i+1]*chebt_cache(order, i+1) + in[i]*chebt_cache(order, i))*(angle_grid[i+1]-angle_grid[i]);
                };
            return s*0.5; 
        }

protected:
    /// Uniform grid between (-1;1) 
    Eigen::VectorXd angle_grid;
    /// Lobatto grid (x[i] = -cos(PI*angle_grid[i]))  
    Eigen::VectorXd lobatto_grid;
    /// Cache of Chebyshev T polynomials at the points defined by grid
    Eigen::MatrixXd chebt_cache;
};

} // end of namespace chebyshev
} // end of namespace fk
