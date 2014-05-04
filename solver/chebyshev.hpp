#pragma once

#include <functional>
#include <string>

namespace fk {

inline double chebyshev_t(double x, size_t n)
{
    assert(std::abs(x)<=1);
    return cos(n*acos(x));
}

struct chebyshev_eval 
{
    chebyshev_eval(int max_moment, int n):
        angle_grid(Eigen::VectorXd::LinSpaced(Eigen::Sequential, n, 0., 1.)),
        lobatto_grid(n),
        chebt_cache(max_moment, n)
    {
        for (int i=0; i<angle_grid.size(); i++) { 
            //angle_grid[i] = (1.*i/n);//(2.*i + 1.) / (2.*n);
            lobatto_grid[i] = -cos(M_PI*angle_grid[i]);
            for (int k=0; k<max_moment; k++) chebt_cache(k,i) = chebyshev_t(lobatto_grid[i], k); 
            }
    }

    template <typename F>
    inline auto moment(const F& op, int order) -> typename std::remove_reference<typename std::result_of<F(double)>::type>::type const { // trapezoidal
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

} // end of namespace fk
