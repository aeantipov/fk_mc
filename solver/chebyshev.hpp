#pragma once

#include <functional>
#include <string>

namespace fk {

inline double chebyshev_t(double x, size_t n)
{
    #ifndef NDEBUG
        if (std::abs(x)>1) throw std::logic_error("Argument " + std::to_string(x) + " of Chebyshev polynomial is not within -1 and 1 bounds");
    #endif
    return cos(n*acos(x));
}

template <typename F, typename Grid>
inline double simpson(F&& f, const Grid& v)
{
    int n = v.size();
    typedef decltype(f(v[0])) value_type; 
    value_type s=0;
    assert(n>=8);
    for (int i=4;i<n-4;i++) { s += 48.*f(v[i]); }
    auto dx = v[1] - v[0];
    return dx/48.*(17.*f(v[0]) + 59.*f(v[1]) + 43.*f(v[2])+49.*f(v[3]) + s + 49. *f(v[n-4]) + 43.*f(v[n-3]) + 59.*f(v[n-2]) + 17.*f(v[n-1]));
}

template <typename F, typename Grid>
inline double chebyshev_moment(const F& op, int order, const Grid& grid)
{
    auto ch = std::bind(chebyshev_t, std::placeholders::_1, order);
    auto f = [&](double w){return 1./M_PI/sqrt(1-w*w)*op(w)*ch(w);};
    auto d = simpson(f,grid);
    return d; 
}
/*
template <typename F, typename Grid>
inline double chebyshev_moment2(const F& op, int order, const Grid& grid)
{
    auto ch = std::bind(chebyshev_t, std::placeholders::_1, order);
    auto f = [&](double w){return 1./M_PI*op(w)*ch(w);};
    double s = 0.0;
    for (size_t i =0; i<grid.size(); ++i) s+=f(grid[i]);
    return s; 
}*/



} // end of namespace fk
