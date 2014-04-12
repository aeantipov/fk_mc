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
inline double simpson(F&& func, const Grid& grid)
{
    int n = grid.size();
    double a = grid[0], b = grid[n-1];
    double deltaX = (b - a)/n;
    double s = 0.0;
    for (int i = 1;i < n;i+=2) {
        double  xi = a + i * deltaX;
        double  y0 = func( xi - deltaX );
        double  y1 = func( xi );
        double  y2 = func( xi + deltaX );
        s += (y0 + 4*y1 + y2); 
        }
    double out = (deltaX / 3)*s;
    return out;
}

template <typename F, typename Grid>
inline double chebyshev_moment(const F& op, int order, const Grid& grid)
{
    auto ch = std::bind(chebyshev_t, std::placeholders::_1, order);
    auto f = [&](double w){return 1./M_PI/sqrt(1-w*w)*op(w)*ch(w);};
    auto d = simpson(f,grid);
    return d; 
}

} // end of namespace fk
