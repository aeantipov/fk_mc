#pragma once

#include <string> 

namespace alps {

/// A small structure to output calculation progress
struct percent_output {
    std::string operator()(double fraction, int step) const {
        int new_percent = fraction * 100.;
        if (new_percent >= last_percent_) {
            last_percent_ = new_percent + step;
            return std::to_string(new_percent) + "% ";
        }
        else { return ""; }
    }
    std::string operator()(double fraction) const { return (*this)(fraction, step_); }
    percent_output(int step = 3, int begin_percentage = 0) : last_percent_(begin_percentage), step_(step) { }
    mutable int last_percent_ = 0;
    int step_ = 3;
};

} // end of namespace bold_hyb

