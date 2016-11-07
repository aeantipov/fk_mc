#include "fk_mc/mc_metropolis.hpp"

namespace alps {

namespace extra {
// sign function
template<typename T>
int sgn(T val) { return (T(0) < val) - (val < T(0)); }
}

mc_metropolis::parameters_type &mc_metropolis::define_parameters(parameters_type &p) {
    alps::mcbase::define_parameters(p);
    p.define<int>("nsweeps", 1024, "Total number of sweeps (1 sweep = #sweep_len moves + 1 measurement)")
        .define<int>("sweep_len", 16, "Number of moves between subsequent measurements")
        .define<int>("ntherm_sweeps", 1, "How many sweeps to do before start measuring")
        .define<bool>("show_output", true, "Show run progress of mc");
    p["nprocs"] = 1;
    return p;
}

mc_metropolis::mc_metropolis(parameters_type const &p, int rank) :
// warning : define_parameters should be called first
    alps::mcbase(p, rank),
#warning SEED is signed (alpscore)
    random(p["SEED"].as<long>() + rank),
    rank_(rank),
    measure_sweeps_(p["nsweeps"]),
    sweep_len_(p["sweep_len"]),
    thermalization_sweeps_(p["ntherm_sweeps"]),
    nprocs_(p["nprocs"]) {
    moves_.reserve(20);
}

void mc_metropolis::update() {
    if (!moves_.size()) {
        std::cout << ALPS_STACKTRACE;
        throw std::logic_error("No registered moves");
    }
    // Select a random move and do it
    for (size_t m = 0; m < sweep_len_; m++) {
        auto move_index = move_distrib_(random);
        mc_weight_t weight = moves_[move_index].attempt();
        if (std::abs(weight) > metropolis_distrib_(random)) {
            weight *= moves_[move_index].accept();
            naccept_++;
            phase_ *= extra::sgn(weight);
            assert(std::abs(std::abs(phase_) - 1.0) < 1e-8);
        }
        else moves_[move_index].reject();
    };
    sweep_count_++;
};

void mc_metropolis::measure() {
    if (measure_count_ >= thermalization_sweeps_) {
        for (auto &measure : measures_) {
            measure.second.accumulate(phase_);
        }
    }
    measure_count_++;
}

void mc_metropolis::collect_results(boost::mpi::communicator const &c)
{
    for (auto &measure : measures_) {
        measure.second.collect_results(c);
    }
}

inline std::string percent(int x, int y, size_t step = 2) {
    static int last_percent = -int(step);
    int new_percent = double(x) / double(y) * 100.;
    if (new_percent != last_percent && new_percent >= last_percent + step) {
        last_percent = new_percent;
        return std::to_string(new_percent) + "% ";
    }
    else { return ""; }
}

double mc_metropolis::fraction_completed() const {
    double f = double(sweep_count_) / double(measure_sweeps_ + thermalization_sweeps_);
    if (!rank_ && parameters["show_output"].as<bool>()) std::cout << percent(f * nprocs_) << std::flush;
    return f;
}

} // end of namespace alps

