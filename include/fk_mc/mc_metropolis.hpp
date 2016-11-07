#pragma once
#include <type_traits>
#include <string>
#include <map>
#include <random>
#include <memory>

#include <boost/mpi.hpp>
#include "percent_output.hpp"
#include <alps/mc/mcbase.hpp>

namespace alps {

typedef double mc_weight_t; // make it a template for future purposes

/// A wrap structure to register MC moves. It wraps any other class, who has to have attempt(), accept() and reject() methods.
struct move_wrap {

    move_wrap(move_wrap &&r) noexcept {
        ptr_.swap(r.ptr_);
        attempt_.swap(r.attempt_);
        accept_.swap(r.accept_);
        reject_.swap(r.reject_);
    };
    move_wrap &operator=(move_wrap &&r);
    move_wrap(move_wrap const &r) = default;
    template<typename MoveType, typename = typename std::enable_if<!std::is_convertible<MoveType, move_wrap>::value, std::true_type>::type>
    move_wrap(MoveType &&in);

    mc_weight_t attempt() { return attempt_(); };
    mc_weight_t accept() { return accept_(); };
    void reject() { reject_(); };

    std::shared_ptr<void> ptr_;
    std::function<mc_weight_t(void)> attempt_;
    std::function<mc_weight_t(void)> accept_;
    std::function<void(void)> reject_;
};

/// A wrap structure to make measures in MC. Wraps any class with measure(sign) method.
struct measure_wrap {
    /// Perform a measurement with a given phase.
    void accumulate(mc_weight_t phase) { accumulate_(phase); };
    //measure_wrap(){};
    measure_wrap(measure_wrap &&r) noexcept {
        ptr_.swap(r.ptr_);
        accumulate_.swap(r.accumulate_);
    }
    measure_wrap(const measure_wrap &r) : ptr_(r.ptr_), accumulate_(r.accumulate_) { };
    template<typename MeasureType, typename = typename std::enable_if<!std::is_convertible<MeasureType, measure_wrap>::value, move_wrap>::type>
    measure_wrap(MeasureType &&in);
    template<typename MeasureType>
    MeasureType const &cast() const { return *(static_cast<MeasureType *>(ptr_.get())); }
    template<typename MeasureType>
    MeasureType &cast() { return *(static_cast<MeasureType *>(ptr_.get())); }

    void collect_results(boost::mpi::communicator const &c) { collect_results_(c); } 

    std::shared_ptr<void> ptr_;
    std::function<void(mc_weight_t)> accumulate_;
    std::function<void(boost::mpi::communicator const&)> collect_results_;
};

/** mc_metropolis is a TRIQS-inspired MC adapter for ALPSCore. 
 * It is alps::mcbase with some extra features :
 * it does Metropolis MC and allows to register moves and measures. */
struct mc_metropolis : public alps::mcbase {

    typedef alps::mcbase::observable_collection_type observable_collection_type;
    typedef std::mt19937 random_generator;
    /// Register a move.
    template<typename Move_t>
    bool add_move(Move_t &&move, std::string name, double move_prob = 1.0);
    /// Register a measure.
    template<typename Measure_t>
    bool add_measure(Measure_t &&measure, std::string name);

    /// Define used parameters
    static parameters_type &define_parameters(parameters_type &p);
    /// Constructor from alps::params and random seed offset (typically an MPI rank). Rank = 0 outputs the progress.
    mc_metropolis(parameters_type const &p, int rank = 0);

    /// Perform an update - do sweep_len_ moves.
    virtual void update();
    /// Perform all measures. 
    void measure();
    /// Return an estimate for a completed number of sweeps. Output progress
    double fraction_completed() const;
    /// Return acceptance rate. 
    double acceptance_rate() const { return double(naccept_) / double(sweep_count_ * sweep_len_); }
    /// Return observables
    observable_collection_type &observables() { return measurements; };
    /// Return random generator
    random_generator &rng() { return random; }

    void collect_results(boost::mpi::communicator const &c);

    void print_moves(std::ostream &out = std::cout) const {
        for (int i = 0; i < moves_.size(); i++) {
            out << "move : " << move_names_[i] << ", weight = " << move_probs_[i] << std::endl;
        };
    }

    template<typename MeasureType>
    MeasureType const &extract_measurement(std::string name) const {
        auto it = measures_.find(name);
        if (it == measures_.end()) throw std::logic_error("mc_metropolis : can not extract nonexistent measurement" + name);
        return it->second.cast<MeasureType>();
    }

protected:
    /// Vector of moves.
    std::vector<move_wrap> moves_;
    /// Vector of names of moves.
    std::vector<std::string> move_names_;
    /// Vector of move probability
    std::vector<double> move_probs_;
    /// Map between the name of the measure and the measure
    std::map<std::string, measure_wrap> measures_;

    /** Random generator - superseeds alps::random01. Should be removed as soon as alps adopts standard generator. */
    // has to overwrite name of alps random generator, so name 'random' is used.
    random_generator random;

    /// Rank of the process - taken from the constructor
    int rank_;

    /// Total number of sweeps to make
    long measure_sweeps_ = 1000;
    /// Number of steps done during each sweep
    long sweep_len_ = 1;
    /// Number of sweeps to thermalize
    long thermalization_sweeps_ = 0;
    /// Cycle counter
    long sweep_count_ = 0;
    /// Measure counter
    long measure_count_ = 0;
    /// Count accepted moves
    long naccept_ = 0;
    /// Distribution to choose the move, based on it's weight assigned when adding a move
    std::discrete_distribution<> move_distrib_;
    /// Metropolis distribution
    std::uniform_real_distribution<> metropolis_distrib_ = std::uniform_real_distribution<>(0, 1);
    /// Accumulated sign (or phase_) from accepted moves
    mc_weight_t phase_ = 1.0;
    /// Total amount of procs running - used only for having instructive output on the progress of calculation for MPI jobs.
    size_t nprocs_ = 1;
    percent_output percent;
};


template<typename MeasureType, typename>
//measure_wrap::measure_wrap(typename std::remove_reference<MeasureType>::type &&in) {
measure_wrap::measure_wrap(MeasureType &&in) {
    typedef typename std::remove_reference<MeasureType>::type m_type;
    m_type *m = new m_type(std::forward<MeasureType>(in));
    ptr_.reset(m);
    accumulate_ = [m](mc_weight_t p) { m->accumulate(p); };
}

template<typename MoveType, typename>
// = typename std::enable_if<!std::is_same<MoveType,move_wrap>::value, move_wrap>::type>
move_wrap::move_wrap(MoveType &&in) {
    static_assert(std::is_same<MoveType, decltype(*this)>::value == 0, "");
    //static_assert(std::is_move_assignable<MoveType>::value == 1,"");
    typedef typename std::remove_reference<MoveType>::type m_type;

    m_type *m = new m_type(std::forward<m_type>(in));
    ptr_.reset(m);

    attempt_ = [m, this]() { return m->attempt(); };
    accept_ = [m, this]() { return m->accept(); };
    reject_ = [m, this]() { m->reject(); };
}

template<typename Move_t>
bool mc_metropolis::add_move(Move_t &&move, std::string name, double move_prob) {
    moves_.emplace_back(std::forward<Move_t>(move));
    move_names_.emplace_back(name);
    move_probs_.emplace_back(move_prob);
    move_distrib_ = std::discrete_distribution<>(move_probs_.begin(), move_probs_.end());
    return true;
};

template<typename Measure_t>
bool mc_metropolis::add_measure(Measure_t &&measure, std::string name) {
//#ifdef BOLD_HYB_INTEL_BUGFIX
//    measures_.insert(std::make_pair(name, measure_wrap(std::forward<Measure_t>(measure))));
//#else
    measures_.emplace(std::forward<std::string>(name), std::forward<Measure_t>(measure));
//#endif
    return true;
};

}; // end of namespace alps
