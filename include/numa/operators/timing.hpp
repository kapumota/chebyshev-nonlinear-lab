#pragma once

// Instrumentación de tiempo por componente para las campañas C1.
// Timed mide F, JVP y D2. TimedPrecond mide construcción dinámica y
// aplicaciones del precondicionador. La construcción estática del proveedor
// se mide en el ejecutor porque debe ocurrir una vez por repetición.

#include "operator.hpp"

#include <chrono>
#include <type_traits>
#include <utility>

namespace numa::op {

struct Timers {
    double t_F = 0.0;
    double t_JVP = 0.0;
    double t_D2 = 0.0;
    double t_prec = 0.0;
    double t_setup_dynamic = 0.0;

    void reset() { *this = Timers{}; }
};

namespace timing_detail {

using Clock = std::chrono::steady_clock;

inline double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

} // namespace timing_detail

template <class Op>
class Timed {
public:
    explicit Timed(const Op& op) : op_(op) {}

    std::size_t size() const { return op_.size(); }

    void eval(const DVec& x, DVec& out) const {
        const auto start = timing_detail::Clock::now();
        op_.eval(x, out);
        timers_.t_F += timing_detail::elapsed(start);
        ++counters_.n_F;
    }

    void jvp(const DVec& x, const DVec& v, DVec& out) const {
        const auto start = timing_detail::Clock::now();
        op_.jvp(x, v, out);
        timers_.t_JVP += timing_detail::elapsed(start);
        ++counters_.n_JVP;
    }

    void d2(const DVec& x, const DVec& v, DVec& out) const
        requires SecondOrderOperator<Op>
    {
        const auto start = timing_detail::Clock::now();
        op_.d2(x, v, out);
        timers_.t_D2 += timing_detail::elapsed(start);
        ++counters_.n_D2;
    }

    const OpCounters& counters() const { return counters_; }
    const Timers& timers() const { return timers_; }
    Timers& timers() { return timers_; }

private:
    const Op& op_;
    mutable OpCounters counters_;
    mutable Timers timers_;
};

template <class Provider>
class TimedPrecond {
    using Inner = std::decay_t<
        decltype(std::declval<const Provider&>().make(
            std::declval<const DVec&>()))>;

public:
    TimedPrecond(const Provider& provider, Timers& timers)
        : provider_(provider), timers_(&timers) {}

    struct Bound {
        Inner inner;
        Timers* timers = nullptr;

        void operator()(const DVec& r, DVec& z) const {
            const auto start = timing_detail::Clock::now();
            inner(r, z);
            timers->t_prec += timing_detail::elapsed(start);
        }
    };

    Bound make(const DVec& x) const {
        const auto start = timing_detail::Clock::now();
        Inner inner = provider_.make(x);
        timers_->t_setup_dynamic += timing_detail::elapsed(start);
        return Bound{std::move(inner), timers_};
    }

private:
    const Provider& provider_;
    Timers* timers_;
};

} // namespace numa::op
