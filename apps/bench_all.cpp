#include <benchmark/benchmark.h>
#include <numa/methods/newton.hpp>
#include <numa/methods/newton_df.hpp>
#include <numa/methods/broyden.hpp>
#include <numa/methods/kantorovich.hpp>
#include <numa/methods/chebyshev.hpp>
#include <numa/methods/newton_omp.hpp>
#include <numa/methods/chebyshev_omp.hpp>
#include <numa/problems/chemical_reactor.hpp>
#include <numa/problems/chandrasekhar.hpp>

using namespace numa;
using namespace numa::problems;

// ---------- 2x2: reactor químico ----------
static void BM_Newton2x2(benchmark::State& state) {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    for (auto _ : state) {
        auto r = newton<2>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Newton2x2);

static void BM_NewtonDF2x2(benchmark::State& state) {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    for (auto _ : state) {
        auto r = newton_df<2>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_NewtonDF2x2);

static void BM_Broyden2x2(benchmark::State& state) {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    for (auto _ : state) {
        auto r = broyden<2>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Broyden2x2);

static void BM_Chebyshev2x2(benchmark::State& state) {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    for (auto _ : state) {
        auto r = chebyshev<2>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Chebyshev2x2);

// ---------- 8x8: Chandrasekhar ----------
static void BM_Newton8(benchmark::State& state) {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    for (auto _ : state) {
        auto r = newton<8>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Newton8);

static void BM_Kantorovich8(benchmark::State& state) {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    for (auto _ : state) {
        auto r = kantorovich<8>(f, x0, 0.53, 0.265);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Kantorovich8);

static void BM_Chebyshev8(benchmark::State& state) {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    for (auto _ : state) {
        auto r = chebyshev<8>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_Chebyshev8);

#ifdef NUMA_HAS_OPENMP
static void BM_NewtonOMP8(benchmark::State& state) {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    for (auto _ : state) {
        auto r = newton_omp<8>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_NewtonOMP8);

static void BM_ChebyshevOMP8(benchmark::State& state) {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    for (auto _ : state) {
        auto r = chebyshev_omp<8>(f, x0);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_ChebyshevOMP8);
#endif

BENCHMARK_MAIN();
