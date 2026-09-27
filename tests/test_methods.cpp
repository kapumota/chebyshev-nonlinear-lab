#include <catch2/catch_test_macros.hpp>
#include <numa/methods/newton.hpp>
#include <numa/methods/newton_df.hpp>
#include <numa/methods/broyden.hpp>
#include <numa/methods/chebyshev.hpp>
#include <numa/methods/newton_omp.hpp>
#include <numa/methods/chebyshev_omp.hpp>
#include <numa/problems/chemical_reactor.hpp>
#include <numa/problems/chandrasekhar.hpp>
#include <cmath>

using namespace numa;
using namespace numa::problems;

TEST_CASE("Newton 2x2 converge", "[methods]") {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = newton<2>(f, x0);
    REQUIRE(r.er_abs < 1e-10);
    REQUIRE(std::abs(r.x[0] - 0.831437753041991) < 1e-6);
    REQUIRE(std::abs(r.x[1] - 0.455654808048899) < 1e-6);
}

TEST_CASE("Newton-DF 2x2 converge", "[methods]") {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = newton_df<2>(f, x0);
    REQUIRE(r.er_abs < 1e-8);
}

TEST_CASE("Broyden 2x2 converge", "[methods]") {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = broyden<2>(f, x0);
    REQUIRE(r.er_abs < 1e-10);
}

TEST_CASE("Chebyshev 2x2 converge", "[methods]") {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = chebyshev<2>(f, x0);
    REQUIRE(r.er_abs < 1e-10);
}

TEST_CASE("Chebyshev 8x8 Chandrasekhar converge", "[methods]") {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;
    auto r = chebyshev<8>(f, x0, 1e-10, 20);
    REQUIRE(r.er_abs < 1e-8);
}

TEST_CASE("Newton-OMP y Chebyshev-OMP coinciden con secuencial",
          "[methods][omp]") {
    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;

    auto r_seq = chebyshev<8>(f, x0, 1e-10, 20);
    auto r_omp = chebyshev_omp<8>(f, x0, 1e-10, 20);

    // Coinciden en al menos 6 dígitos
    for (std::size_t i = 0; i < 8; ++i)
        REQUIRE(std::abs(r_seq.x[i] - r_omp.x[i]) < 1e-6);
}
