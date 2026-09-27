#include <catch2/catch_test_macros.hpp>
#include <numa/methods/newton.hpp>
#include <numa/methods/broyden.hpp>
#include <numa/problems/chemical_reactor.hpp>
#include <cmath>

using namespace numa;
using namespace numa::problems;

// Golden test: Tabla 1.1 de la tesis
// Newton 2x2, x0 = (1.5, 2) converge a (1,1) en 6 iteraciones
TEST_CASE("Golden Newton 2x2 desde (1.5,2)", "[golden]") {
    ChemicalReactor f;
    Vec<2> x0 = {1.5, 2.0};
    auto r = newton<2>(f, x0, 1e-15, 10);

    REQUIRE(std::abs(r.x[0] - 1.0) < 1e-10);
    REQUIRE(std::abs(r.x[1] - 1.0) < 1e-10);
    REQUIRE(r.iterations >= 5);
    REQUIRE(r.iterations <= 8);
}

// Golden test: Tabla 1.3 (Broyden)
TEST_CASE("Golden Broyden 2x2 desde (1.5,2)", "[golden]") {
    ChemicalReactor f;
    Vec<2> x0 = {1.5, 2.0};
    auto r = broyden<2>(f, x0, 1e-15, 20);

    REQUIRE(std::abs(r.x[0] - 1.0) < 1e-6);
    REQUIRE(std::abs(r.x[1] - 1.0) < 1e-6);
}

// Verificar que la solución del reactor químico satisface F(x*) = 0
TEST_CASE("Reactor quimico F(x*) ~ 0", "[golden]") {
    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = newton<2>(f, x0, 1e-14, 20);
    auto F = f.eval(r.x);
    REQUIRE(std::abs(F[0]) < 1e-10);
    REQUIRE(std::abs(F[1]) < 1e-10);
}
