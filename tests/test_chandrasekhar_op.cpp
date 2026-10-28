#include <catch2/catch_test_macros.hpp>

#include <numa/functions/gauss_legendre_dynamic.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/operators/operator.hpp>
#include <numa/problems/chandrasekhar_op.hpp>

#include <cmath>
#include <stdexcept>

using namespace numa::op;

namespace {

DVec deterministic_vector(std::size_t n, double scale) {
    DVec values(n);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] =
            1.0 + scale *
            std::sin(0.37 * static_cast<double>(i) + 0.5);
    }
    return values;
}

double relative_error(const DVec& a, const DVec& b) {
    DVec difference = a;
    axpy(-1.0, b, difference);
    const double denominator = norm2(b);
    return denominator > 0.0
        ? norm2(difference) / denominator
        : norm2(difference);
}

} // namespace

TEST_CASE(
    "Gauss-Legendre dinámico integra polinomios hasta grado 15 con N igual a 8",
    "[quadrature]")
{
    const auto rule = numa::quadrature::gauss_legendre_unit(8);
    REQUIRE(rule.nodes.size() == 8);
    REQUIRE(rule.weights.size() == 8);

    for (std::size_t degree = 0; degree <= 15; ++degree) {
        double numerical = 0.0;
        for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
            numerical +=
                rule.weights[i] *
                std::pow(rule.nodes[i], static_cast<int>(degree));
        }
        const double exact =
            1.0 / static_cast<double>(degree + 1);
        REQUIRE(std::abs(numerical - exact) < 5e-14);
    }
}

TEST_CASE(
    "ChandrasekharOp valida parámetros y dimensiones",
    "[operators]")
{
    REQUIRE_THROWS_AS(
        ChandrasekharOp(0, 0.5), std::invalid_argument);
    REQUIRE_THROWS_AS(
        ChandrasekharOp(10, -0.1), std::invalid_argument);
    REQUIRE_THROWS_AS(
        ChandrasekharOp(10, 1.1), std::invalid_argument);

    ChandrasekharOp problem(10, 0.5);
    DVec out;
    REQUIRE_THROWS_AS(
        problem.eval(DVec(9, 1.0), out),
        std::invalid_argument);
}

TEST_CASE(
    "ChandrasekharOp deriva JVP y D2 de forma consistente",
    "[operators]")
{
    ChandrasekharOp problem(24, 0.8);
    STATIC_REQUIRE(SecondOrderOperator<ChandrasekharOp>);

    const DVec x = deterministic_vector(24, 0.2);
    DVec direction(24);
    for (std::size_t i = 0; i < direction.size(); ++i) {
        direction[i] =
            0.5 * std::cos(0.23 * static_cast<double>(i) + 0.1);
    }

    DVec analytical;
    DVec finite_difference;

    problem.jvp(x, direction, analytical);
    jvp_fd(problem, x, direction, finite_difference);
    REQUIRE(relative_error(finite_difference, analytical) < 2e-6);

    problem.d2(x, direction, analytical);
    d2_fd(problem, x, direction, finite_difference);
    REQUIRE(relative_error(finite_difference, analytical) < 2e-4);
}

TEST_CASE(
    "ChandrasekharOp reproduce el caso N igual a 8 y lambda igual a un cuarto de la tesis",
    "[golden][methods]")
{
    ChandrasekharOp problem(8, 0.5);
    const DVec expected{
        1.021719731462,
        1.073186381734,
        1.125724893657,
        1.169753312169,
        1.203071751305,
        1.226490874633,
        1.241524600594,
        1.249448516693
    };

    DVec residual;
    problem.eval(expected, residual);
    REQUIRE(norm2(residual) < 2e-10);

    const DVec x0(8, 1.0);
    const auto nk =
        newton_krylov(problem, x0, KrylovParams{});
    const auto cheb =
        chebyshev_inexact(problem, x0, KrylovParams{});

    REQUIRE(nk.converged);
    REQUIRE(cheb.converged);
    REQUIRE(relative_error(nk.x, expected) < 2e-7);
    REQUIRE(relative_error(cheb.x, expected) < 2e-7);
}
