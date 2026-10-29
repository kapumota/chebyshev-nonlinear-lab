#include <catch2/catch_test_macros.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/burgers2d.hpp>
#include <cmath>

using namespace numa::op;

namespace {

DVec point(std::size_t n, double scale, double phase = 0.5) {
    DVec x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = scale * std::sin(0.37 * static_cast<double>(i) + phase);
    }
    return x;
}

double relative_error(const DVec& a, const DVec& b) {
    DVec d = a;
    axpy(-1.0, b, d);
    const double nb = norm2(b);
    return nb > 0.0 ? norm2(d) / nb : norm2(d);
}

} // namespace

TEST_CASE("Burgers2D reproduce la raíz manufacturada", "[burgers]") {
    Burgers2D f(15, 0.05, 0.1, 1.0);
    DVec residual;
    f.eval(f.exact(), residual);
    REQUIRE(norm2(residual) < 1e-10 * (1.0 + norm2(f.exact())));
}

TEST_CASE("Burgers2D JVP y D2 concuerdan con diferencias finitas", "[burgers]") {
    Burgers2D f(15, 0.05, 0.1, 1.0);
    STATIC_REQUIRE(SecondOrderOperator<Burgers2D>);
    const DVec x = point(f.size(), 0.5);
    const DVec v = point(f.size(), 0.7, 1.1);
    DVec analytic, finite_difference;

    f.jvp(x, v, analytic);
    jvp_fd(f, x, v, finite_difference);
    REQUIRE(relative_error(finite_difference, analytic) < 1e-6);

    f.d2(x, v, analytic);
    d2_fd(f, x, v, finite_difference);
    REQUIRE(relative_error(finite_difference, analytic) < 2e-4);
}

TEST_CASE("HelmholtzMG2D construye jerarquía completa en n igual a 31", "[burgers]") {
    HelmholtzMG2D mg(31, 10.0, 0.05);
    REQUIRE(mg.levels() == 5);
}

TEST_CASE("Burgers2D converge con ambos métodos desde perturbación controlada", "[burgers][methods]") {
    Burgers2D f(31, 0.05, 0.1, 1.0);
    BurgersMGProvider pc(31, 0.05, 0.1);
    DVec x0 = f.exact();
    for (std::size_t i = 0; i < x0.size(); ++i) {
        x0[i] += 0.05 * std::sin(1.3 * static_cast<double>(i));
    }
    KrylovParams p;
    const auto nk = newton_krylov(f, x0, p, pc);
    const auto ch = chebyshev_inexact(f, x0, p, pc);
    REQUIRE(nk.converged);
    REQUIRE(ch.converged);
}
