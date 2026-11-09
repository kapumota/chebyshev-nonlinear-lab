#include <catch2/catch_test_macros.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/brusselator1d.hpp>

#include <cmath>

using namespace numa::op;

namespace {

DVec direction(std::size_t n, double scale, double phase) {
    DVec x(n);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = scale * std::sin(0.31 * static_cast<double>(i) + phase);
    }
    return x;
}

double relative_error(const DVec& a, const DVec& b) {
    DVec delta = a;
    axpy(-1.0, b, delta);
    return norm2(delta) / (1.0 + norm2(b));
}

} // namespace

TEST_CASE("Brusselator1DStep tiene estado inicial público no manufacturado", "[brusselator]") {
    Brusselator1DStep f(31, 0.01);
    REQUIRE(f.size() == 62);
    REQUIRE(f.previous().size() == f.size());
    DVec residual;
    f.eval(f.previous(), residual);
    REQUIRE(std::isfinite(norm2(residual)));
    REQUIRE(norm2(residual) > 0.0);
}

TEST_CASE("Brusselator1DStep JVP y D2 concuerdan con diferencias finitas", "[brusselator]") {
    Brusselator1DStep f(31, 0.01);
    STATIC_REQUIRE(SecondOrderOperator<Brusselator1DStep>);
    DVec x = f.previous();
    const DVec shift = direction(f.size(), 0.05, 0.4);
    axpy(1.0, shift, x);
    const DVec d = direction(f.size(), 0.4, 1.2);
    DVec analytic, finite_difference;

    f.jvp(x, d, analytic);
    jvp_fd(f, x, d, finite_difference);
    REQUIRE(relative_error(finite_difference, analytic) < 2e-6);

    f.d2(x, d, analytic);
    d2_fd(f, x, d, finite_difference);
    REQUIRE(relative_error(finite_difference, analytic) < 5e-4);
}

TEST_CASE("Brusselator1DStep converge desde el estado anterior", "[brusselator][methods]") {
    Brusselator1DStep f(31, 0.01);
    BrusselatorJacobi pc{f};
    KrylovParams p;
    p.rtol_nl = 1e-9;
    p.gmres_restart = 30;
    p.gmres_max_iters = 200;

    const auto nk = newton_krylov(f, f.previous(), p, pc);
    const auto ch = chebyshev_inexact(f, f.previous(), p, pc);
    REQUIRE(nk.converged);
    REQUIRE(ch.converged);
}
