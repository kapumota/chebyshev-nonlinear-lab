#include <catch2/catch_test_macros.hpp>
#include <numa/methods/tf_chebyshev_2018.hpp>
#include <numa/problems/bratu3d.hpp>

#include <cmath>

using namespace numa::op;

namespace {

double relative_error(const DVec& a, const DVec& b) {
    DVec delta = a;
    axpy(-1.0, b, delta);
    return norm2(delta) / (1.0 + norm2(b));
}

} // namespace

TEST_CASE("Curvatura tensor-free 2018 aproxima D2 analítica en Bratu", "[tf2018]") {
    Bratu3D f(4, 1.0);
    DVec x(f.size()), d(f.size());
    for (std::size_t i = 0; i < f.size(); ++i) {
        x[i] = 0.1 + 0.03 * std::sin(0.2 * static_cast<double>(i));
        d[i] = 0.5 * std::cos(0.17 * static_cast<double>(i));
    }
    DVec tensor_free, exact;
    double h = 0.0;
    tensor_free_curvature_2018(f, x, d, tensor_free, &h);
    f.d2(x, d, exact);
    REQUIRE(h > 0.0);
    REQUIRE(h <= 0.5);
    REQUIRE(relative_error(tensor_free, exact) < 5e-6);
}

TEST_CASE("TF-Chebyshev 2018 converge en Bratu pequeño", "[tf2018][methods]") {
    Bratu3D f(4, 1.0);
    BratuJacobi pc{f};
    KrylovParams p;
    p.rtol_nl = 1e-8;
    p.max_nl = 20;
    p.gmres_restart = 30;
    p.gmres_max_iters = 300;
    const DVec x0(f.size(), 0.0);
    const auto result = tf_chebyshev_2018(f, x0, p, pc);
    REQUIRE(result.converged);
    REQUIRE(result.n_D2 == 0);
    REQUIRE(result.n_JVP > 0);
}
