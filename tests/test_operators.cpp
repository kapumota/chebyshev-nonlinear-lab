#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numa/operators/operator.hpp>
#include <numa/operators/baseline_adapter.hpp>
#include <numa/linalg/gmres.hpp>
#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/methods/chebyshev.hpp>
#include <numa/problems/bratu3d.hpp>
#include <numa/problems/thesis_example.hpp>
#include <cmath>

using namespace numa;
using namespace numa::op;

namespace {

// Punto y dirección deterministas, sin generador aleatorio.
DVec make_point(std::size_t n, double scale) {
    DVec x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = scale * std::sin(0.37 * static_cast<double>(i) + 0.5);
    return x;
}

double relative_error(const DVec& a, const DVec& b) {
    DVec d = a;
    axpy(-1.0, b, d);
    return norm2(d) / norm2(b);
}

} // namespace

TEST_CASE("Bratu3D concepto de operador de segundo orden", "[operators]") {
    STATIC_REQUIRE(NonlinearOperator<Bratu3D>);
    STATIC_REQUIRE(SecondOrderOperator<Bratu3D>);
    STATIC_REQUIRE(SecondOrderOperator<Counted<Bratu3D>>);
}

TEST_CASE("Bratu3D JVP analítico coincide con diferencias finitas", "[operators]") {
    Bratu3D f(6, 1.0);
    const DVec x = make_point(f.size(), 0.3);
    const DVec v = make_point(f.size(), 1.0);
    DVec a, b;
    f.jvp(x, v, a);
    jvp_fd(f, x, v, b);
    REQUIRE(relative_error(b, a) < 1e-6);
}

TEST_CASE("Bratu3D segunda derivada analítica coincide con diferencias finitas", "[operators]") {
    Bratu3D f(6, 1.0);
    const DVec x = make_point(f.size(), 0.3);
    const DVec v = make_point(f.size(), 0.5);
    DVec a, b;
    f.d2(x, v, a);
    d2_fd(f, x, v, b);
    // El error de la diferencia central está dominado por la cancelación de
    // los términos lineales del Laplaciano en F, que son mucho mayores que
    // F''[v,v]. Por eso la tolerancia es mayor que la raíz de epsilon. Este
    // efecto es el objeto de la campaña C2.
    REQUIRE(relative_error(b, a) < 1e-4);
}

TEST_CASE("Adaptador del baseline reproduce J v y F''[v,v]", "[operators]") {
    problems::ThesisExample base;
    BaselineAdapter<2> f(base);
    const DVec x = {1.5, 2.0};
    const DVec v = {0.3, -0.7};
    DVec jv, d2v;
    f.jvp(x, v, jv);
    f.d2(x, v, d2v);
    const auto J = base.jacobian({1.5, 2.0});
    REQUIRE_THAT(jv[0], Catch::Matchers::WithinAbs(J[0][0] * 0.3 + J[0][1] * -0.7, 1e-14));
    REQUIRE_THAT(jv[1], Catch::Matchers::WithinAbs(J[1][0] * 0.3 + J[1][1] * -0.7, 1e-14));
    // Primera ecuación, F1 = x^2 + y^2 - 2, F1''[v,v] = 2 vx^2 + 2 vy^2.
    REQUIRE_THAT(d2v[0], Catch::Matchers::WithinAbs(2 * 0.09 + 2 * 0.49, 1e-14));
}

TEST_CASE("GMRES resuelve el sistema del Jacobiano de Bratu3D", "[linalg]") {
    Bratu3D f(8, 1.0);
    const DVec x0(f.size(), 0.0);
    const DVec b = make_point(f.size(), 1.0);
    DVec sol(f.size(), 0.0);
    auto apply = [&](const DVec& v, DVec& o) { f.jvp(x0, v, o); };
    const auto r = gmres(apply, IdentityPrecond{}, b, sol, 1e-10, 40, 400);
    REQUIRE(r.converged);
    DVec Ax;
    f.jvp(x0, sol, Ax);
    REQUIRE(relative_error(Ax, b) < 1e-9);
}

TEST_CASE("Newton-Krylov converge en Bratu3D", "[methods]") {
    Bratu3D f(8, 1.0);
    const DVec x0(f.size(), 0.0);
    const auto r = newton_krylov(f, x0, KrylovParams{}, BratuJacobi{f});
    REQUIRE(r.converged);
    DVec F;
    f.eval(r.x, F);
    REQUIRE(norm2(F) <= 1e-8 * r.res0 * 1.0001);
}

TEST_CASE("Chebyshev inexacto converge en Bratu3D", "[methods]") {
    Bratu3D f(8, 1.0);
    const DVec x0(f.size(), 0.0);
    const auto r = chebyshev_inexact(f, x0, KrylovParams{}, BratuJacobi{f});
    REQUIRE(r.converged);
    REQUIRE(r.n_D2 >= 1);
    DVec F;
    f.eval(r.x, F);
    REQUIRE(norm2(F) <= 1e-8 * r.res0 * 1.0001);
}

TEST_CASE("Contabilidad coherente entre solver y envoltorio Counted", "[operators]") {
    Bratu3D base(8, 1.0);
    Counted<Bratu3D> f(base);
    const DVec x0(base.size(), 0.0);
    const auto r = chebyshev_inexact(f, x0, KrylovParams{}, BratuJacobi{base});
    REQUIRE(r.n_F == f.counters().n_F);
    REQUIRE(r.n_JVP == f.counters().n_JVP);
    REQUIRE(r.n_D2 == f.counters().n_D2);
}

TEST_CASE("Regresión frente al baseline, Chebyshev 2x2 desde (1.5,2)", "[golden][methods]") {
    problems::ThesisExample base;
    BaselineAdapter<2> f(base);
    KrylovParams p;
    p.rtol_nl = 1e-13;
    p.eta_init = 1e-12;
    p.eta_max = 1e-12;
    p.eta2_factor = 1.0;
    p.eta2_max = 1e-12;
    p.gmres_restart = 2;
    const auto r = chebyshev_inexact(f, {1.5, 2.0}, p);
    REQUIRE(r.converged);

    const auto b = chebyshev<2>(base, {1.5, 2.0}, 1e-14, 20);
    REQUIRE_THAT(r.x[0], Catch::Matchers::WithinAbs(b.x[0], 1e-9));
    REQUIRE_THAT(r.x[1], Catch::Matchers::WithinAbs(b.x[1], 1e-9));
    REQUIRE_THAT(r.x[0], Catch::Matchers::WithinAbs(1.0, 1e-9));
}

TEST_CASE("GMRES informa una ruptura numérica ante un operador nulo", "[linalg]") {
    const DVec b = {1.0, -2.0, 3.0};
    DVec x(b.size(), 0.0);
    auto apply = [](const DVec& v, DVec& out) { out.assign(v.size(), 0.0); };
    const auto r = gmres(apply, IdentityPrecond{}, b, x, 1e-10, 5, 20);
    REQUIRE_FALSE(r.converged);
    REQUIRE(r.breakdown);
}

TEST_CASE("GMRES rechaza reinicio igual a cero sin entrar en un bucle", "[linalg]") {
    const DVec b = {1.0, 2.0};
    DVec x(b.size(), 0.0);
    auto apply = [](const DVec& v, DVec& out) { out = v; };
    const auto r = gmres(apply, IdentityPrecond{}, b, x, 1e-10, 0, 20);
    REQUIRE_FALSE(r.converged);
    REQUIRE(r.breakdown);
}

namespace {

// Operador -Laplaciano_h para verificar el multigrid de forma independiente.
void apply_poisson(std::size_t n, const DVec& v, DVec& out) {
    const double h2inv =
        (static_cast<double>(n) + 1.0) *
        (static_cast<double>(n) + 1.0);

    out.assign(n * n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t id = (i * n + j) * n + k;
                double sum = 6.0 * v[id];

                if (i > 0) sum -= v[id - n * n];
                if (i + 1 < n) sum -= v[id + n * n];
                if (j > 0) sum -= v[id - n];
                if (j + 1 < n) sum -= v[id + n];
                if (k > 0) sum -= v[id - 1];
                if (k + 1 < n) sum -= v[id + 1];

                out[id] = h2inv * sum;
            }
        }
    }
}

} // namespace

TEST_CASE(
    "Multigrid construye la jerarquía y reduce el residual de un ciclo V",
    "[linalg]")
{
    PoissonMG mg(15);
    REQUIRE(mg.levels() == 4); // 15, 7, 3, 1

    const DVec b = make_point(15 * 15 * 15, 1.0);
    DVec z;
    DVec Az;

    mg(b, z);
    apply_poisson(15, z, Az);
    axpy(-1.0, b, Az);

    REQUIRE(norm2(Az) < 0.5 * norm2(b));
}

TEST_CASE(
    "GMRES con multigrid mantiene acotadas las iteraciones al refinar",
    "[linalg]")
{
    std::size_t iterations[2] = {0, 0};
    const std::size_t grids[2] = {15, 31};

    for (int grid = 0; grid < 2; ++grid) {
        const std::size_t n = grids[grid];
        PoissonMG mg(n);
        const DVec b = make_point(n * n * n, 1.0);
        DVec x(n * n * n, 0.0);

        auto apply = [&](const DVec& v, DVec& out) {
            apply_poisson(n, v, out);
        };

        const auto result =
            gmres(apply, mg, b, x, 1e-8, 50, 200);

        REQUIRE(result.converged);
        iterations[grid] = result.iters;
    }

    // Con multigrid el conteo debe permanecer aproximadamente acotado.
    REQUIRE(iterations[1] <= iterations[0] + 3);
}

TEST_CASE(
    "La segunda resolución reciclada converge sin productos ocultos",
    "[linalg]")
{
    Bratu3D base(15, 1.0);
    Counted<Bratu3D> counted(base);
    const DVec x(base.size(), 0.1);
    PoissonMG mg(15);

    auto apply = [&](const DVec& v, DVec& out) {
        counted.jvp(x, v, out);
    };

    const DVec b1 = make_point(base.size(), 1.0);
    DVec b2(base.size());
    for (std::size_t i = 0; i < b2.size(); ++i) {
        b2[i] = b1[i] * b1[i];
    }

    DVec z(base.size(), 0.0);
    KrylovSpace space;
    const auto first =
        gmres(apply, mg, b1, z, 1e-8, 50, 100, true, &space);

    REQUIRE(first.converged);
    REQUIRE(space.k > 0);
    REQUIRE(space.k <= first.iters);

    counted.reset();

    DVec u;
    const auto second =
        gmres_recycled(apply, mg, b2, u, 1e-8, 50, 100, space);

    REQUIRE(second.converged);
    REQUIRE(counted.counters().n_JVP == second.matvecs);

    DVec Ju;
    base.jvp(x, u, Ju);
    REQUIRE(relative_error(Ju, b2) < 2e-8);
}

TEST_CASE(
    "Chebyshev inexacto con reciclaje converge y mantiene la contabilidad",
    "[methods]")
{
    Bratu3D base(15, 1.0);
    Counted<Bratu3D> counted(base);
    const DVec x0(base.size(), 0.0);

    KrylovParams params;
    params.recycle = true;

    const auto result =
        chebyshev_inexact(
            counted, x0, params, PoissonMGProvider(15));

    REQUIRE(result.converged);
    REQUIRE(result.n_JVP == counted.counters().n_JVP);
    REQUIRE(result.n_F == counted.counters().n_F);
    REQUIRE(result.n_D2 == counted.counters().n_D2);
    REQUIRE(result.n_prec > 0);
}
