#include <catch2/catch_test_macros.hpp>
#include <numa/recurrence/sigma.hpp>
#include <numa/recurrence/pol.hpp>
#include <cmath>

using namespace numa;

TEST_CASE("sigma en (0,2)", "[recurrence]") {
    for (double a : {0.05, 0.1, 0.2, 0.3, 0.4}) {
        double s = sigma(a);
        REQUIRE(s > 0.0);
        REQUIRE(s < 2.0);
    }
}

TEST_CASE("pol converge a r1", "[recurrence]") {
    double L = 0.530394219033;
    double a = 0.265197109517;
    double t = 0.0;
    for (int i = 0; i < 30; ++i) t = pol(t, L, a).t;
    double r1 = (1.0 - std::sqrt(1.0 - 2.0*L*a)) / L;
    REQUIRE(std::abs(t - r1) < 1e-10);
}

TEST_CASE("pol raiz de p(t)=0", "[recurrence]") {
    double L = 0.530394219033;
    double a = 0.265197109517;
    double t = 0.0;
    for (int i = 0; i < 50; ++i) t = pol(t, L, a).t;
    // p(t) = (L/2)t^2 - t + a
    double pt = 0.5 * L * t * t - t + a;
    REQUIRE(std::abs(pt) < 1e-12);
}
