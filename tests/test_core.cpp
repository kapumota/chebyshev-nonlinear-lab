#include <catch2/catch_test_macros.hpp>
#include <numa/core/vector.hpp>
#include <numa/core/matrix.hpp>

using namespace numa;

TEST_CASE("norm_inf vector", "[core]") {
    Vec<3> v = {1.0, -2.5, 3.0};
    REQUIRE(norm_inf(v) == 3.0);
}

TEST_CASE("norm_2 vector", "[core]") {
    Vec<3> v = {3.0, 4.0, 0.0};
    REQUIRE(norm_2(v) == 5.0);
}

TEST_CASE("operadores vector", "[core]") {
    Vec<2> a = {1.0, 2.0}, b = {3.0, 4.0};
    auto c = a + b;
    REQUIRE(c[0] == 4.0);
    REQUIRE(c[1] == 6.0);
    auto d = a - b;
    REQUIRE(d[0] == -2.0);
    REQUIRE(d[1] == -2.0);
    auto e = 2.0 * a;
    REQUIRE(e[0] == 2.0);
    REQUIRE(e[1] == 4.0);
}

TEST_CASE("matvec identidad", "[core]") {
    auto I = identity<3>();
    Vec<3> x = {1.0, 2.0, 3.0};
    auto y = matvec(I, x);
    REQUIRE(y[0] == 1.0);
    REQUIRE(y[1] == 2.0);
    REQUIRE(y[2] == 3.0);
}

TEST_CASE("norm_inf matriz", "[core]") {
    Matrix<2,2> A = {{{1.0, -2.0}, {3.0, 4.0}}};
    REQUIRE(norm_inf(A) == 7.0);
}
