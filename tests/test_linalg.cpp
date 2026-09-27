#include <catch2/catch_test_macros.hpp>
#include <numa/linalg/gauss_seidel.hpp>
#include <numa/linalg/tri_sup.hpp>
#include <numa/linalg/jacobi_omp.hpp>

using namespace numa;

TEST_CASE("gauss_seidel 2x2", "[linalg]") {
    Matrix<2,2> D = {{{2.0, 1.0}, {4.0, 3.0}}};
    Vec<2> d = {3.0, 7.0};
    gauss_seidel(D, d);
    REQUIRE(D[1][0] == 0.0);
    REQUIRE(d[1] == 1.0);
}

TEST_CASE("tri_sup 2x2", "[linalg]") {
    Matrix<2,2> D = {{{2.0, 1.0}, {0.0, 1.0}}};
    Vec<2> d = {3.0, 1.0};
    auto x = tri_sup(D, d);
    REQUIRE(x[1] == 1.0);
    REQUIRE(x[0] == 1.0);
}

TEST_CASE("sistema 3x3 con gauss_seidel+tri_sup", "[linalg]") {
    Matrix<3,3> A = {{{3,2,-1},{2,-2,4},{-1,0.5,-1}}};
    Vec<3> b = {1,-2,0};
    auto D = A; auto d = b;
    gauss_seidel(D, d);
    auto x = tri_sup(D, d);
    auto r = matvec(A, x);
    REQUIRE(std::abs(r[0] - 1.0) < 1e-12);
    REQUIRE(std::abs(r[1] + 2.0) < 1e-12);
    REQUIRE(std::abs(r[2] - 0.0) < 1e-12);
}

TEST_CASE("jacobi_omp converge en diagonal dominante", "[linalg][omp]") {
    Matrix<3,3> A = {{{5,1,0},{1,5,1},{0,1,5}}};
    Vec<3> b = {6,7,6};
    auto x = jacobi_omp(A, b, 500, 1e-12);
    auto r = matvec(A, x);
    REQUIRE(std::abs(r[0] - b[0]) < 1e-8);
    REQUIRE(std::abs(r[1] - b[1]) < 1e-8);
    REQUIRE(std::abs(r[2] - b[2]) < 1e-8);
}
