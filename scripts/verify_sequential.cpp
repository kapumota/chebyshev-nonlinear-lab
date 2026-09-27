// Verifica que las correcciones no rompieron nada:
// 1. GaussLegendre<8> reproduce la tabla de la tesis.
// 2. row()/jacobian_row()/hessian_row() coinciden con
//    eval()/jacobian()/hessian() completos.
// 3. chebyshev_mpi con size=1 converge igual que la version secuencial.
//
// Compilar: g++ -std=c++20 -Iinclude -O2 scripts/verify_sequential.cpp -o /tmp/verify_seq
#include <cstdio>
#include <cmath>
#include <numa/functions/gauss_legendre.hpp>
#include <numa/problems/chandrasekhar.hpp>
#include <numa/problems/chemical_reactor.hpp>
#include <numa/methods/newton.hpp>
#include <numa/methods/chebyshev.hpp>
#include <numa/mpi/chebyshev_mpi.hpp>

using namespace numa;

static int g_fails = 0;

void check(bool cond, const char* msg) {
    if (!cond) { std::printf("[FALLA] %s\n", msg); ++g_fails; }
    else       { std::printf("[ok]    %s\n", msg); }
}

int main() {
    // 1. Tabla de Gauss-Legendre para N=8, contra los valores del Anexo.
    const double t_ref[8] = {0.0198550718, 0.1016667613, 0.2372337950, 0.4082826788,
                              0.5917173212, 0.7627662050, 0.8983332387, 0.9801449282};
    const double w_ref[8] = {0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
                              0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363};
    const auto& t = GaussLegendre<8>::t();
    const auto& w = GaussLegendre<8>::w();
    bool gl_ok = true;
    for (int i = 0; i < 8; ++i) {
        if (std::abs(t[i] - t_ref[i]) > 1e-9) gl_ok = false;
        if (std::abs(w[i] - w_ref[i]) > 1e-9) gl_ok = false;
    }
    check(gl_ok, "GaussLegendre<8> reproduce la tabla del Anexo (tol 1e-9)");

    // 2. Compatibilidad row() vs eval() completo, para varios x.
    numa::problems::Chandrasekhar<8> ch;
    Vec<8> x0{}; for (auto& v : x0) v = 1.0;
    Vec<8> x1{}; for (std::size_t i = 0; i < 8; ++i) x1[i] = 0.5 + 0.1 * static_cast<double>(i);

    for (const auto& x : {x0, x1}) {
        auto F = ch.eval(x);
        auto J = ch.jacobian(x);
        auto H = ch.hessian(x);
        bool row_ok = true;
        for (std::size_t i = 0; i < 8; ++i) {
            if (std::abs(ch.row(i, x) - F[i]) > 1e-14) row_ok = false;
            auto Ji = ch.jacobian_row(i, x);
            for (std::size_t j = 0; j < 8; ++j)
                if (std::abs(Ji[j] - J[i][j]) > 1e-14) row_ok = false;
            auto Hi = ch.hessian_row(i, x);
            for (std::size_t jk = 0; jk < 64; ++jk)
                if (std::abs(Hi[jk] - H[i][jk]) > 1e-14) row_ok = false;
        }
        check(row_ok, "row()/jacobian_row()/hessian_row() == version completa");
    }

    // 3. chebyshev_mpi con size=1 debe converger igual que newton/chebyshev secuenciales.
    auto rn = numa::newton<8>(ch, x0);
    auto rm = numa::mpi::chebyshev_mpi<8>(ch, x0, /*rank=*/0, /*size=*/1);
    bool same = true;
    for (std::size_t i = 0; i < 8; ++i)
        if (std::abs(rn.x[i] - rm.x[i]) > 1e-8) same = false;
    check(same, "chebyshev_mpi(size=1) converge al mismo punto que newton()");

    std::printf("\n%d falla(s)\n", g_fails);
    return g_fails == 0 ? 0 : 1;
}
