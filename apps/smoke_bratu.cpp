// Prueba de humo de R1 sobre Bratu3D. NO es evidencia experimental.
// No sigue un preregistro, usa un solo hilo y una sola repetición, y el
// precondicionador de Jacobi es débil. Sirve para comprobar que el código
// corre a N de 10^4 a 10^5 y para anticipar riesgos, no para comparar métodos.

#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/bratu3d.hpp>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace numa::op;

namespace {

template <class Fn>
double measure(Fn&& fn) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now() - t0)
        .count();
}

bool parse_grid_size(const char* text, std::size_t& value) {
    if (text == nullptr || *text == '\0' || *text == '-' || *text == '+') {
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(text, &end, 10);

    if (errno == ERANGE || end == text || *end != '\0' || parsed == 0) {
        return false;
    }

    if (parsed > std::numeric_limits<std::size_t>::max()) {
        return false;
    }

    const std::size_t n = static_cast<std::size_t>(parsed);
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();

    // Bratu3D usa N = n^3. Se valida antes de construir el problema para
    // evitar desbordamiento de size_t con entradas absurdamente grandes.
    if (n > max_size / n) {
        return false;
    }
    const std::size_t n2 = n * n;
    if (n2 > max_size / n) {
        return false;
    }

    value = n;
    return true;
}

bool parse_lambda(const char* text, double& value) {
    if (text == nullptr || *text == '\0') {
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(text, &end);

    if (errno == ERANGE || end == text || *end != '\0' ||
        !std::isfinite(parsed) || parsed <= 0.0) {
        return false;
    }

    value = parsed;
    return true;
}

void print_usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s n [n ...]\n"
        "Cada n debe ser un entero positivo y define una malla n x n x n.\n",
        program);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 2;
    }

    double lambda = 1.0;
    if (const char* text = std::getenv("NUMA_LAMBDA")) {
        if (!parse_lambda(text, lambda)) {
            std::fprintf(
                stderr,
                "Error: NUMA_LAMBDA debe ser un número positivo y finito.\n");
            return 2;
        }
    }

    std::printf(
        "PRUEBA DE HUMO, sin evidencia experimental. "
        "Un hilo, una repetición, Jacobi, sin preregistro.\n");
    std::printf("Parámetro lambda de Bratu: %.17g.\n", lambda);
    std::printf(
        "%8s %-10s %4s %5s %6s %8s %6s %9s %8s\n",
        "N", "método", "ok", "iter", "n_F", "n_JVP", "n_D2",
        "t_total", "res_rel");

    for (int a = 1; a < argc; ++a) {
        std::size_t n = 0;
        if (!parse_grid_size(argv[a], n)) {
            std::fprintf(
                stderr,
                "Error: '%s' no es un tamaño de malla válido.\n",
                argv[a]);
            return 2;
        }

        Bratu3D f(n, lambda);
        const DVec x0(f.size(), 0.0);
        KrylovParams params;
        KrylovResult nk;
        KrylovResult chebyshev;

        const double t_nk = measure([&] {
            nk = newton_krylov(f, x0, params, BratuJacobi{f});
        });
        const double t_chebyshev = measure([&] {
            chebyshev = chebyshev_inexact(f, x0, params, BratuJacobi{f});
        });

        auto print_row = [&](const char* method,
                             const KrylovResult& result,
                             double elapsed) {
            const double relative_residual =
                result.res0 > 0.0 ? result.res / result.res0 : result.res;
            std::printf(
                "%8zu %-10s %4s %5zu %6zu %8zu %6zu %9.3f %8.1e\n",
                f.size(), method, result.converged ? "si" : "no",
                result.iters, result.n_F, result.n_JVP, result.n_D2,
                elapsed, relative_residual);
        };

        print_row("NK", nk, t_nk);
        print_row("Chebyshev", chebyshev, t_chebyshev);
    }

    return 0;
}

