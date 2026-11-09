// Prueba de humo de R1 sobre B5, ecuación H de Chandrasekhar.
// No constituye evidencia experimental: usa una repetición y no aplica
// afinación confirmatoria. El número de hilos depende del entorno OpenMP.
//
// Uso:
//   NUMA_C=0.5 ./smoke_chandrasekhar N1 [N2 ...]
//
// El valor por defecto c = 0.5 reproduce el caso lambda = 1/4 de la tesis.

#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/chandrasekhar_op.hpp>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

using namespace numa::op;

namespace {

bool parse_size(const char* text, std::size_t& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const unsigned long long raw = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || raw == 0) {
        return false;
    }
    if (raw >
        static_cast<unsigned long long>(
            std::numeric_limits<std::size_t>::max())) {
        return false;
    }
    value = static_cast<std::size_t>(raw);
    return true;
}

bool parse_c(const char* text, double& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const double raw = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' ||
        !std::isfinite(raw) || raw < 0.0 || raw > 1.0) {
        return false;
    }
    value = raw;
    return true;
}

template <class Solve>
KrylovResult measure(Solve&& solve, double& elapsed) {
    const auto t0 = std::chrono::steady_clock::now();
    KrylovResult result = solve();
    elapsed =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();
    return result;
}

void print_row(
    std::size_t n,
    const char* method,
    const KrylovResult& result,
    double elapsed)
{
    const double relative_residual =
        result.res0 > 0.0 ? result.res / result.res0 : result.res;

    std::printf(
        "%8zu %-13s %3s %5zu %5zu %7zu %5zu %7zu %7zu %7zu %9.3f %9.1e\n",
        n, method, result.converged ? "si" : "no",
        result.iters, result.n_F, result.n_JVP, result.n_D2,
        result.n_prec, result.n_lin_iters, result.n_lin_failures,
        elapsed, relative_residual);
}

} // namespace

int main(int argc, char** argv) {
    double c = 0.5;
    if (const char* env = std::getenv("NUMA_C"); env != nullptr) {
        if (!parse_c(env, c)) {
            std::fprintf(
                stderr,
                "Error: NUMA_C debe ser un número finito en [0,1].\n");
            return 2;
        }
    }

    if (argc < 2) {
        std::fprintf(stderr, "Uso: %s N1 [N2 ...]\n", argv[0]);
        return 2;
    }

    std::printf(
        "PRUEBA DE HUMO, sin evidencia experimental. "
        "Una repetición, sin afinación confirmatoria.\n");
    std::printf(
        "Chandrasekhar H: c = %.17g. "
        "El caso c = 0.5 corresponde a lambda = 1/4 de la tesis.\n",
        c);
    std::printf(
        "%8s %-13s %3s %5s %5s %7s %5s %7s %7s %7s %9s %9s\n",
        "N", "método", "ok", "iter", "n_F", "n_JVP",
        "n_D2", "n_prec", "n_lin", "fallos",
        "t_total", "res_rel");

    for (int arg = 1; arg < argc; ++arg) {
        std::size_t n = 0;
        if (!parse_size(argv[arg], n)) {
            std::fprintf(
                stderr, "Error: tamaño inválido '%s'.\n", argv[arg]);
            return 2;
        }

        ChandrasekharOp problem(n, c);
        const DVec x0(n, 1.0);

        KrylovParams base;
        KrylovParams recycled = base;
        recycled.recycle = true;

        double elapsed = 0.0;
        auto result = measure(
            [&] { return newton_krylov(problem, x0, base); },
            elapsed);
        print_row(n, "NK", result, elapsed);

        result = measure(
            [&] { return chebyshev_inexact(problem, x0, base); },
            elapsed);
        print_row(n, "Cheb", result, elapsed);

        result = measure(
            [&] { return chebyshev_inexact(problem, x0, recycled); },
            elapsed);
        print_row(n, "Cheb+rec", result, elapsed);
    }

    return 0;
}
