// Prueba de humo de R1 sobre Bratu3D. NO es evidencia experimental.
// No sigue un preregistro, usa un solo hilo y una sola repetición, y los
// hiperparámetros no están afinados. Sirve para comprobar que el código corre
// a N de 10^4 a 10^6 y para anticipar riesgos, no para comparar métodos.
//
// Uso: smoke_bratu n1 n2 ...
// Variables:
//   NUMA_LAMBDA, por defecto 1.0.
//   NUMA_PRECOND, valores permitidos: mg o jacobi; por defecto mg.

#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/bratu3d.hpp>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

using namespace numa::op;

namespace {

enum class PreconditionerKind {
    Multigrid,
    Jacobi
};

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

bool parse_preconditioner(
    const char* text, PreconditionerKind& kind)
{
    if (text == nullptr || *text == '\0' ||
        std::strcmp(text, "mg") == 0) {
        kind = PreconditionerKind::Multigrid;
        return true;
    }

    if (std::strcmp(text, "jacobi") == 0) {
        kind = PreconditionerKind::Jacobi;
        return true;
    }

    return false;
}

const char* preconditioner_name(PreconditionerKind kind) {
    return kind == PreconditionerKind::Jacobi
        ? "jacobi"
        : "multigrid";
}

void print_usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s n [n ...]\n"
        "Cada n debe ser un entero positivo y define una malla n x n x n.\n"
        "NUMA_PRECOND admite 'mg' o 'jacobi'.\n",
        program);
}

template <class Provider>
void run_campaign(const Bratu3D& problem, const Provider& preconditioner) {
    const DVec x0(problem.size(), 0.0);

    KrylovParams base_params;
    KrylovParams recycled_params;
    recycled_params.recycle = true;

    KrylovResult newton;
    KrylovResult chebyshev;
    KrylovResult recycled;

    const double t_newton = measure([&] {
        newton = newton_krylov(
            problem, x0, base_params, preconditioner);
    });
    const double t_chebyshev = measure([&] {
        chebyshev = chebyshev_inexact(
            problem, x0, base_params, preconditioner);
    });
    const double t_recycled = measure([&] {
        recycled = chebyshev_inexact(
            problem, x0, recycled_params, preconditioner);
    });

    auto print_row = [&](const char* method,
                         const KrylovResult& result,
                         double elapsed) {
        const double relative_residual =
            result.res0 > 0.0
                ? result.res / result.res0
                : result.res;

        std::printf(
            "%8zu %-15s %3s %5zu %5zu %7zu %5zu %7zu %9.3f %8.1e\n",
            problem.size(), method,
            result.converged ? "si" : "no",
            result.iters, result.n_F, result.n_JVP,
            result.n_D2, result.n_prec, elapsed,
            relative_residual);
    };

    print_row("NK", newton, t_newton);
    print_row("Chebyshev", chebyshev, t_chebyshev);
    print_row("Cheb+reciclaje", recycled, t_recycled);
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

    PreconditionerKind preconditioner_kind =
        PreconditionerKind::Multigrid;
    if (!parse_preconditioner(
            std::getenv("NUMA_PRECOND"), preconditioner_kind)) {
        std::fprintf(
            stderr,
            "Error: NUMA_PRECOND debe ser 'mg' o 'jacobi'.\n");
        return 2;
    }

    std::printf(
        "PRUEBA DE HUMO, sin evidencia experimental. "
        "Un hilo, una repetición, sin preregistro ni afinación.\n");
    std::printf(
        "Parámetro lambda de Bratu: %.17g. Precondicionador: %s.\n",
        lambda, preconditioner_name(preconditioner_kind));
    std::printf(
        "%8s %-15s %3s %5s %5s %7s %5s %7s %9s %8s\n",
        "N", "método", "ok", "iter", "n_F", "n_JVP",
        "n_D2", "n_prec", "t_total", "res_rel");

    for (int arg = 1; arg < argc; ++arg) {
        std::size_t n = 0;
        if (!parse_grid_size(argv[arg], n)) {
            std::fprintf(
                stderr,
                "Error: '%s' no es un tamaño de malla válido.\n",
                argv[arg]);
            return 2;
        }

        Bratu3D problem(n, lambda);

        if (preconditioner_kind == PreconditionerKind::Jacobi) {
            run_campaign(problem, BratuJacobi{problem});
        } else {
            run_campaign(problem, PoissonMGProvider(n));
        }
    }

    return 0;
}

