// Continuación espectral para Bratu3D en C1.
//
// Esta herramienta NO define el fold por el fallo de Newton-Krylov. Newton-Krylov
// se usa únicamente como corrector de continuación sobre la rama inferior.
// Cada punto aceptado debe tener una solución convergida y un estimador espectral
// convergido. Si el corrector falla, el paso se reduce; si no puede recuperarse,
// la herramienta termina con error en lugar de declarar un fold.
//
// Uso:
//   calibrate_b6 --n 31 --lambda-start 0.0 --lambda-stop 9.8 --step 0.25 \
//                --min-step 1e-4 --spectral-stop 5e-2 --out bratu_fold.csv

#include <numa/analysis/spectral.hpp>
#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/bratu3d.hpp>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace numa::op;

namespace {

struct Args {
    std::size_t n = 31;
    double lambda_start = 0.0;
    double lambda_stop = 9.8;
    double step = 0.25;
    double min_step = 1e-4;
    double spectral_stop = 5e-2;
    double eigen_rtol = 1e-7;
    double linear_rtol = 1e-10;
    std::string out = "bratu_fold.csv";
};

bool parse_size(const char* text, std::size_t& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const unsigned long long raw = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || raw == 0) return false;
    if (raw > static_cast<unsigned long long>(
                  std::numeric_limits<std::size_t>::max())) {
        return false;
    }
    value = static_cast<std::size_t>(raw);
    return true;
}

bool parse_double(const char* text, double& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const double raw = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(raw)) {
        return false;
    }
    value = raw;
    return true;
}

void usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s [--n N] [--lambda-start L0] [--lambda-stop L1] "
        "[--step H] [--min-step HMIN] [--spectral-stop S] "
        "[--eigen-rtol R] [--linear-rtol R] [--out archivo.csv]\n",
        program);
}

bool parse_args(int argc, char** argv, Args& args) {
    if ((argc - 1) % 2 != 0) return false;
    for (int i = 1; i < argc; i += 2) {
        const std::string key = argv[i];
        const char* value = argv[i + 1];
        if (key == "--n") {
            if (!parse_size(value, args.n)) return false;
        } else if (key == "--lambda-start") {
            if (!parse_double(value, args.lambda_start)) return false;
        } else if (key == "--lambda-stop") {
            if (!parse_double(value, args.lambda_stop)) return false;
        } else if (key == "--step") {
            if (!parse_double(value, args.step)) return false;
        } else if (key == "--min-step") {
            if (!parse_double(value, args.min_step)) return false;
        } else if (key == "--spectral-stop") {
            if (!parse_double(value, args.spectral_stop)) return false;
        } else if (key == "--eigen-rtol") {
            if (!parse_double(value, args.eigen_rtol)) return false;
        } else if (key == "--linear-rtol") {
            if (!parse_double(value, args.linear_rtol)) return false;
        } else if (key == "--out") {
            args.out = value;
        } else {
            return false;
        }
    }

    return args.lambda_start >= 0.0 &&
           args.lambda_stop > args.lambda_start &&
           args.step > 0.0 &&
           args.min_step > 0.0 &&
           args.min_step <= args.step &&
           args.spectral_stop > 0.0 &&
           args.eigen_rtol > 0.0 &&
           args.linear_rtol > 0.0 &&
           !args.out.empty();
}

struct Point {
    double lambda = 0.0;
    double spectral_gap = 0.0;
    double eigen_residual = 0.0;
    std::size_t continuation_iters = 0;
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_prec = 0;
};

bool append_point(std::ofstream& out, const Point& p) {
    out << p.lambda << ','
        << p.spectral_gap << ','
        << p.eigen_residual << ','
        << p.continuation_iters << ','
        << p.n_F << ','
        << p.n_JVP << ','
        << p.n_prec << '\n';
    return static_cast<bool>(out);
}

bool estimate_fold_proxy(
    const std::vector<Point>& points,
    double& lambda_fold_proxy)
{
    if (points.size() < 4) return false;
    const std::size_t first = points.size() - 4;

    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (std::size_t i = first; i < points.size(); ++i) {
        const double x = points[i].lambda;
        const double y =
            points[i].spectral_gap * points[i].spectral_gap;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }

    constexpr double count = 4.0;
    const double denominator = count * sxx - sx * sx;
    if (std::abs(denominator) < 1e-14) return false;

    const double slope = (count * sxy - sx * sy) / denominator;
    const double intercept = (sy - slope * sx) / count;
    if (!(slope < 0.0)) return false;

    const double root = -intercept / slope;
    if (!std::isfinite(root) || root <= points.back().lambda) return false;

    lambda_fold_proxy = root;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) {
        usage(argv[0]);
        return 2;
    }

    PoissonMGProvider provider(args.n);

    KrylovParams continuation_params;
    continuation_params.rtol_nl = 1e-11;
    continuation_params.max_nl = 100;
    continuation_params.eta_init = 1e-8;
    continuation_params.eta_max = 1e-6;
    continuation_params.gmres_restart = 100;
    continuation_params.gmres_max_iters = 2000;

    DVec previous(args.n * args.n * args.n, 0.0);
    double current_lambda = args.lambda_start;
    double current_step = args.step;

    std::ofstream csv(args.out);
    if (!csv) {
        std::fprintf(stderr, "No se pudo abrir %s.\n", args.out.c_str());
        return 2;
    }
    csv << "lambda,spectral_gap,eigen_residual,continuation_iters,"
           "n_F,n_JVP,n_prec\n";

    std::vector<Point> accepted;

    while (current_lambda <= args.lambda_stop + 1e-14) {
        Bratu3D problem(args.n, current_lambda);
        const auto solution =
            newton_krylov(
                problem, previous, continuation_params, provider);

        if (!solution.converged) {
            current_step *= 0.5;
            if (current_step < args.min_step) {
                std::fprintf(
                    stderr,
                    "El corrector de continuación falló en lambda %.17g "
                    "y el paso cayó por debajo de %.17g. "
                    "No se declara un fold.\n",
                    current_lambda, args.min_step);
                return 3;
            }
            current_lambda =
                accepted.empty()
                    ? args.lambda_start
                    : accepted.back().lambda + current_step;
            continue;
        }

        auto apply = [&](const DVec& v, DVec& out) {
            problem.jvp(solution.x, v, out);
        };

        const auto eig =
            numa::analysis::inverse_smallest_eigenpair(
                apply, *provider.mg, problem.size(),
                args.eigen_rtol, args.linear_rtol,
                30, 100, 2000);

        if (!eig.converged || !std::isfinite(eig.eigenvalue) ||
            eig.eigenvalue <= 0.0) {
            std::fprintf(
                stderr,
                "El estimador espectral no convergió de forma válida "
                "en lambda %.17g. No se declara un fold.\n",
                current_lambda);
            return 4;
        }

        Point point;
        point.lambda = current_lambda;
        point.spectral_gap = eig.eigenvalue;
        point.eigen_residual = eig.residual_rel;
        point.continuation_iters = solution.iters;
        point.n_F = solution.n_F;
        point.n_JVP =
            solution.n_JVP + eig.linear_matvecs;
        point.n_prec =
            solution.n_prec + eig.preconditioner_applications;
        accepted.push_back(point);

        if (!append_point(csv, point)) {
            std::fprintf(stderr, "No se pudo escribir %s.\n", args.out.c_str());
            return 2;
        }

        std::printf(
            "lambda=%.8f gap=%.6e eig_res=%.3e cont_it=%zu\n",
            point.lambda, point.spectral_gap,
            point.eigen_residual, point.continuation_iters);

        previous = solution.x;

        if (point.spectral_gap <= args.spectral_stop) break;

        if (solution.iters <= 5 && current_step < args.step) {
            current_step = std::min(args.step, 1.5 * current_step);
        }
        current_lambda += current_step;
    }

    double proxy = 0.0;
    if (estimate_fold_proxy(accepted, proxy)) {
        std::printf(
            "FOLD_PROXY_SPECTRAL=%.10f "
            "(ajuste lineal de gap^2 en los cuatro últimos puntos; "
            "no es freeze)\n",
            proxy);
    } else {
        std::printf(
            "FOLD_PROXY_SPECTRAL=NA "
            "(datos insuficientes o ajuste no admisible)\n");
    }

    std::printf(
        "CALIBRACION_SOLAMENTE. El fold no se define por fallo del solver. "
        "Archivo: %s\n",
        args.out.c_str());
    return 0;
}
