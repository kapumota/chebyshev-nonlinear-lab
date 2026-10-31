// Continuación pseudo-arclength para Bratu3D en C1.
//
// Esta herramienta se activa como contingencia preregistrada cuando la continuación
// natural no logra acercarse al fold sin violar min_step. El fallo del corrector
// nunca define el fold. La pseudo-arclength sólo permite aproximar la rama inferior
// hasta satisfacer el gate espectral fijado antes del freeze.

#include <numa/analysis/spectral.hpp>
#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/problems/bratu3d.hpp>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace numa::op;

namespace {

struct Args {
    std::size_t n = 31;
    double lambda_start = 0.0;
    double lambda_stop = 10.2;
    double step = 0.25;
    double min_step = 1e-4;
    double spectral_stop = 5e-2;
    double eigen_rtol = 1e-7;
    double linear_rtol = 1e-10;
    std::string out = "bratu_fold.csv";
};

constexpr double arc_min_step = 1e-6;
constexpr std::size_t arc_max_attempts = 80;

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
    std::string mode;
};

bool append_point(std::ofstream& out, const Point& p) {
    out << p.lambda << ','
        << p.spectral_gap << ','
        << p.eigen_residual << ','
        << p.continuation_iters << ','
        << p.n_F << ','
        << p.n_JVP << ','
        << p.n_prec << ','
        << p.mode << '\n';
    return static_cast<bool>(out);
}

double weighted_dot(const DVec& a, const DVec& b, std::size_t primal_size) {
    double value = 0.0;
    for (std::size_t i = 0; i < primal_size; ++i) {
        value += a[i] * b[i];
    }
    value /= static_cast<double>(primal_size);
    value += a[primal_size] * b[primal_size];
    return value;
}

bool normalized_secant(
    const DVec& previous,
    const DVec& current,
    std::size_t primal_size,
    DVec& tangent,
    double& secant_norm)
{
    if (previous.size() != primal_size + 1 ||
        current.size() != primal_size + 1) {
        return false;
    }

    tangent.resize(primal_size + 1);
    for (std::size_t i = 0; i <= primal_size; ++i) {
        tangent[i] = current[i] - previous[i];
    }

    const double norm2_weighted = weighted_dot(tangent, tangent, primal_size);
    if (!(norm2_weighted > 0.0) || !std::isfinite(norm2_weighted)) {
        return false;
    }

    secant_norm = std::sqrt(norm2_weighted);
    for (double& value : tangent) value /= secant_norm;

    if (tangent[primal_size] < 0.0) {
        for (double& value : tangent) value = -value;
    }
    return true;
}

struct ArcCorrectResult {
    DVec state;
    bool converged = false;
    std::size_t iters = 0;
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_prec = 0;
    double residual = std::numeric_limits<double>::infinity();
};

double arc_residual_norm(
    std::size_t grid,
    const DVec& state,
    const DVec& tangent,
    const DVec& predictor,
    DVec& fx,
    double& constraint)
{
    const std::size_t primal_size = grid * grid * grid;
    const double lambda = state[primal_size];
    DVec x(
        state.begin(),
        state.begin() + static_cast<std::ptrdiff_t>(primal_size));

    Bratu3D problem(grid, lambda);
    problem.eval(x, fx);

    DVec delta(primal_size + 1);
    for (std::size_t i = 0; i <= primal_size; ++i) {
        delta[i] = state[i] - predictor[i];
    }
    constraint = weighted_dot(tangent, delta, primal_size);

    const double rms_f = norm2(fx) / std::sqrt(static_cast<double>(primal_size));
    return std::hypot(rms_f, constraint);
}

ArcCorrectResult correct_arclength(
    std::size_t grid,
    const DVec& tangent,
    const DVec& predictor,
    const PoissonMGProvider& provider,
    double linear_rtol)
{
    const std::size_t primal_size = grid * grid * grid;
    constexpr std::size_t max_newton = 30;
    constexpr double residual_tol = 1e-10;
    constexpr double armijo = 1e-4;
    constexpr double min_line = 1e-6;

    ArcCorrectResult result;
    result.state = predictor;

    DVec fx;
    double constraint = 0.0;

    for (std::size_t iter = 0; iter < max_newton; ++iter) {
        const double residual = arc_residual_norm(
            grid, result.state, tangent, predictor, fx, constraint);
        ++result.n_F;
        result.residual = residual;

        if (!std::isfinite(residual)) return result;
        if (residual <= residual_tol) {
            result.converged = true;
            result.iters = iter;
            return result;
        }

        const double lambda = result.state[primal_size];
        DVec x(
            result.state.begin(),
            result.state.begin() + static_cast<std::ptrdiff_t>(primal_size));
        Bratu3D problem(grid, lambda);
        const auto precond = provider.make(x);
        auto apply = [&](const DVec& v, DVec& out) {
            problem.jvp(x, v, out);
        };

        DVec rhs_a = fx;
        for (double& value : rhs_a) value = -value;
        DVec a(primal_size, 0.0);
        const auto solve_a = gmres(
            apply, precond, rhs_a, a, linear_rtol, 100, 2000);
        result.n_JVP += solve_a.matvecs;
        result.n_prec += solve_a.precs;
        if (!solve_a.converged || solve_a.breakdown) return result;

        DVec rhs_b(primal_size);
        for (std::size_t i = 0; i < primal_size; ++i) {
            rhs_b[i] = std::exp(x[i]);
        }
        DVec b(primal_size, 0.0);
        const auto solve_b = gmres(
            apply, precond, rhs_b, b, linear_rtol, 100, 2000);
        result.n_JVP += solve_b.matvecs;
        result.n_prec += solve_b.precs;
        if (!solve_b.converged || solve_b.breakdown) return result;

        double tangent_a = 0.0;
        double tangent_b = 0.0;
        for (std::size_t i = 0; i < primal_size; ++i) {
            tangent_a += tangent[i] * a[i];
            tangent_b += tangent[i] * b[i];
        }
        tangent_a /= static_cast<double>(primal_size);
        tangent_b /= static_cast<double>(primal_size);

        const double denominator = tangent[primal_size] + tangent_b;
        if (!std::isfinite(denominator) || std::abs(denominator) < 1e-12) {
            return result;
        }

        const double delta_lambda =
            (-constraint - tangent_a) / denominator;
        if (!std::isfinite(delta_lambda)) return result;

        DVec direction(primal_size + 1);
        for (std::size_t i = 0; i < primal_size; ++i) {
            direction[i] = a[i] + b[i] * delta_lambda;
        }
        direction[primal_size] = delta_lambda;

        bool accepted = false;
        double line = 1.0;
        while (line >= min_line) {
            DVec trial = result.state;
            for (std::size_t i = 0; i <= primal_size; ++i) {
                trial[i] += line * direction[i];
            }

            DVec trial_fx;
            double trial_constraint = 0.0;
            const double trial_residual = arc_residual_norm(
                grid, trial, tangent, predictor, trial_fx, trial_constraint);
            ++result.n_F;

            if (std::isfinite(trial_residual) &&
                trial_residual <= (1.0 - armijo * line) * residual) {
                result.state = std::move(trial);
                result.residual = trial_residual;
                accepted = true;
                break;
            }
            line *= 0.5;
        }

        if (!accepted) return result;
        result.iters = iter + 1;
    }

    DVec final_fx;
    double final_constraint = 0.0;
    result.residual = arc_residual_norm(
        grid, result.state, tangent, predictor, final_fx, final_constraint);
    ++result.n_F;
    result.converged =
        std::isfinite(result.residual) && result.residual <= residual_tol;
    return result;
}

bool spectral_point(
    std::size_t n,
    const DVec& solution,
    double lambda,
    const PoissonMGProvider& provider,
    const Args& args,
    std::size_t continuation_iters,
    std::size_t n_F,
    std::size_t n_JVP,
    std::size_t n_prec,
    const std::string& mode,
    Point& point)
{
    Bratu3D problem(n, lambda);
    auto apply = [&](const DVec& v, DVec& out) {
        problem.jvp(solution, v, out);
    };

    const auto eig =
        numa::analysis::inverse_smallest_eigenpair(
            apply, *provider.mg, problem.size(),
            args.eigen_rtol, args.linear_rtol,
            30, 100, 2000);

    if (!eig.converged || !std::isfinite(eig.eigenvalue) ||
        !std::isfinite(eig.residual_rel)) {
        return false;
    }

    point.lambda = lambda;
    point.spectral_gap = eig.eigenvalue;
    point.eigen_residual = eig.residual_rel;
    point.continuation_iters = continuation_iters;
    point.n_F = n_F;
    point.n_JVP = n_JVP + eig.linear_matvecs;
    point.n_prec = n_prec + eig.preconditioner_applications;
    point.mode = mode;
    return true;
}

DVec augmented_state(const DVec& x, double lambda) {
    DVec y = x;
    y.push_back(lambda);
    return y;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) {
        usage(argv[0]);
        return 2;
    }

    const std::size_t primal_size = args.n * args.n * args.n;
    PoissonMGProvider provider(args.n);

    KrylovParams continuation_params;
    continuation_params.rtol_nl = 1e-11;
    continuation_params.max_nl = 100;
    continuation_params.eta_init = 1e-8;
    continuation_params.eta_max = 1e-6;
    continuation_params.gmres_restart = 100;
    continuation_params.gmres_max_iters = 2000;

    std::ofstream csv(args.out);
    if (!csv) {
        std::fprintf(stderr, "No se pudo abrir %s.\n", args.out.c_str());
        return 2;
    }
    csv << "lambda,spectral_gap,eigen_residual,continuation_iters,"
           "n_F,n_JVP,n_prec,continuation_mode\n";

    DVec latest_x(primal_size, 0.0);
    DVec older_x;
    double latest_lambda = std::numeric_limits<double>::quiet_NaN();
    double older_lambda = std::numeric_limits<double>::quiet_NaN();
    bool have_latest = false;
    bool have_older = false;

    double current_lambda = args.lambda_start;
    double current_step = args.step;
    bool natural_completed = false;
    bool fallback_required = false;

    while (current_lambda <= args.lambda_stop + 1e-14) {
        Bratu3D problem(args.n, current_lambda);
        const DVec guess = have_latest ? latest_x : DVec(primal_size, 0.0);
        const auto solution =
            newton_krylov(problem, guess, continuation_params, provider);

        if (!solution.converged) {
            current_step *= 0.5;
            if (current_step < args.min_step) {
                fallback_required = true;
                std::printf(
                    "PSEUDO_ARCLENGTH_ACTIVATED lambda_failed=%.17g "
                    "min_step=%.17g\n",
                    current_lambda, args.min_step);
                break;
            }
            current_lambda =
                have_latest
                    ? latest_lambda + current_step
                    : args.lambda_start;
            continue;
        }

        Point point;
        if (!spectral_point(
                args.n, solution.x, current_lambda, provider, args,
                solution.iters, solution.n_F, solution.n_JVP, solution.n_prec,
                "natural", point) ||
            point.spectral_gap <= 0.0) {
            std::fprintf(
                stderr,
                "El estimador espectral natural no produjo un autovalor positivo "
                "en lambda %.17g.\n",
                current_lambda);
            return 4;
        }

        if (!append_point(csv, point)) {
            std::fprintf(stderr, "No se pudo escribir %s.\n", args.out.c_str());
            return 2;
        }

        std::printf(
            "mode=natural lambda=%.10f gap=%.6e eig_res=%.3e cont_it=%zu\n",
            point.lambda, point.spectral_gap,
            point.eigen_residual, point.continuation_iters);

        if (have_latest) {
            older_x = latest_x;
            older_lambda = latest_lambda;
            have_older = true;
        }
        latest_x = solution.x;
        latest_lambda = current_lambda;
        have_latest = true;

        if (point.spectral_gap <= args.spectral_stop) {
            natural_completed = true;
            break;
        }

        if (solution.iters <= 5 && current_step < args.step) {
            current_step = std::min(args.step, 1.5 * current_step);
        }
        current_lambda = latest_lambda + current_step;
    }

    if (natural_completed) {
        std::printf(
            "BRATU_FOLD_CALIBRATION=PASS mode=natural "
            "Archivo: %s\n",
            args.out.c_str());
        return 0;
    }

    if (!fallback_required || !have_older || !have_latest) {
        std::fprintf(
            stderr,
            "No hay dos estados aceptados para iniciar pseudo-arclength.\n");
        return 5;
    }

    DVec y_previous = augmented_state(older_x, older_lambda);
    DVec y_current = augmented_state(latest_x, latest_lambda);
    DVec tangent;
    double secant_norm = 0.0;
    if (!normalized_secant(
            y_previous, y_current, primal_size, tangent, secant_norm)) {
        std::fprintf(stderr, "No se pudo construir la secante inicial.\n");
        return 5;
    }

    double arc_step = secant_norm;
    if (!(arc_step > arc_min_step) || !std::isfinite(arc_step)) {
        std::fprintf(stderr, "Paso pseudo-arclength inicial inválido.\n");
        return 5;
    }

    for (std::size_t attempt = 0; attempt < arc_max_attempts; ++attempt) {
        double local_secant = 0.0;
        if (!normalized_secant(
                y_previous, y_current, primal_size, tangent, local_secant)) {
            std::fprintf(stderr, "Secante pseudo-arclength degenerada.\n");
            return 5;
        }

        DVec predictor = y_current;
        for (std::size_t i = 0; i <= primal_size; ++i) {
            predictor[i] += arc_step * tangent[i];
        }

        const auto corrected =
            correct_arclength(
                args.n, tangent, predictor, provider, args.linear_rtol);

        if (!corrected.converged) {
            std::printf(
                "PSEUDO_ARCLENGTH_RETRY arc_step=%.6e residual=%.6e\n",
                arc_step, corrected.residual);
            arc_step *= 0.5;
            if (arc_step < arc_min_step) {
                std::fprintf(
                    stderr,
                    "Pseudo-arclength no convergió antes de arc_min_step %.17g.\n",
                    arc_min_step);
                return 5;
            }
            continue;
        }

        const double candidate_lambda = corrected.state[primal_size];

        // C1 congela el fold desde el lado de la rama inferior. Si el corrector
        // cruza el turning point, se reduce el paso y se reintenta desde el último
        // punto inferior aceptado. El cruce no se usa como definición del fold.
        if (!(candidate_lambda > y_current[primal_size]) ||
            candidate_lambda > args.lambda_stop) {
            arc_step *= 0.5;
            if (arc_step < arc_min_step) {
                std::fprintf(
                    stderr,
                    "No se pudo aproximar el fold desde la rama inferior "
                    "antes de arc_min_step %.17g.\n",
                    arc_min_step);
                return 5;
            }
            continue;
        }

        DVec candidate_x(
            corrected.state.begin(),
            corrected.state.begin() + static_cast<std::ptrdiff_t>(primal_size));

        Point point;
        if (!spectral_point(
                args.n, candidate_x, candidate_lambda, provider, args,
                corrected.iters, corrected.n_F, corrected.n_JVP,
                corrected.n_prec, "pseudo_arclength", point) ||
            point.spectral_gap <= 0.0) {
            arc_step *= 0.5;
            if (arc_step < arc_min_step) {
                std::fprintf(
                    stderr,
                    "El estimador espectral no pudo mantenerse en la rama "
                    "inferior antes de arc_min_step %.17g.\n",
                    arc_min_step);
                return 5;
            }
            continue;
        }

        if (!append_point(csv, point)) {
            std::fprintf(stderr, "No se pudo escribir %s.\n", args.out.c_str());
            return 2;
        }

        std::printf(
            "mode=pseudo_arclength lambda=%.10f gap=%.6e "
            "eig_res=%.3e cont_it=%zu arc_step=%.6e\n",
            point.lambda, point.spectral_gap,
            point.eigen_residual, point.continuation_iters, arc_step);

        y_previous = y_current;
        y_current = corrected.state;

        if (point.spectral_gap <= args.spectral_stop) {
            std::printf(
                "BRATU_FOLD_CALIBRATION=PASS mode=pseudo_arclength "
                "Archivo: %s\n",
                args.out.c_str());
            return 0;
        }
    }

    std::fprintf(
        stderr,
        "Pseudo-arclength alcanzó el máximo de intentos sin satisfacer "
        "spectral_stop.\n");
    return 5;
}
