#define main embedded_calibrate_b6_arclength_main
#include "../../../apps/calibrate_b6_arclength.cpp"
#undef main

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

bool write_start(
    const std::string& path,
    const DVec& values)
{
    const std::filesystem::path p(path);

    if (!p.parent_path().empty()) {
        std::filesystem::create_directories(
            p.parent_path());
    }

    std::ofstream out(path);

    if (!out) {
        return false;
    }

    out << std::setprecision(
        std::numeric_limits<double>::max_digits10);

    for (double value : values) {
        out << value << '\n';
    }

    return static_cast<bool>(out);
}

DVec primal_part(
    const DVec& state,
    std::size_t n)
{
    return DVec(
        state.begin(),
        state.begin()
            + static_cast<std::ptrdiff_t>(n));
}

}

int main(int argc, char** argv)
{
    if (argc != 7) {
        std::cerr
            << "Uso: materialize_bratu_pair_arclength "
            << "n lambda_fold rho "
            << "reference.txt continuation.txt label\n";
        return 2;
    }

    const std::size_t n =
        static_cast<std::size_t>(
            std::stoull(argv[1]));

    const double lambda_fold =
        std::stod(argv[2]);

    const double rho =
        std::stod(argv[3]);

    const std::string reference_path = argv[4];
    const std::string continuation_path = argv[5];
    const std::string label = argv[6];

    if (n == 0 ||
        !(lambda_fold > 0.0) ||
        !(rho > 0.0 && rho < 1.0)) {

        std::cerr << "Argumentos inválidos\n";
        return 2;
    }

    const std::size_t primal_size =
        n * n * n;

    const double target =
        rho * lambda_fold;

    constexpr double nominal_step = 0.125;
    constexpr double min_step = 1e-4;

    PoissonMGProvider provider(n);

    KrylovParams params;
    params.rtol_nl = 1e-11;
    params.max_nl = 100;
    params.eta_init = 1e-8;
    params.eta_max = 1e-6;
    params.gmres_restart = 100;
    params.gmres_max_iters = 2000;

    DVec latest_x(
        primal_size,
        0.0);

    DVec older_x;

    double latest_lambda = 0.0;
    double older_lambda = 0.0;

    bool have_older = false;

    double current_step = nominal_step;
    bool fallback_required = false;

    while (latest_lambda < target) {
        const double next_lambda =
            std::min(
                latest_lambda + current_step,
                target);

        Bratu3D problem(
            n,
            next_lambda);

        const auto solution =
            newton_krylov(
                problem,
                latest_x,
                params,
                provider);

        if (!solution.converged) {
            current_step *= 0.5;

            if (current_step < min_step) {
                fallback_required = true;
                break;
            }

            continue;
        }

        const double tolerance =
            32.0
            * std::numeric_limits<double>::epsilon()
            * std::max(1.0, std::abs(target));

        if (std::abs(
                next_lambda - target)
            <= tolerance) {

            if (!write_start(
                    reference_path,
                    latest_x) ||
                !write_start(
                    continuation_path,
                    latest_x)) {

                std::cerr
                    << "Falló escritura\n";
                return 4;
            }

            std::cout
                << std::setprecision(17)
                << "label=" << label << '\n'
                << "n=" << n << '\n'
                << "rho=" << rho << '\n'
                << "lambda_target="
                << target << '\n'
                << "continuation_lambda="
                << latest_lambda << '\n'
                << "mode=natural\n"
                << "BRATU_PAIR_ARCLENGTH=PASS\n";

            return 0;
        }

        older_x = latest_x;
        older_lambda = latest_lambda;
        have_older = true;

        latest_x = solution.x;
        latest_lambda = next_lambda;

        if (solution.iters <= 5 &&
            current_step < nominal_step) {

            current_step =
                std::min(
                    nominal_step,
                    1.5 * current_step);
        }
    }

    if (!fallback_required ||
        !have_older) {

        std::cerr
            << "No hay estados suficientes "
            << "para pseudo-arclength\n";
        return 5;
    }

    DVec y_previous =
        augmented_state(
            older_x,
            older_lambda);

    DVec y_current =
        augmented_state(
            latest_x,
            latest_lambda);

    DVec tangent;
    double secant_norm = 0.0;

    if (!normalized_secant(
            y_previous,
            y_current,
            primal_size,
            tangent,
            secant_norm)) {

        std::cerr
            << "No se pudo construir "
            << "la secante inicial\n";
        return 5;
    }

    double arc_step = secant_norm;

    if (!(arc_step > arc_min_step) ||
        !std::isfinite(arc_step)) {

        std::cerr
            << "Paso pseudo-arclength inválido\n";
        return 5;
    }

    Args gate_args;
    gate_args.n = n;
    gate_args.eigen_rtol = 1e-7;
    gate_args.linear_rtol = 1e-10;

    for (std::size_t attempt = 0;
         attempt < arc_max_attempts;
         ++attempt) {

        DVec predictor = y_current;

        for (std::size_t i = 0;
             i <= primal_size;
             ++i) {

            predictor[i] +=
                arc_step * tangent[i];
        }

        const auto corrected =
            correct_arclength(
                n,
                tangent,
                predictor,
                provider,
                1e-10);

        if (!corrected.converged) {
            arc_step *= 0.5;

            if (arc_step < arc_min_step) {
                std::cerr
                    << "Pseudo-arclength agotó "
                    << "arc_min_step\n";
                return 5;
            }

            continue;
        }

        const double candidate_lambda =
            corrected.state[primal_size];

        if (!(candidate_lambda >
              y_current[primal_size]) ||
            candidate_lambda >= lambda_fold) {

            arc_step *= 0.5;

            if (arc_step < arc_min_step) {
                std::cerr
                    << "No se pudo permanecer "
                    << "en la rama inferior\n";
                return 5;
            }

            continue;
        }

        const DVec candidate_x =
            primal_part(
                corrected.state,
                primal_size);

        Point point;

        if (!spectral_point(
                n,
                candidate_x,
                candidate_lambda,
                provider,
                gate_args,
                corrected.iters,
                corrected.n_F,
                corrected.n_JVP,
                corrected.n_prec,
                "pseudo_arclength",
                point) ||
            !(point.spectral_gap > 0.0)) {

            arc_step *= 0.5;

            if (arc_step < arc_min_step) {
                std::cerr
                    << "Falló gate espectral\n";
                return 5;
            }

            continue;
        }

        if (candidate_lambda >= target) {
            const double reference_lambda =
                y_current[primal_size];

            const DVec reference =
                primal_part(
                    y_current,
                    primal_size);

            if (!write_start(
                    reference_path,
                    reference) ||
                !write_start(
                    continuation_path,
                    reference)) {

                std::cerr
                    << "Falló escritura\n";
                return 4;
            }

            std::cout
                << std::setprecision(17)
                << "label=" << label << '\n'
                << "n=" << n << '\n'
                << "rho=" << rho << '\n'
                << "lambda_target="
                << target << '\n'
                << "reference_lambda="
                << reference_lambda << '\n'
                << "candidate_lambda="
                << candidate_lambda << '\n'
                << "spectral_gap="
                << point.spectral_gap << '\n'
                << "mode=pseudo_arclength\n"
                << "BRATU_PAIR_ARCLENGTH=PASS\n";

            return 0;
        }

        y_previous = y_current;
        y_current = corrected.state;
    }

    std::cerr
        << "Pseudo-arclength agotó intentos\n";

    return 5;
}
