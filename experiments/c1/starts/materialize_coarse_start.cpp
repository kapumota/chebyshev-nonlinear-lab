#define main embedded_calibrate_b6_arclength_main
#include "../../../apps/calibrate_b6_arclength.cpp"
#undef main

#include <numa/problems/burgers2d.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

bool write_vector(
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

DVec primal_from_state(
    const DVec& state,
    std::size_t primal_size)
{
    return DVec(
        state.begin(),
        state.begin()
            + static_cast<std::ptrdiff_t>(
                primal_size));
}

void interpolation_terms(
    std::size_t fine_index,
    std::size_t coarse_size,
    std::size_t coarse_indices[2],
    double weights[2],
    int& count)
{
    count = 0;

    if (fine_index % 2 == 1) {
        coarse_indices[count] =
            (fine_index - 1) / 2;

        weights[count] = 1.0;
        ++count;
        return;
    }

    const std::size_t right =
        fine_index / 2;

    if (right >= 1) {
        coarse_indices[count] =
            right - 1;

        weights[count] = 0.5;
        ++count;
    }

    if (right < coarse_size) {
        coarse_indices[count] =
            right;

        weights[count] = 0.5;
        ++count;
    }
}

std::size_t index3(
    std::size_t n,
    std::size_t i,
    std::size_t j,
    std::size_t k)
{
    return (i * n + j) * n + k;
}

DVec prolong3d(
    const DVec& coarse,
    std::size_t nc,
    std::size_t nf)
{
    if (nf != 2 * nc + 1 ||
        coarse.size() != nc * nc * nc) {

        throw std::runtime_error(
            "Dimensiones incompatibles en "
            "prolongación 3D");
    }

    DVec fine(
        nf * nf * nf,
        0.0);

    for (std::size_t i = 0;
         i < nf;
         ++i) {

        std::size_t ci[2];
        double wi[2];
        int ni = 0;

        interpolation_terms(
            i,
            nc,
            ci,
            wi,
            ni);

        for (std::size_t j = 0;
             j < nf;
             ++j) {

            std::size_t cj[2];
            double wj[2];
            int nj = 0;

            interpolation_terms(
                j,
                nc,
                cj,
                wj,
                nj);

            for (std::size_t k = 0;
                 k < nf;
                 ++k) {

                std::size_t ck[2];
                double wk[2];
                int nk = 0;

                interpolation_terms(
                    k,
                    nc,
                    ck,
                    wk,
                    nk);

                double sum = 0.0;

                for (int a = 0;
                     a < ni;
                     ++a) {

                    for (int b = 0;
                         b < nj;
                         ++b) {

                        for (int c = 0;
                             c < nk;
                             ++c) {

                            sum +=
                                wi[a]
                                * wj[b]
                                * wk[c]
                                * coarse[
                                    index3(
                                        nc,
                                        ci[a],
                                        cj[b],
                                        ck[c])];
                        }
                    }
                }

                fine[
                    index3(
                        nf,
                        i,
                        j,
                        k)] = sum;
            }
        }
    }

    return fine;
}

DVec prolong2d_components(
    const DVec& coarse,
    std::size_t nc,
    std::size_t nf)
{
    const std::size_t coarse_nn =
        nc * nc;

    const std::size_t fine_nn =
        nf * nf;

    if (nf != 2 * nc + 1 ||
        coarse.size() !=
            2 * coarse_nn) {

        throw std::runtime_error(
            "Dimensiones incompatibles en "
            "prolongación 2D");
    }

    DVec fine(
        2 * fine_nn,
        0.0);

    for (std::size_t component = 0;
         component < 2;
         ++component) {

        const std::size_t coarse_offset =
            component * coarse_nn;

        const std::size_t fine_offset =
            component * fine_nn;

        for (std::size_t i = 0;
             i < nf;
             ++i) {

            std::size_t ci[2];
            double wi[2];
            int ni = 0;

            interpolation_terms(
                i,
                nc,
                ci,
                wi,
                ni);

            for (std::size_t j = 0;
                 j < nf;
                 ++j) {

                std::size_t cj[2];
                double wj[2];
                int nj = 0;

                interpolation_terms(
                    j,
                    nc,
                    cj,
                    wj,
                    nj);

                double sum = 0.0;

                for (int a = 0;
                     a < ni;
                     ++a) {

                    for (int b = 0;
                         b < nj;
                         ++b) {

                        sum +=
                            wi[a]
                            * wj[b]
                            * coarse[
                                coarse_offset
                                + ci[a] * nc
                                + cj[b]];
                    }
                }

                fine[
                    fine_offset
                    + i * nf
                    + j] = sum;
            }
        }
    }

    return fine;
}

bool certify_bratu_reference(
    std::size_t n,
    const DVec& x,
    double lambda,
    const PoissonMGProvider& provider,
    double& spectral_gap)
{
    Args gate_args;
    gate_args.n = n;
    gate_args.eigen_rtol = 1e-7;
    gate_args.linear_rtol = 1e-10;

    Point point;

    const bool ok =
        spectral_point(
            n,
            x,
            lambda,
            provider,
            gate_args,
            0,
            0,
            0,
            0,
            "coarse_reference",
            point);

    spectral_gap =
        ok
            ? point.spectral_gap
            : std::numeric_limits<double>
                  ::quiet_NaN();

    return ok &&
           spectral_gap > 0.0;
}

int materialize_b1(
    std::size_t fine_n,
    double target,
    const std::string& output,
    const std::string& label)
{
    if (fine_n <= 1 ||
        fine_n % 2 == 0 ||
        !(target > 0.0) ||
        !std::isfinite(target)) {

        std::cerr
            << "Argumentos B1 inválidos\n";
        return 2;
    }

    const std::size_t coarse_n =
        (fine_n - 1) / 2;

    const std::size_t primal_size =
        coarse_n
        * coarse_n
        * coarse_n;

    constexpr double nominal_step =
        0.125;

    constexpr double min_step =
        1e-4;

    PoissonMGProvider provider(
        coarse_n);

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
    bool fallback_required = false;

    double current_step =
        nominal_step;

    while (latest_lambda < target) {
        const double next_lambda =
            std::min(
                latest_lambda
                    + current_step,
                target);

        Bratu3D problem(
            coarse_n,
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
            * std::numeric_limits<double>
                  ::epsilon()
            * std::max(
                1.0,
                std::abs(target));

        if (std::abs(
                next_lambda - target)
            <= tolerance) {

            double reference_gap = 0.0;

            if (!(latest_lambda > 0.0) ||
                !certify_bratu_reference(
                    coarse_n,
                    latest_x,
                    latest_lambda,
                    provider,
                    reference_gap)) {

                std::cerr
                    << "Falló certificación "
                    << "de referencia B1\n";
                return 5;
            }

            const DVec fine =
                prolong3d(
                    latest_x,
                    coarse_n,
                    fine_n);

            if (!write_vector(
                    output,
                    fine)) {

                std::cerr
                    << "Falló escritura B1\n";
                return 4;
            }

            std::cout
                << std::setprecision(17)
                << "family=B1\n"
                << "label=" << label << '\n'
                << "fine_n=" << fine_n << '\n'
                << "coarse_n="
                << coarse_n << '\n'
                << "lambda_target="
                << target << '\n'
                << "reference_lambda="
                << latest_lambda << '\n'
                << "reference_spectral_gap="
                << reference_gap << '\n'
                << "mode=natural\n"
                << "N=" << fine.size() << '\n'
                << "COARSE_START=PASS\n";

            return 0;
        }

        older_x = latest_x;
        older_lambda =
            latest_lambda;

        have_older = true;

        latest_x = solution.x;
        latest_lambda =
            next_lambda;

        if (solution.iters <= 5 &&
            current_step
                < nominal_step) {

            current_step =
                std::min(
                    nominal_step,
                    1.5
                        * current_step);
        }
    }

    if (!fallback_required ||
        !have_older) {

        std::cerr
            << "No hay estados suficientes "
            << "para pseudo-arclength B1\n";
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
            << "Secante inicial inválida\n";
        return 5;
    }

    double arc_step =
        secant_norm;

    if (!(arc_step > arc_min_step) ||
        !std::isfinite(arc_step)) {

        std::cerr
            << "Paso pseudo-arclength "
            << "inicial inválido\n";
        return 5;
    }

    Args gate_args;
    gate_args.n = coarse_n;
    gate_args.eigen_rtol = 1e-7;
    gate_args.linear_rtol = 1e-10;

    for (std::size_t attempt = 0;
         attempt < arc_max_attempts;
         ++attempt) {

        double local_secant = 0.0;

        if (!normalized_secant(
                y_previous,
                y_current,
                primal_size,
                tangent,
                local_secant)) {

            std::cerr
                << "Secante "
                << "pseudo-arclength "
                << "degenerada\n";
            return 5;
        }

        DVec predictor =
            y_current;

        for (std::size_t i = 0;
             i <= primal_size;
             ++i) {

            predictor[i] +=
                arc_step
                * tangent[i];
        }

        const auto corrected =
            correct_arclength(
                coarse_n,
                tangent,
                predictor,
                provider,
                1e-10);

        if (!corrected.converged) {
            arc_step *= 0.5;

            if (arc_step <
                arc_min_step) {

                std::cerr
                    << "Pseudo-arclength "
                    << "agotó arc_min_step\n";
                return 5;
            }

            continue;
        }

        const double candidate_lambda =
            corrected.state[
                primal_size];

        if (!(candidate_lambda >
              y_current[
                  primal_size])) {

            arc_step *= 0.5;

            if (arc_step <
                arc_min_step) {

                std::cerr
                    << "No se pudo avanzar "
                    << "en lambda\n";
                return 5;
            }

            continue;
        }

        const DVec candidate_x =
            primal_from_state(
                corrected.state,
                primal_size);

        Point candidate_point;

        if (!spectral_point(
                coarse_n,
                candidate_x,
                candidate_lambda,
                provider,
                gate_args,
                corrected.iters,
                corrected.n_F,
                corrected.n_JVP,
                corrected.n_prec,
                "pseudo_arclength",
                candidate_point) ||
            !(candidate_point
                  .spectral_gap > 0.0)) {

            arc_step *= 0.5;

            if (arc_step <
                arc_min_step) {

                std::cerr
                    << "Falló gate "
                    << "espectral B1\n";
                return 5;
            }

            continue;
        }

        if (candidate_lambda >=
            target) {

            const double reference_lambda =
                y_current[
                    primal_size];

            const DVec reference =
                primal_from_state(
                    y_current,
                    primal_size);

            double reference_gap = 0.0;

            if (!certify_bratu_reference(
                    coarse_n,
                    reference,
                    reference_lambda,
                    provider,
                    reference_gap)) {

                std::cerr
                    << "Falló certificación "
                    << "de referencia B1\n";
                return 5;
            }

            const DVec fine =
                prolong3d(
                    reference,
                    coarse_n,
                    fine_n);

            if (!write_vector(
                    output,
                    fine)) {

                std::cerr
                    << "Falló escritura B1\n";
                return 4;
            }

            std::cout
                << std::setprecision(17)
                << "family=B1\n"
                << "label=" << label << '\n'
                << "fine_n=" << fine_n << '\n'
                << "coarse_n="
                << coarse_n << '\n'
                << "lambda_target="
                << target << '\n'
                << "reference_lambda="
                << reference_lambda << '\n'
                << "candidate_lambda="
                << candidate_lambda << '\n'
                << "reference_spectral_gap="
                << reference_gap << '\n'
                << "candidate_spectral_gap="
                << candidate_point
                       .spectral_gap
                << '\n'
                << "mode=pseudo_arclength\n"
                << "N=" << fine.size() << '\n'
                << "COARSE_START=PASS\n";

            return 0;
        }

        y_previous =
            y_current;

        y_current =
            corrected.state;
    }

    std::cerr
        << "Pseudo-arclength agotó "
        << "intentos B1\n";

    return 5;
}

int materialize_b3a(
    std::size_t fine_n,
    double nu,
    double dt,
    double amp,
    const std::string& output,
    const std::string& label)
{
    if (fine_n <= 1 ||
        fine_n % 2 == 0 ||
        !(nu > 0.0) ||
        !(dt > 0.0) ||
        !std::isfinite(nu) ||
        !std::isfinite(dt) ||
        !std::isfinite(amp)) {

        std::cerr
            << "Argumentos B3a inválidos\n";
        return 2;
    }

    const std::size_t coarse_n =
        (fine_n - 1) / 2;

    Burgers2D coarse_problem(
        coarse_n,
        nu,
        dt,
        amp);

    const DVec fine =
        prolong2d_components(
            coarse_problem.exact(),
            coarse_n,
            fine_n);

    if (!write_vector(
            output,
            fine)) {

        std::cerr
            << "Falló escritura B3a\n";
        return 4;
    }

    std::cout
        << std::setprecision(17)
        << "family=B3a\n"
        << "label=" << label << '\n'
        << "fine_n=" << fine_n << '\n'
        << "coarse_n="
        << coarse_n << '\n'
        << "nu=" << nu << '\n'
        << "dt=" << dt << '\n'
        << "amp=" << amp << '\n'
        << "N=" << fine.size() << '\n'
        << "COARSE_START=PASS\n";

    return 0;
}

}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr
            << "Uso: materialize_coarse_start "
            << "B1|B3a ...\n";
        return 2;
    }

    const std::string family =
        argv[1];

    if (family == "B1") {
        if (argc != 6) {
            std::cerr
                << "Uso B1: B1 fine_n "
                << "lambda_target "
                << "salida.txt label\n";
            return 2;
        }

        return materialize_b1(
            static_cast<std::size_t>(
                std::stoull(argv[2])),
            std::stod(argv[3]),
            argv[4],
            argv[5]);
    }

    if (family == "B3a") {
        if (argc != 8) {
            std::cerr
                << "Uso B3a: B3a fine_n "
                << "nu dt amp "
                << "salida.txt label\n";
            return 2;
        }

        return materialize_b3a(
            static_cast<std::size_t>(
                std::stoull(argv[2])),
            std::stod(argv[3]),
            std::stod(argv[4]),
            std::stod(argv[5]),
            argv[6],
            argv[7]);
    }

    std::cerr
        << "Familia desconocida\n";

    return 2;
}
