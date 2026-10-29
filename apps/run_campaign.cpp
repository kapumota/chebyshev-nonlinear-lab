// Ejecutor reproducible para calibración y campaña C1.
//
// El ejecutor no calcula una "solución de referencia" con uno de los métodos
// comparados. Los puntos natural, continuation y coarse se suministran de forma
// explícita. Para perturbed se exige un vector de referencia externo y se usan
// los niveles rho_x preregistrados.
//
// Ejemplo:
//   NUMA_COMMIT=$(git rev-parse HEAD) OMP_NUM_THREADS=8 \
//   run_campaign --phase calibration --problem bratu --n 31 --param 8.0 \
//     --precond mg --methods nk,cheb,rec --start perturbed \
//     --start-file start.txt --dist near,mid,far --seeds 12 \
//     --warmup 1 --reps 7 --out calibration.csv

#include <numa/instrument/csv.hpp>
#include <numa/linalg/poisson_mg.hpp>
#include <numa/methods/krylov_newton.hpp>
#include <numa/operators/timing.hpp>
#include <numa/problems/bratu3d.hpp>
#include <numa/problems/burgers2d.hpp>
#include <numa/problems/chandrasekhar_op.hpp>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <vector>

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

using namespace numa::op;
using namespace numa::instrument;

namespace {

struct Args {
    std::string phase = "exploratory";
    std::string problem = "bratu";
    std::string family = "B1";
    std::string precond = "mg";
    std::string methods = "nk,cheb,rec";
    std::string start = "natural";
    std::string start_file;
    std::string start_id;
    std::string dist = "near,mid,far";
    std::string out = "campaign.csv";
    std::string commit;
    std::size_t n = 31;
    double param = 1.0;
    double dt = 1.0;
    double amp = 1.0;
    double rtol = 1e-8;
    std::size_t max_nl = 50;
    double eta_init = 0.5;
    double eta_max = 0.9;
    double gamma = 0.9;
    double alpha = 2.0;
    double eta2_factor = 10.0;
    double eta2_max = 0.5;
    unsigned seeds = 12;
    int warmup = 1;
    int reps = 7;
    std::size_t restart = 50;
    std::size_t maxlin = 500;
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

bool parse_unsigned(const char* text, unsigned& value) {
    std::size_t parsed = 0;
    if (!parse_size(text, parsed) ||
        parsed > std::numeric_limits<unsigned>::max()) {
        return false;
    }
    value = static_cast<unsigned>(parsed);
    return true;
}

bool parse_nonnegative_int(const char* text, int& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const long raw = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        raw < 0 || raw > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(raw);
    return true;
}

bool parse_positive_int(const char* text, int& value) {
    if (!parse_nonnegative_int(text, value)) return false;
    return value > 0;
}

bool parse_double(const char* text, double& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const double raw = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' ||
        !std::isfinite(raw)) {
        return false;
    }
    value = raw;
    return true;
}

std::vector<std::string> split(const std::string& text) {
    std::vector<std::string> values;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (!token.empty()) values.push_back(token);
    }
    return values;
}

bool valid_method(const std::string& method) {
    return method == "nk" || method == "cheb" || method == "rec";
}

bool valid_phase(const std::string& phase) {
    return phase == "exploratory" ||
           phase == "calibration" ||
           phase == "confirmation";
}

bool valid_start(const std::string& start) {
    return start == "natural" ||
           start == "perturbed" ||
           start == "continuation" ||
           start == "coarse";
}

double rho_for(const std::string& distance) {
    if (distance == "near") return 0.05;
    if (distance == "mid") return 0.20;
    if (distance == "far") return 0.50;
    return -1.0;
}

void usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s --phase exploratory|calibration|confirmation "
        "--problem bratu|chandra|burgers --n N --param P "
        "--precond mg|jacobi|none --methods nk,cheb,rec "
        "--start natural|perturbed|continuation|coarse "
        "[--dt DT] [--amp A] [--start-file vector.txt] [--start-id ID] [--dist near,mid,far] "
        "[--seeds 12] [--warmup 1] [--reps 7] "
        "[--rtol 1e-8] [--max-nl 50] [--eta-init 0.5] "
        "[--eta-max 0.9] [--gamma 0.9] [--alpha 2.0] "
        "[--eta2-factor 10] [--eta2-max 0.5] "
        "[--restart 50] [--maxlin 500] "
        "--out archivo.csv [--commit SHA]\n",
        program);
}

bool parse_args(int argc, char** argv, Args& args) {
    if ((argc - 1) % 2 != 0) return false;

    if (const char* env = std::getenv("NUMA_COMMIT");
        env != nullptr && *env != '\0') {
        args.commit = env;
    }

    for (int i = 1; i < argc; i += 2) {
        const std::string key = argv[i];
        const char* value = argv[i + 1];

        if (key == "--phase") args.phase = value;
        else if (key == "--problem") args.problem = value;
        else if (key == "--family") args.family = value;
        else if (key == "--precond") args.precond = value;
        else if (key == "--methods") args.methods = value;
        else if (key == "--start") args.start = value;
        else if (key == "--start-file") args.start_file = value;
        else if (key == "--start-id") args.start_id = value;
        else if (key == "--dist") args.dist = value;
        else if (key == "--out") args.out = value;
        else if (key == "--commit") args.commit = value;
        else if (key == "--n") {
            if (!parse_size(value, args.n)) return false;
        } else if (key == "--param") {
            if (!parse_double(value, args.param)) return false;
        } else if (key == "--dt") {
            if (!parse_double(value, args.dt)) return false;
        } else if (key == "--amp") {
            if (!parse_double(value, args.amp)) return false;
        } else if (key == "--rtol") {
            if (!parse_double(value, args.rtol)) return false;
        } else if (key == "--max-nl") {
            if (!parse_size(value, args.max_nl)) return false;
        } else if (key == "--eta-init") {
            if (!parse_double(value, args.eta_init)) return false;
        } else if (key == "--eta-max") {
            if (!parse_double(value, args.eta_max)) return false;
        } else if (key == "--gamma") {
            if (!parse_double(value, args.gamma)) return false;
        } else if (key == "--alpha") {
            if (!parse_double(value, args.alpha)) return false;
        } else if (key == "--eta2-factor") {
            if (!parse_double(value, args.eta2_factor)) return false;
        } else if (key == "--eta2-max") {
            if (!parse_double(value, args.eta2_max)) return false;
        } else if (key == "--seeds") {
            if (!parse_unsigned(value, args.seeds)) return false;
        } else if (key == "--warmup") {
            if (!parse_nonnegative_int(value, args.warmup)) return false;
        } else if (key == "--reps") {
            if (!parse_positive_int(value, args.reps)) return false;
        } else if (key == "--restart") {
            if (!parse_size(value, args.restart)) return false;
        } else if (key == "--maxlin") {
            if (!parse_size(value, args.maxlin)) return false;
        } else {
            return false;
        }
    }

    if (!valid_phase(args.phase) || !valid_start(args.start) ||
        args.out.empty() || args.commit.empty() ||
        !(args.rtol > 0.0) || args.max_nl == 0 ||
        !(args.eta_init > 0.0) || !(args.eta_init <= args.eta_max) ||
        !(args.eta_max > 0.0 && args.eta_max < 1.0) ||
        !(args.gamma > 0.0) || !(args.alpha > 0.0) ||
        !(args.eta2_factor > 0.0) ||
        !(args.eta2_max > 0.0 && args.eta2_max < 1.0) ||
        args.restart == 0 || args.maxlin == 0) {
        return false;
    }

    if (args.problem != "bratu" && args.problem != "chandra" &&
        args.problem != "burgers") return false;
    if (args.problem == "burgers" && (!(args.param > 0.0) || !(args.dt > 0.0))) return false;
    if (args.precond != "mg" &&
        args.precond != "jacobi" &&
        args.precond != "none") {
        return false;
    }

    const auto methods = split(args.methods);
    if (methods.empty()) return false;
    for (const auto& method : methods) {
        if (!valid_method(method)) return false;
    }

    if (args.start == "perturbed") {
        if (args.start_file.empty() || args.start_id.empty() || args.seeds == 0) return false;
        for (const auto& distance : split(args.dist)) {
            if (rho_for(distance) <= 0.0) return false;
        }
    } else if ((args.start == "continuation" || args.start == "coarse") &&
               (args.start_file.empty() || args.start_id.empty())) {
        return false;
    }

    return true;
}

bool load_vector(
    const std::string& path,
    std::size_t expected,
    DVec& values)
{
    std::ifstream input(path);
    if (!input) return false;

    values.clear();
    double value = 0.0;
    while (input >> value) {
        if (!std::isfinite(value)) return false;
        values.push_back(value);
    }
    return input.eof() && values.size() == expected;
}

DVec perturb(
    const DVec& reference,
    double rho,
    unsigned seed)
{
    std::mt19937_64 generator(
        static_cast<unsigned long long>(seed) * 7919ULL + 13ULL);
    std::normal_distribution<double> normal(0.0, 1.0);

    DVec direction(reference.size());
    for (double& value : direction) value = normal(generator);

    const double direction_norm = norm2(direction);
    const double reference_norm = norm2(reference);
    if (!(direction_norm > 0.0) || !(reference_norm > 0.0)) {
        return {};
    }

    const double scale = rho * reference_norm / direction_norm;
    DVec result = reference;
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] += scale * direction[i];
    }
    return result;
}

int thread_count() {
#ifdef NUMA_HAS_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

long max_rss_kb() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
    return usage.ru_maxrss;
}

template <class Problem, class ProviderFactory>
bool run_method(
    const Problem& problem,
    const DVec& x0,
    const Args& args,
    const std::string& method,
    const std::string& distance,
    double rho_x,
    unsigned seed,
    ProviderFactory&& provider_factory)
{
    auto solve_once = [&](int repetition, bool record) -> bool {
        const auto total_start = timing_detail::Clock::now();

        const auto setup_start = timing_detail::Clock::now();
        auto provider = provider_factory();
        const double static_setup =
            timing_detail::elapsed(setup_start);

        Timed<Problem> timed_problem(problem);
        TimedPrecond<std::decay_t<decltype(provider)>>
            timed_provider(provider, timed_problem.timers());

        KrylovParams params;
        params.rtol_nl = args.rtol;
        params.max_nl = args.max_nl;
        params.eta_init = args.eta_init;
        params.eta_max = args.eta_max;
        params.gamma = args.gamma;
        params.alpha = args.alpha;
        params.eta2_factor = args.eta2_factor;
        params.eta2_max = args.eta2_max;
        params.gmres_restart = args.restart;
        params.gmres_max_iters = args.maxlin;
        params.recycle = method == "rec";

        KrylovResult result;
        if (method == "nk") {
            result = newton_krylov(
                timed_problem, x0, params, timed_provider);
        } else {
            result = chebyshev_inexact(
                timed_problem, x0, params, timed_provider);
        }

        const double total =
            timing_detail::elapsed(total_start);

        if (!record) return true;

        CsvRow row;
        row.phase = args.phase;
        row.problem = args.problem;
        row.family = args.family;
        row.n = args.n;
        row.N = problem.size();
        row.param = args.param;
        if (args.problem == "bratu") {
            row.instance_id = "lambda=" + std::to_string(args.param);
        } else if (args.problem == "chandra") {
            row.instance_id = "c=" + std::to_string(args.param);
        } else {
            row.instance_id = "nu=" + std::to_string(args.param) +
                              ";dt=" + std::to_string(args.dt) +
                              ";amp=" + std::to_string(args.amp);
        }
        row.method = method;
        row.precond = args.precond;
        row.start_kind = args.start;
        row.start_id =
            args.start == "natural" ? "natural" : args.start_id;
        row.distance = distance;
        row.rho_x = rho_x;
        row.seed = seed;
        row.success = result.converged;
        row.nonlinear_iters = result.iters;
        row.n_F = result.n_F;
        row.n_JVP = result.n_JVP;
        row.n_D2 = result.n_D2;
        row.n_lin = result.n_lin_iters;
        row.n_lin_failures = result.n_lin_failures;
        row.n_prec = result.n_prec;
        row.t_setup =
            static_setup + timed_problem.timers().t_setup_dynamic;
        row.t_F = timed_problem.timers().t_F;
        row.t_JVP = timed_problem.timers().t_JVP;
        row.t_D2 = timed_problem.timers().t_D2;
        row.t_prec = timed_problem.timers().t_prec;
        row.t_total = total;
        row.max_rss_kb = max_rss_kb();
        row.residual_rel =
            result.res0 > 0.0
                ? result.res / result.res0
                : result.res;
        row.threads = thread_count();
        row.repetition = repetition;
        row.rtol = args.rtol;
        row.max_nl = args.max_nl;
        row.eta_init = args.eta_init;
        row.eta_max = args.eta_max;
        row.gamma = args.gamma;
        row.alpha = args.alpha;
        row.eta2_factor = args.eta2_factor;
        row.eta2_max = args.eta2_max;
        row.restart = args.restart;
        row.maxlin = args.maxlin;
        row.commit = args.commit;

        if (!append_csv(args.out, row)) {
            std::fprintf(
                stderr,
                "No se pudo escribir %s.\n",
                args.out.c_str());
            return false;
        }
        return true;
    };

    for (int warmup = 0; warmup < args.warmup; ++warmup) {
        if (!solve_once(0, false)) return false;
    }
    for (int repetition = 1; repetition <= args.reps; ++repetition) {
        if (!solve_once(repetition, true)) return false;
    }
    return true;
}

template <class Problem, class ProviderFactory>
bool run_starts(
    const Problem& problem,
    const DVec& natural,
    const Args& args,
    ProviderFactory&& provider_factory)
{
    DVec reference;
    if (!args.start_file.empty() &&
        !load_vector(args.start_file, problem.size(), reference)) {
        std::fprintf(
            stderr,
            "No se pudo leer un vector de %zu elementos desde %s.\n",
            problem.size(), args.start_file.c_str());
        return false;
    }

    for (const auto& method : split(args.methods)) {
        if (args.start == "natural") {
            if (!run_method(
                    problem, natural, args, method,
                    "natural", 0.0, 0, provider_factory)) {
                return false;
            }
            continue;
        }

        if (args.start == "continuation" || args.start == "coarse") {
            if (!run_method(
                    problem, reference, args, method,
                    args.start, 0.0, 0, provider_factory)) {
                return false;
            }
            continue;
        }

        for (const auto& distance : split(args.dist)) {
            const double rho = rho_for(distance);
            for (unsigned seed = 1; seed <= args.seeds; ++seed) {
                const DVec x0 = perturb(reference, rho, seed);
                if (x0.empty()) {
                    std::fprintf(
                        stderr,
                        "No se pudo construir la perturbación %s, semilla %u.\n",
                        distance.c_str(), seed);
                    return false;
                }
                if (!run_method(
                        problem, x0, args, method,
                        distance, rho, seed, provider_factory)) {
                    return false;
                }
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) {
        usage(argv[0]);
        return 2;
    }

    if (args.phase == "confirmation") {
        std::printf(
            "CONFIRMATION solicitado. Verifique externamente que el preregistro "
            "esté frozen y que el árbol Git esté limpio.\n");
    }

    if (args.problem == "bratu") {
        if (!(args.param >= 0.0)) {
            std::fprintf(stderr, "lambda de Bratu debe ser no negativo.\n");
            return 2;
        }

        Bratu3D problem(args.n, args.param);
        const DVec natural(problem.size(), 0.0);

        if (args.precond == "mg") {
            return run_starts(
                       problem, natural, args,
                       [&] { return PoissonMGProvider(args.n); })
                ? 0 : 3;
        }
        if (args.precond == "jacobi") {
            return run_starts(
                       problem, natural, args,
                       [&] { return BratuJacobi{problem}; })
                ? 0 : 3;
        }
        return run_starts(
                   problem, natural, args,
                   [] { return NoPrecond{}; })
            ? 0 : 3;
    }

    if (args.problem == "burgers") {
        if (args.precond != "mg" && args.precond != "none") {
            std::fprintf(stderr, "Burgers admite --precond mg o none.\n");
            return 2;
        }
        Burgers2D problem(args.n, args.param, args.dt, args.amp);
        const DVec natural(problem.size(), 0.0);
        if (args.precond == "mg") {
            return run_starts(
                       problem, natural, args,
                       [&] { return BurgersMGProvider(args.n, args.param, args.dt); })
                ? 0 : 3;
        }
        return run_starts(
                   problem, natural, args,
                   [] { return NoPrecond{}; })
            ? 0 : 3;
    }

    if (!(args.param >= 0.0 && args.param <= 1.0)) {
        std::fprintf(stderr, "c de Chandrasekhar debe estar en [0,1].\n");
        return 2;
    }
    if (args.precond != "none") {
        std::fprintf(
            stderr,
            "Chandrasekhar usa --precond none en esta entrega.\n");
        return 2;
    }

    ChandrasekharOp problem(args.n, args.param);
    const DVec natural(problem.size(), 1.0);
    return run_starts(
               problem, natural, args,
               [] { return NoPrecond{}; })
        ? 0 : 3;
}
