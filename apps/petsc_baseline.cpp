// Baseline externo PETSc para C1. Usa el mismo F, Jv y precondicionador
// de la implementación interna para B1, B5, B3a y Bext.

#include <numa/instrument/csv.hpp>
#include <numa/operators/operator.hpp>
#include <numa/problems/bratu3d.hpp>
#include <numa/problems/brusselator1d.hpp>
#include <numa/problems/burgers2d.hpp>
#include <numa/problems/chandrasekhar_op.hpp>

#include <petscsnes.h>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <utility>
#include <vector>

#ifdef NUMA_HAS_OPENMP
#include <omp.h>
#endif

using numa::instrument::CsvRow;
using numa::instrument::append_csv;
using numa::op::Bratu3D;
using numa::op::BratuJacobi;
using numa::op::Brusselator1DStep;
using numa::op::BrusselatorJacobi;
using numa::op::Burgers2D;
using numa::op::BurgersMGProvider;
using numa::op::ChandrasekharOp;
using numa::op::DVec;
using numa::op::norm2;

namespace {

using Clock = std::chrono::steady_clock;

struct Args {
    std::string phase = "exploratory";
    std::string problem = "bratu";
    std::string family = "B1";
    std::string precond = "jacobi";
    std::string start = "natural";
    std::string start_file;
    std::string start_id;
    std::string distance = "natural";
    std::string out = "petsc_campaign.csv";
    std::string commit;
    std::size_t n = 7;
    double param = 0.5;
    double dt = 0.1;
    double amp = 1.0;
    double rho_x = 0.0;
    unsigned seed = 0;
    double rtol = 1e-8;
    std::size_t max_nl = 50;
    std::size_t restart = 30;
    std::size_t maxlin = 200;
    int warmup = 0;
    int reps = 1;
};

struct PetscContext {
    std::size_t size = 0;
    DVec initial;
    DVec x_current;
    bool use_precond = false;
    std::function<void(const DVec&, DVec&)> eval;
    std::function<void(const DVec&, const DVec&, DVec&)> jvp;
    std::function<void(const DVec&)> setup_precond;
    std::function<void(const DVec&, DVec&)> apply_precond;
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_prec = 0;
    double t_F = 0.0;
    double t_JVP = 0.0;
    double t_prec = 0.0;
    double t_setup = 0.0;
};

struct EffectiveConfig {
    std::string petsc_version;
    std::string snes_type;
    std::string linesearch_type;
    std::string ksp_type;
    std::string pc_type;
    PetscReal snes_atol = 0.0;
    PetscReal snes_rtol = 0.0;
    PetscReal snes_stol = 0.0;
    PetscInt snes_max_it = 0;
    PetscInt snes_max_funcs = 0;
    PetscReal ksp_rtol = 0.0;
    PetscReal ksp_atol = 0.0;
    PetscReal ksp_dtol = 0.0;
    PetscInt ksp_max_it = 0;
    PetscInt gmres_restart = 0;
};

double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

long max_rss_kb() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
    return usage.ru_maxrss;
}

int thread_count() {
#ifdef NUMA_HAS_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

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
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const unsigned long raw = std::strtoul(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        raw > std::numeric_limits<unsigned>::max()) {
        return false;
    }
    value = static_cast<unsigned>(raw);
    return true;
}

bool parse_positive_int(const char* text, int& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const long raw = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || raw <= 0 ||
        raw > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(raw);
    return true;
}

bool parse_nonnegative_int(const char* text, int& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const long raw = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || raw < 0 ||
        raw > std::numeric_limits<int>::max()) {
        return false;
    }
    value = static_cast<int>(raw);
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

bool valid_phase(const std::string& phase) {
    return phase == "exploratory" || phase == "calibration" ||
           phase == "confirmation";
}

bool valid_start(const std::string& start) {
    return start == "natural" || start == "perturbed" ||
           start == "continuation" || start == "coarse";
}

void usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s --phase exploratory|calibration|confirmation "
        "--problem bratu|chandra|burgers|brusselator "
        "--family B1|B5|B3a|Bext --n N --param P "
        "--precond jacobi|mg|none --start natural|perturbed|continuation|coarse "
        "[--dt DT] [--amp A] [--start-file FILE] [--start-id ID] "
        "[--distance natural|near|mid|far] [--rho-x RHO] [--seed N] "
        "--rtol TOL --max-nl N --restart N --maxlin N "
        "--warmup N --reps N --out FILE [opciones PETSc]\n",
        program);
}

bool parse_args(
    int argc,
    char** argv,
    Args& args,
    std::vector<std::string>& petsc_strings)
{
    if (const char* commit = std::getenv("NUMA_COMMIT")) args.commit = commit;
    petsc_strings.clear();
    petsc_strings.emplace_back(argv[0]);

    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--help") return false;

        const bool custom =
            key == "--phase" || key == "--problem" || key == "--family" ||
            key == "--precond" || key == "--start" || key == "--start-file" ||
            key == "--start-id" || key == "--distance" || key == "--out" ||
            key == "--commit" || key == "--n" || key == "--param" ||
            key == "--dt" || key == "--amp" || key == "--rho-x" ||
            key == "--seed" || key == "--rtol" || key == "--max-nl" ||
            key == "--restart" || key == "--maxlin" || key == "--warmup" ||
            key == "--reps";

        if (!custom) {
            petsc_strings.push_back(key);
            continue;
        }
        if (i + 1 >= argc) return false;
        const char* value = argv[++i];

        if (key == "--phase") args.phase = value;
        else if (key == "--problem") args.problem = value;
        else if (key == "--family") args.family = value;
        else if (key == "--precond") args.precond = value;
        else if (key == "--start") args.start = value;
        else if (key == "--start-file") args.start_file = value;
        else if (key == "--start-id") args.start_id = value;
        else if (key == "--distance") args.distance = value;
        else if (key == "--out") args.out = value;
        else if (key == "--commit") args.commit = value;
        else if (key == "--n" && !parse_size(value, args.n)) return false;
        else if (key == "--param" && !parse_double(value, args.param)) return false;
        else if (key == "--dt" && !parse_double(value, args.dt)) return false;
        else if (key == "--amp" && !parse_double(value, args.amp)) return false;
        else if (key == "--rho-x" && !parse_double(value, args.rho_x)) return false;
        else if (key == "--seed" && !parse_unsigned(value, args.seed)) return false;
        else if (key == "--rtol" && !parse_double(value, args.rtol)) return false;
        else if (key == "--max-nl" && !parse_size(value, args.max_nl)) return false;
        else if (key == "--restart" && !parse_size(value, args.restart)) return false;
        else if (key == "--maxlin" && !parse_size(value, args.maxlin)) return false;
        else if (key == "--warmup" && !parse_nonnegative_int(value, args.warmup)) return false;
        else if (key == "--reps" && !parse_positive_int(value, args.reps)) return false;
    }

    if (!valid_phase(args.phase) || !valid_start(args.start)) return false;
    if (!(args.rtol > 0.0) || args.restart == 0 || args.maxlin == 0 ||
        args.max_nl == 0) {
        return false;
    }
    if (args.commit.empty()) {
        std::fprintf(stderr, "Falta NUMA_COMMIT o --commit para trazabilidad.\n");
        return false;
    }

    if (args.start == "natural") {
        if (!args.start_file.empty()) {
            std::fprintf(stderr, "start natural no admite --start-file.\n");
            return false;
        }
        args.start_id = "natural";
        args.distance = "natural";
        args.rho_x = 0.0;
        args.seed = 0;
    } else if (args.start_file.empty() || args.start_id.empty()) {
        std::fprintf(stderr, "Los starts no naturales requieren --start-file y --start-id.\n");
        return false;
    }

    if (args.problem == "bratu") {
        return args.family == "B1" && args.param >= 0.0 &&
               (args.precond == "jacobi" || args.precond == "none");
    }
    if (args.problem == "chandra") {
        return args.family == "B5" && args.param >= 0.0 && args.param <= 1.0 &&
               args.precond == "none";
    }
    if (args.problem == "burgers") {
        return args.family == "B3a" && args.param > 0.0 && args.dt > 0.0 &&
               std::isfinite(args.amp) &&
               (args.precond == "mg" || args.precond == "none");
    }
    if (args.problem == "brusselator") {
        return args.family == "Bext" && args.param > 0.0 &&
               (args.precond == "jacobi" || args.precond == "none");
    }
    return false;
}

bool load_vector(const std::string& path, std::size_t expected, DVec& values) {
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

std::string instance_id(const Args& args) {
    char buffer[256];
    if (args.family == "B1") {
        std::snprintf(buffer, sizeof(buffer), "lambda=%.17g", args.param);
    } else if (args.family == "B5") {
        std::snprintf(buffer, sizeof(buffer), "c=%.17g", args.param);
    } else if (args.family == "B3a") {
        std::snprintf(
            buffer, sizeof(buffer), "nu=%.17g;dt=%.17g;amp=%.17g",
            args.param, args.dt, args.amp);
    } else {
        std::snprintf(
            buffer, sizeof(buffer), "dt=%.17g;A=1;B=3;alpha=0.02",
            args.param);
    }
    return buffer;
}

std::unique_ptr<PetscContext> make_context(const Args& args) {
    auto ctx = std::make_unique<PetscContext>();

    if (args.family == "B1") {
        auto problem = std::make_shared<Bratu3D>(args.n, args.param);
        ctx->size = problem->size();
        ctx->initial.assign(ctx->size, 0.0);
        ctx->eval = [problem](const DVec& x, DVec& out) { problem->eval(x, out); };
        ctx->jvp = [problem](const DVec& x, const DVec& v, DVec& out) {
            problem->jvp(x, v, out);
        };
        if (args.precond == "jacobi") {
            auto bound = std::make_shared<BratuJacobi::Bound>();
            ctx->use_precond = true;
            ctx->setup_precond = [problem, bound](const DVec& x) {
                *bound = BratuJacobi{*problem}.make(x);
            };
            ctx->apply_precond = [bound](const DVec& r, DVec& z) {
                (*bound)(r, z);
            };
        }
    } else if (args.family == "B5") {
        auto problem = std::make_shared<ChandrasekharOp>(args.n, args.param);
        ctx->size = problem->size();
        ctx->initial.assign(ctx->size, 1.0);
        ctx->eval = [problem](const DVec& x, DVec& out) { problem->eval(x, out); };
        ctx->jvp = [problem](const DVec& x, const DVec& v, DVec& out) {
            problem->jvp(x, v, out);
        };
    } else if (args.family == "B3a") {
        auto problem = std::make_shared<Burgers2D>(args.n, args.param, args.dt, args.amp);
        ctx->size = problem->size();
        ctx->initial.assign(ctx->size, 0.0);
        ctx->eval = [problem](const DVec& x, DVec& out) { problem->eval(x, out); };
        ctx->jvp = [problem](const DVec& x, const DVec& v, DVec& out) {
            problem->jvp(x, v, out);
        };
        if (args.precond == "mg") {
            auto provider = std::make_shared<BurgersMGProvider>(
                args.n, args.param, args.dt);
            auto bound = std::make_shared<BurgersMGProvider::Bound>();
            ctx->use_precond = true;
            ctx->setup_precond = [provider, bound](const DVec& x) {
                *bound = provider->make(x);
            };
            ctx->apply_precond = [bound](const DVec& r, DVec& z) {
                (*bound)(r, z);
            };
        }
    } else {
        auto problem = std::make_shared<Brusselator1DStep>(args.n, args.param);
        ctx->size = problem->size();
        ctx->initial = problem->previous();
        ctx->eval = [problem](const DVec& x, DVec& out) { problem->eval(x, out); };
        ctx->jvp = [problem](const DVec& x, const DVec& v, DVec& out) {
            problem->jvp(x, v, out);
        };
        if (args.precond == "jacobi") {
            auto bound = std::make_shared<BrusselatorJacobi::Bound>();
            ctx->use_precond = true;
            ctx->setup_precond = [problem, bound](const DVec& x) {
                *bound = BrusselatorJacobi{*problem}.make(x);
            };
            ctx->apply_precond = [bound](const DVec& r, DVec& z) {
                (*bound)(r, z);
            };
        }
    }

    if (args.start != "natural") {
        if (!load_vector(args.start_file, ctx->size, ctx->initial)) {
            throw std::runtime_error(
                "No se pudo leer el start materializado con la dimensión esperada");
        }
    }
    ctx->x_current = ctx->initial;
    return ctx;
}

PetscErrorCode vec_to_dvec(Vec source, DVec& target) {
    PetscInt n = 0;
    const PetscScalar* values = nullptr;
    PetscFunctionBeginUser;
    PetscCall(VecGetLocalSize(source, &n));
    target.resize(static_cast<std::size_t>(n));
    PetscCall(VecGetArrayRead(source, &values));
    for (PetscInt i = 0; i < n; ++i) {
        target[static_cast<std::size_t>(i)] = PetscRealPart(values[i]);
    }
    PetscCall(VecRestoreArrayRead(source, &values));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode dvec_to_vec(const DVec& source, Vec target) {
    PetscInt n = 0;
    PetscScalar* values = nullptr;
    PetscFunctionBeginUser;
    PetscCall(VecGetLocalSize(target, &n));
    PetscCheck(
        static_cast<std::size_t>(n) == source.size(),
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_SIZ,
        "Dimensión PETSc incompatible con el operador C1");
    PetscCall(VecGetArray(target, &values));
    for (PetscInt i = 0; i < n; ++i) {
        values[i] = source[static_cast<std::size_t>(i)];
    }
    PetscCall(VecRestoreArray(target, &values));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_function(SNES, Vec x, Vec f, void* raw) {
    auto* ctx = static_cast<PetscContext*>(raw);
    DVec xd, fd;
    PetscFunctionBeginUser;
    PetscCall(vec_to_dvec(x, xd));
    const auto start = Clock::now();
    ctx->eval(xd, fd);
    ctx->t_F += elapsed(start);
    ++ctx->n_F;
    PetscCall(dvec_to_vec(fd, f));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode mat_mult(Mat A, Vec v, Vec y) {
    void* raw = nullptr;
    DVec vd, yd;
    PetscFunctionBeginUser;
    PetscCall(MatShellGetContext(A, &raw));
    auto* ctx = static_cast<PetscContext*>(raw);
    PetscCall(vec_to_dvec(v, vd));
    const auto start = Clock::now();
    ctx->jvp(ctx->x_current, vd, yd);
    ctx->t_JVP += elapsed(start);
    ++ctx->n_JVP;
    PetscCall(dvec_to_vec(yd, y));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_jacobian(SNES, Vec x, Mat, Mat, void* raw) {
    auto* ctx = static_cast<PetscContext*>(raw);
    PetscFunctionBeginUser;
    PetscCall(vec_to_dvec(x, ctx->x_current));
    if (ctx->use_precond) {
        const auto start = Clock::now();
        ctx->setup_precond(ctx->x_current);
        ctx->t_setup += elapsed(start);
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode pc_apply(PC pc, Vec r, Vec z) {
    void* raw = nullptr;
    DVec rd, zd;
    PetscFunctionBeginUser;
    PetscCall(PCShellGetContext(pc, &raw));
    auto* ctx = static_cast<PetscContext*>(raw);
    PetscCheck(
        ctx->use_precond && static_cast<bool>(ctx->apply_precond),
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_WRONGSTATE,
        "PCSHELL se aplicó sin precondicionador de aplicación");
    PetscCall(vec_to_dvec(r, rd));
    const auto start = Clock::now();
    ctx->apply_precond(rd, zd);
    ctx->t_prec += elapsed(start);
    ++ctx->n_prec;
    PetscCall(dvec_to_vec(zd, z));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode read_effective_config(SNES snes, EffectiveConfig& out) {
    char version[256] = {};
    SNESType snes_type = nullptr;
    SNESLineSearch line_search = nullptr;
    SNESLineSearchType line_type = nullptr;
    KSP ksp = nullptr;
    KSPType ksp_type = nullptr;
    PC pc = nullptr;
    PCType pc_type = nullptr;

    PetscFunctionBeginUser;
    PetscCall(PetscGetVersion(version, sizeof(version)));
    out.petsc_version = version;
    PetscCall(SNESGetType(snes, &snes_type));
    PetscCall(SNESGetLineSearch(snes, &line_search));
    PetscCall(SNESLineSearchGetType(line_search, &line_type));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetType(ksp, &ksp_type));
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCGetType(pc, &pc_type));
    out.snes_type = snes_type != nullptr ? snes_type : "unset";
    out.linesearch_type = line_type != nullptr ? line_type : "unset";
    out.ksp_type = ksp_type != nullptr ? ksp_type : "unset";
    out.pc_type = pc_type != nullptr ? pc_type : "unset";
    PetscCall(SNESGetTolerances(
        snes,
        &out.snes_atol,
        &out.snes_rtol,
        &out.snes_stol,
        &out.snes_max_it,
        &out.snes_max_funcs));
    PetscCall(KSPGetTolerances(
        ksp,
        &out.ksp_rtol,
        &out.ksp_atol,
        &out.ksp_dtol,
        &out.ksp_max_it));
    if (out.ksp_type == KSPGMRES) {
        PetscCall(KSPGMRESGetRestart(ksp, &out.gmres_restart));
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

void print_effective_config(const Args& args, const EffectiveConfig& cfg) {
    std::printf("PETSC_BASELINE_BEGIN\n");
    std::printf("problem=%s\n", args.problem.c_str());
    std::printf("family=%s\n", args.family.c_str());
    std::printf("application_preconditioner=%s\n", args.precond.c_str());
    std::printf("start_kind=%s\n", args.start.c_str());
    std::printf("start_id=%s\n", args.start_id.c_str());
    std::printf("distance=%s\n", args.distance.c_str());
    std::printf("rho_x=%.17g\n", args.rho_x);
    std::printf("seed=%u\n", args.seed);
    std::printf("petsc_version=%s\n", cfg.petsc_version.c_str());
    std::printf("snes_type=%s\n", cfg.snes_type.c_str());
    std::printf("linesearch_type=%s\n", cfg.linesearch_type.c_str());
    std::printf("ksp_type=%s\n", cfg.ksp_type.c_str());
    std::printf("pc_type=%s\n", cfg.pc_type.c_str());
    std::printf("matrix_mode=MatShell analytic Jv\n");
    std::printf("snes_rtol=%.17g\n", static_cast<double>(cfg.snes_rtol));
    std::printf("snes_atol=%.17g\n", static_cast<double>(cfg.snes_atol));
    std::printf("snes_stol=%.17g\n", static_cast<double>(cfg.snes_stol));
    std::printf("snes_max_it=%" PetscInt_FMT "\n", cfg.snes_max_it);
    std::printf("snes_max_funcs=%" PetscInt_FMT "\n", cfg.snes_max_funcs);
    std::printf("ksp_rtol=%.17g\n", static_cast<double>(cfg.ksp_rtol));
    std::printf("ksp_atol=%.17g\n", static_cast<double>(cfg.ksp_atol));
    std::printf("ksp_dtol=%.17g\n", static_cast<double>(cfg.ksp_dtol));
    std::printf("ksp_max_it=%" PetscInt_FMT "\n", cfg.ksp_max_it);
    std::printf("gmres_restart=%" PetscInt_FMT "\n", cfg.gmres_restart);
    std::printf("PETSC_BASELINE_END\n");
}

PetscErrorCode solve_once(
    const Args& args,
    int repetition,
    bool record,
    EffectiveConfig* effective)
{
    std::unique_ptr<PetscContext> ctx_owner;
    try {
        ctx_owner = make_context(args);
    } catch (const std::exception& error) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "%s", error.what());
    }
    auto& ctx = *ctx_owner;
    const PetscInt N = static_cast<PetscInt>(ctx.size);
    Vec x = nullptr, f = nullptr;
    Mat J = nullptr;
    SNES snes = nullptr;
    KSP ksp = nullptr;
    PC pc = nullptr;
    SNESLineSearch line_search = nullptr;

    PetscFunctionBeginUser;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &x));
    PetscCall(VecDuplicate(x, &f));
    PetscCall(dvec_to_vec(ctx.initial, x));
    PetscCall(MatCreateShell(PETSC_COMM_SELF, N, N, N, N, &ctx, &J));
    PetscCall(MatShellSetOperation(
        J,
        MATOP_MULT,
        reinterpret_cast<void (*)(void)>(mat_mult)));

    PetscCall(SNESCreate(PETSC_COMM_SELF, &snes));
    PetscCall(SNESSetType(snes, SNESNEWTONLS));
    PetscCall(SNESSetFunction(snes, f, form_function, &ctx));
    PetscCall(SNESSetJacobian(snes, J, J, form_jacobian, &ctx));
    PetscCall(SNESSetTolerances(
        snes,
        PETSC_CURRENT,
        args.rtol,
        PETSC_CURRENT,
        static_cast<PetscInt>(args.max_nl),
        PETSC_CURRENT));

    PetscCall(SNESGetLineSearch(snes, &line_search));
    PetscCall(SNESLineSearchSetType(line_search, SNESLINESEARCHBT));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPSetType(ksp, KSPGMRES));
    PetscCall(KSPGMRESSetRestart(ksp, static_cast<PetscInt>(args.restart)));
    PetscCall(KSPSetTolerances(
        ksp,
        PETSC_CURRENT,
        PETSC_CURRENT,
        PETSC_CURRENT,
        static_cast<PetscInt>(args.maxlin)));

    PetscCall(SNESSetFromOptions(snes));

    // F, Jv y el precondicionador de aplicación no son sustituibles por
    // opciones PETSc. Sólo se permiten ajustes del solver/tolerancias.
    PetscCall(SNESSetJacobian(snes, J, J, form_jacobian, &ctx));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    if (ctx.use_precond) {
        PetscCall(PCSetType(pc, PCSHELL));
        PetscCall(PCShellSetContext(pc, &ctx));
        PetscCall(PCShellSetApply(pc, pc_apply));
    } else {
        PetscCall(PCSetType(pc, PCNONE));
    }

    EffectiveConfig cfg;
    PetscCall(read_effective_config(snes, cfg));
    PetscCheck(
        cfg.snes_type == SNESNEWTONLS,
        PETSC_COMM_SELF,
        PETSC_ERR_SUP,
        "El baseline C1 exige SNESNEWTONLS; tipo efectivo: %s",
        cfg.snes_type.c_str());
    PetscCheck(
        cfg.linesearch_type == SNESLINESEARCHBT,
        PETSC_COMM_SELF,
        PETSC_ERR_SUP,
        "El baseline C1 exige line search bt; tipo efectivo: %s",
        cfg.linesearch_type.c_str());
    PetscCheck(
        cfg.ksp_type == KSPGMRES,
        PETSC_COMM_SELF,
        PETSC_ERR_SUP,
        "El baseline C1 exige KSPGMRES; tipo efectivo: %s",
        cfg.ksp_type.c_str());
    const std::string expected_pc = ctx.use_precond ? PCSHELL : PCNONE;
    PetscCheck(
        cfg.pc_type == expected_pc,
        PETSC_COMM_SELF,
        PETSC_ERR_SUP,
        "PC efectivo inesperado: %s",
        cfg.pc_type.c_str());
    if (effective != nullptr) *effective = cfg;

    DVec F0;
    ctx.eval(ctx.initial, F0);
    const double res0 = norm2(F0);

    const auto total_start = Clock::now();
    PetscCall(SNESSolve(snes, nullptr, x));
    const double t_total = elapsed(total_start);

    SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
    PetscInt nonlinear_iters = 0;
    PetscInt linear_iters = 0;
    PetscInt linear_failures = 0;
    PetscInt petsc_function_evals = 0;
    PetscReal residual = 0.0;
    PetscCall(SNESGetConvergedReason(snes, &reason));
    PetscCall(SNESGetIterationNumber(snes, &nonlinear_iters));
    PetscCall(SNESGetLinearSolveIterations(snes, &linear_iters));
    PetscCall(SNESGetLinearSolveFailures(snes, &linear_failures));
    PetscCall(SNESGetNumberFunctionEvals(snes, &petsc_function_evals));
    PetscCall(SNESGetFunctionNorm(snes, &residual));

    if (record) {
        CsvRow row;
        row.phase = args.phase;
        row.problem = args.problem;
        row.family = args.family;
        row.n = args.n;
        row.N = ctx.size;
        row.param = args.param;
        row.instance_id = instance_id(args);
        row.method = "petsc";
        row.precond = args.precond;
        row.start_kind = args.start;
        row.start_id = args.start_id;
        row.distance = args.distance;
        row.rho_x = args.rho_x;
        row.seed = args.seed;
        row.success = reason > 0;
        row.nonlinear_iters = static_cast<std::size_t>(nonlinear_iters);
        row.n_F = ctx.n_F;
        row.n_JVP = ctx.n_JVP;
        row.n_D2 = 0;
        row.n_lin = static_cast<std::size_t>(linear_iters);
        row.n_lin_failures = static_cast<std::size_t>(linear_failures);
        row.n_prec = ctx.n_prec;
        row.t_setup = ctx.t_setup;
        row.t_F = ctx.t_F;
        row.t_JVP = ctx.t_JVP;
        row.t_D2 = 0.0;
        row.t_prec = ctx.t_prec;
        row.t_total = t_total;
        row.max_rss_kb = max_rss_kb();
        row.residual_rel = res0 > 0.0
            ? static_cast<double>(residual) / res0
            : static_cast<double>(residual);
        row.threads = thread_count();
        row.repetition = repetition;
        row.rtol = static_cast<double>(cfg.snes_rtol);
        row.max_nl = static_cast<std::size_t>(cfg.snes_max_it);
        row.eta_init = 0.0;
        row.eta_max = 0.0;
        row.gamma = 0.0;
        row.alpha = 0.0;
        row.eta2_factor = 0.0;
        row.eta2_max = 0.0;
        row.restart = static_cast<std::size_t>(cfg.gmres_restart);
        row.maxlin = static_cast<std::size_t>(cfg.ksp_max_it);
        row.commit = args.commit;
        PetscCheck(
            append_csv(args.out, row),
            PETSC_COMM_SELF,
            PETSC_ERR_FILE_WRITE,
            "No se pudo escribir %s",
            args.out.c_str());

        std::printf(
            "result family=%s repetition=%d success=%d nonlinear_iters=%" PetscInt_FMT
            " linear_iters=%" PetscInt_FMT " n_F=%zu petsc_n_F=%" PetscInt_FMT
            " n_JVP=%zu n_prec=%zu residual_rel=%.9e t_total=%.9e\n",
            args.family.c_str(),
            repetition,
            reason > 0 ? 1 : 0,
            nonlinear_iters,
            linear_iters,
            ctx.n_F,
            petsc_function_evals,
            ctx.n_JVP,
            ctx.n_prec,
            row.residual_rel,
            t_total);
    }

    PetscCall(SNESDestroy(&snes));
    PetscCall(MatDestroy(&J));
    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&x));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode run(const Args& args) {
    PetscMPIInt size = 0;
    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCheck(
        size == 1,
        PETSC_COMM_WORLD,
        PETSC_ERR_SUP,
        "0020c-v1 es un baseline PETSc CPU de un solo proceso");

    for (int i = 0; i < args.warmup; ++i) {
        PetscCall(solve_once(args, 0, false, nullptr));
    }

    EffectiveConfig effective;
    for (int repetition = 1; repetition <= args.reps; ++repetition) {
        PetscCall(solve_once(
            args,
            repetition,
            true,
            repetition == 1 ? &effective : nullptr));
    }
    print_effective_config(args, effective);
    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    std::vector<std::string> petsc_strings;
    if (!parse_args(argc, argv, args, petsc_strings)) {
        usage(argv[0]);
        return 2;
    }

    std::vector<char*> petsc_argv;
    petsc_argv.reserve(petsc_strings.size());
    for (auto& item : petsc_strings) petsc_argv.push_back(item.data());
    int petsc_argc = static_cast<int>(petsc_argv.size());
    char** petsc_argv_raw = petsc_argv.data();

    PetscErrorCode ierr = PetscInitialize(
        &petsc_argc,
        &petsc_argv_raw,
        nullptr,
        "Baseline externo PETSc de C1");
    if (ierr != PETSC_SUCCESS) return static_cast<int>(ierr);

    ierr = run(args);
    const PetscErrorCode finalize_ierr = PetscFinalize();
    if (ierr != PETSC_SUCCESS) return static_cast<int>(ierr);
    return finalize_ierr == PETSC_SUCCESS ? 0 : static_cast<int>(finalize_ierr);
}
