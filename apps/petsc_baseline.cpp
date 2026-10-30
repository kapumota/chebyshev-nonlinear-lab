// Baseline externo PETSc para C1. Esta primera entrega cubre B1 Bratu
// en un proceso, con el mismo F, Jv y Jacobi del código interno.

#include <numa/instrument/csv.hpp>
#include <numa/operators/operator.hpp>
#include <numa/problems/bratu3d.hpp>

#include <petscsnes.h>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
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
    std::string out = "petsc_campaign.csv";
    std::string commit;
    std::size_t n = 7;
    double param = 0.5;
    double rtol = 1e-8;
    std::size_t max_nl = 50;
    std::size_t restart = 30;
    std::size_t maxlin = 200;
    int warmup = 0;
    int reps = 1;
};

struct PetscContext {
    PetscContext(std::size_t n, double lambda, bool jacobi)
        : problem(n, lambda), x_current(problem.size(), 0.0), use_jacobi(jacobi) {}

    Bratu3D problem;
    DVec x_current;
    DVec jacobi_diag;
    bool use_jacobi = false;
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_prec = 0;
    double t_F = 0.0;
    double t_JVP = 0.0;
    double t_prec = 0.0;
    double t_setup = 0.0;
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

void usage(const char* program) {
    std::fprintf(
        stderr,
        "Uso: %s --problem bratu --family B1 --n N --param LAMBDA "
        "--precond jacobi|none --start natural --rtol TOL --max-nl N "
        "--restart N --maxlin N --warmup N --reps N --out FILE "
        "[opciones PETSc]\n",
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
            key == "--precond" || key == "--start" || key == "--out" ||
            key == "--commit" || key == "--n" || key == "--param" ||
            key == "--rtol" || key == "--max-nl" || key == "--restart" ||
            key == "--maxlin" || key == "--warmup" || key == "--reps";

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
        else if (key == "--out") args.out = value;
        else if (key == "--commit") args.commit = value;
        else if (key == "--n" && !parse_size(value, args.n)) return false;
        else if (key == "--param" && !parse_double(value, args.param)) return false;
        else if (key == "--rtol" && !parse_double(value, args.rtol)) return false;
        else if (key == "--max-nl" && !parse_size(value, args.max_nl)) return false;
        else if (key == "--restart" && !parse_size(value, args.restart)) return false;
        else if (key == "--maxlin" && !parse_size(value, args.maxlin)) return false;
        else if (key == "--warmup" && !parse_nonnegative_int(value, args.warmup)) return false;
        else if (key == "--reps" && !parse_positive_int(value, args.reps)) return false;
    }

    if (args.problem != "bratu" || args.family != "B1" || args.start != "natural") {
        std::fprintf(stderr, "0020b-v1 admite sólo B1 Bratu con start natural.\n");
        return false;
    }
    if (args.precond != "jacobi" && args.precond != "none") {
        std::fprintf(stderr, "B1 PETSc admite --precond jacobi o none.\n");
        return false;
    }
    if (!(args.param >= 0.0) || !(args.rtol > 0.0) || args.restart == 0 ||
        args.maxlin == 0 || args.max_nl == 0) {
        return false;
    }
    if (args.commit.empty()) {
        std::fprintf(stderr, "Falta NUMA_COMMIT o --commit para trazabilidad.\n");
        return false;
    }
    return true;
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
    ctx->problem.eval(xd, fd);
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
    ctx->problem.jvp(ctx->x_current, vd, yd);
    ctx->t_JVP += elapsed(start);
    ++ctx->n_JVP;
    PetscCall(dvec_to_vec(yd, y));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode form_jacobian(SNES, Vec x, Mat, Mat, void* raw) {
    auto* ctx = static_cast<PetscContext*>(raw);
    PetscFunctionBeginUser;
    PetscCall(vec_to_dvec(x, ctx->x_current));
    const auto start = Clock::now();
    if (ctx->use_jacobi) {
        ctx->problem.jacobi_diag(ctx->x_current, ctx->jacobi_diag);
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
    PetscCall(vec_to_dvec(r, rd));
    PetscCheck(
        ctx->jacobi_diag.size() == rd.size(),
        PETSC_COMM_SELF,
        PETSC_ERR_ARG_WRONGSTATE,
        "Jacobi PETSc se aplicó antes de fijar la diagonal");

    const auto start = Clock::now();
    zd.resize(rd.size());
    for (std::size_t i = 0; i < rd.size(); ++i) {
        zd[i] = rd[i] / ctx->jacobi_diag[i];
    }
    ctx->t_prec += elapsed(start);
    ++ctx->n_prec;
    PetscCall(dvec_to_vec(zd, z));
    PetscFunctionReturn(PETSC_SUCCESS);
}

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

void print_effective_config(const EffectiveConfig& cfg) {
    std::printf("PETSC_BASELINE_BEGIN\n");
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
    PetscContext ctx(args.n, args.param, args.precond == "jacobi");
    const PetscInt N = static_cast<PetscInt>(ctx.problem.size());
    Vec x = nullptr, f = nullptr;
    Mat J = nullptr;
    SNES snes = nullptr;
    KSP ksp = nullptr;
    PC pc = nullptr;
    SNESLineSearch line_search = nullptr;

    PetscFunctionBeginUser;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, N, &x));
    PetscCall(VecDuplicate(x, &f));
    PetscCall(VecSet(x, 0.0));
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

    // MatShell y el precondicionador de la aplicación son invariantes del
    // baseline; las opciones PETSc pueden cambiar SNES/KSP/tolerancias, no F/Jv.
    PetscCall(SNESSetJacobian(snes, J, J, form_jacobian, &ctx));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetPC(ksp, &pc));
    if (args.precond == "jacobi") {
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
    const std::string expected_pc = args.precond == "jacobi" ? PCSHELL : PCNONE;
    PetscCheck(
        cfg.pc_type == expected_pc,
        PETSC_COMM_SELF,
        PETSC_ERR_SUP,
        "PC efectivo inesperado: %s",
        cfg.pc_type.c_str());
    if (effective != nullptr) *effective = cfg;

    DVec x0(ctx.problem.size(), 0.0), F0;
    ctx.problem.eval(x0, F0);
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
        row.problem = "bratu";
        row.family = "B1";
        row.n = args.n;
        row.N = ctx.problem.size();
        row.param = args.param;
        row.instance_id = "lambda=" + std::to_string(args.param);
        row.method = "petsc";
        row.precond = args.precond;
        row.start_kind = "natural";
        row.start_id = "natural";
        row.distance = "natural";
        row.rho_x = 0.0;
        row.seed = 0;
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
            "result repetition=%d success=%d nonlinear_iters=%" PetscInt_FMT
            " linear_iters=%" PetscInt_FMT " n_F=%zu petsc_n_F=%" PetscInt_FMT
            " n_JVP=%zu n_prec=%zu residual_rel=%.9e t_total=%.9e\n",
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
        "0020b-v1 es un baseline PETSc CPU de un solo proceso");

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
    print_effective_config(effective);
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
