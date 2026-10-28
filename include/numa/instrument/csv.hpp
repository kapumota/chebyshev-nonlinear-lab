#pragma once

// Esquema CSV crudo para C1. No calcula W_eq: los costos primitivos se
// congelan en una fase posterior y W_eq se deriva sin reescribir datos crudos.

#include <cstdio>
#include <string>

namespace numa::instrument {

struct CsvRow {
    std::string phase;
    std::string problem;
    std::string family;
    std::size_t n = 0;
    std::size_t N = 0;
    double param = 0.0;
    std::string method;
    std::string precond;
    std::string start_kind;
    std::string start_id;
    std::string distance;
    double rho_x = 0.0;
    unsigned seed = 0;
    bool success = false;
    std::size_t nonlinear_iters = 0;
    std::size_t n_F = 0;
    std::size_t n_JVP = 0;
    std::size_t n_D2 = 0;
    std::size_t n_lin = 0;
    std::size_t n_lin_failures = 0;
    std::size_t n_prec = 0;
    double t_setup = 0.0;
    double t_F = 0.0;
    double t_JVP = 0.0;
    double t_D2 = 0.0;
    double t_prec = 0.0;
    double t_total = 0.0;
    long max_rss_kb = 0;
    double residual_rel = 0.0;
    int threads = 1;
    int repetition = 0;
    std::string commit;
};

inline const char* csv_header() {
    return
        "phase,problem,family,n,N,param,method,precond,start_kind,start_id,"
        "distance,rho_x,seed,success,nonlinear_iters,n_F,n_JVP,n_D2,n_lin,"
        "n_lin_failures,n_prec,"
        "t_setup,t_F,t_JVP,t_D2,t_prec,t_total,max_rss_kb,residual_rel,"
        "threads,repetition,commit";
}

inline std::string csv_line(const CsvRow& row) {
    char buffer[2048];
    std::snprintf(
        buffer, sizeof(buffer),
        "%s,%s,%s,%zu,%zu,%.17g,%s,%s,%s,%s,%s,%.17g,%u,%d,"
        "%zu,%zu,%zu,%zu,%zu,%zu,%zu,"
        "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%ld,%.9e,%d,%d,%s",
        row.phase.c_str(),
        row.problem.c_str(),
        row.family.c_str(),
        row.n,
        row.N,
        row.param,
        row.method.c_str(),
        row.precond.c_str(),
        row.start_kind.c_str(),
        row.start_id.c_str(),
        row.distance.c_str(),
        row.rho_x,
        row.seed,
        row.success ? 1 : 0,
        row.nonlinear_iters,
        row.n_F,
        row.n_JVP,
        row.n_D2,
        row.n_lin,
        row.n_lin_failures,
        row.n_prec,
        row.t_setup,
        row.t_F,
        row.t_JVP,
        row.t_D2,
        row.t_prec,
        row.t_total,
        row.max_rss_kb,
        row.residual_rel,
        row.threads,
        row.repetition,
        row.commit.c_str());
    return buffer;
}

inline bool append_csv(const std::string& path, const CsvRow& row) {
    bool write_header = true;
    if (FILE* file = std::fopen(path.c_str(), "r")) {
        std::fseek(file, 0, SEEK_END);
        write_header = std::ftell(file) == 0;
        std::fclose(file);
    }

    FILE* file = std::fopen(path.c_str(), "a");
    if (file == nullptr) return false;
    if (write_header) {
        std::fprintf(file, "%s\n", csv_header());
    }
    std::fprintf(file, "%s\n", csv_line(row).c_str());
    std::fclose(file);
    return true;
}

} // namespace numa::instrument
