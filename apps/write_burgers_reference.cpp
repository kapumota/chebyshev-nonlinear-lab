#include <numa/problems/burgers2d.hpp>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>

using namespace numa::op;

namespace {

bool parse_size(const char* text, std::size_t& value) {
    if (text == nullptr || *text == '\0') return false;
    errno = 0;
    char* end = nullptr;
    const unsigned long long raw = std::strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || raw == 0 ||
        raw > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
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
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(raw)) return false;
    value = raw;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr, "Uso: %s n nu dt amp salida.txt\n", argv[0]);
        return 2;
    }

    std::size_t n = 0;
    double nu = 0.0, dt = 0.0, amp = 0.0;
    if (!parse_size(argv[1], n) || !parse_double(argv[2], nu) ||
        !parse_double(argv[3], dt) || !parse_double(argv[4], amp) ||
        !(nu > 0.0) || !(dt > 0.0)) {
        std::fprintf(stderr, "Parámetros inválidos.\n");
        return 2;
    }

    Burgers2D problem(n, nu, dt, amp);
    std::ofstream output(argv[5]);
    if (!output) {
        std::fprintf(stderr, "No se pudo abrir %s.\n", argv[5]);
        return 3;
    }
    output.precision(17);
    for (double value : problem.exact()) output << value << '\n';
    return output ? 0 : 3;
}
