#include <numa/problems/brusselator1d.hpp>
#include <numa/problems/chandrasekhar_op.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace numa::op;

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

double norm2_local(
    const DVec& values)
{
    double sum = 0.0;

    for (double value : values) {
        sum += value * value;
    }

    return std::sqrt(sum);
}

}

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr
            << "Uso: materialize_simple_reference "
            << "B5|Bext n parametro salida.txt\n";
        return 2;
    }

    const std::string family = argv[1];

    const std::size_t n =
        static_cast<std::size_t>(
            std::stoull(argv[2]));

    const double parameter =
        std::stod(argv[3]);

    const std::string output = argv[4];

    DVec reference;

    if (family == "B5") {
        ChandrasekharOp problem(
            n,
            parameter);

        reference.assign(
            problem.size(),
            1.0);
    } else if (family == "Bext") {
        Brusselator1DStep problem(
            n,
            parameter);

        reference =
            problem.previous();

        if (Brusselator1DStep::A() != 1.0 ||
            Brusselator1DStep::B() != 3.0 ||
            Brusselator1DStep::diffusion() !=
                1.0 / 50.0) {

            std::cerr
                << "Parámetros Bext inesperados\n";
            return 3;
        }
    } else {
        std::cerr
            << "Familia desconocida\n";
        return 2;
    }

    if (!write_vector(
            output,
            reference)) {

        std::cerr
            << "Falló escritura\n";
        return 4;
    }

    std::cout
        << std::setprecision(17)
        << "family=" << family << '\n'
        << "n=" << n << '\n'
        << "parameter=" << parameter << '\n'
        << "N=" << reference.size() << '\n'
        << "norm="
        << norm2_local(reference) << '\n'
        << "SIMPLE_REFERENCE=PASS\n";

    return 0;
}
