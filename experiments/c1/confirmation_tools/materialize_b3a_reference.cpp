// Exporta exact() de la implementación existente, sin resolver un sistema.
#include <numa/problems/burgers2d.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Uso: materialize_b3a_reference n nu dt amp salida.txt\n";
        return 2;
    }
    numa::op::Burgers2D problem(std::stoull(argv[1]), std::stod(argv[2]),
                              std::stod(argv[3]), std::stod(argv[4]));
    std::ofstream out(argv[5]);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (double value : problem.exact()) out << value << '\n';
    return out ? 0 : 2;
}
