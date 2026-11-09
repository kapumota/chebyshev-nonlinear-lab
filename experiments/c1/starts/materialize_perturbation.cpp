#define main embedded_run_campaign_main
#include "../../../apps/run_campaign.cpp"
#undef main

#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr
            << "Uso: materialize_perturbation "
            << "reference.txt rho seed output.txt\n";
        return 2;
    }

    const std::string reference_path = argv[1];
    const double rho = std::stod(argv[2]);
    const unsigned seed =
        static_cast<unsigned>(
            std::stoul(argv[3])
        );
    const std::string output_path = argv[4];

    if (!(rho > 0.0)) {
        std::cerr << "rho debe ser positivo\n";
        return 2;
    }

    DVec reference;

    {
        std::ifstream input(reference_path);

        if (!input) {
            std::cerr
                << "No se pudo abrir "
                << reference_path
                << '\n';
            return 2;
        }

        double value = 0.0;

        while (input >> value) {
            if (!std::isfinite(value)) {
                std::cerr
                    << "Valor no finito\n";
                return 2;
            }

            reference.push_back(value);
        }

        if (!input.eof() || reference.empty()) {
            std::cerr
                << "Vector de referencia inválido\n";
            return 2;
        }
    }

    const DVec perturbed =
        perturb(reference, rho, seed);

    if (perturbed.size() != reference.size()) {
        std::cerr
            << "Falló la construcción "
            << "de la perturbación\n";
        return 2;
    }

    {
        std::ofstream output(output_path);

        if (!output) {
            std::cerr
                << "No se pudo crear "
                << output_path
                << '\n';
            return 2;
        }

        output
            << std::setprecision(
                std::numeric_limits<double>::max_digits10
            );

        for (double value : perturbed) {
            output << value << '\n';
        }

        if (!output) {
            std::cerr
                << "Falló escritura\n";
            return 2;
        }
    }

    DVec roundtrip;

    {
        std::ifstream input(output_path);

        double value = 0.0;

        while (input >> value) {
            roundtrip.push_back(value);
        }
    }

    if (roundtrip.size() != perturbed.size()) {
        std::cerr
            << "Falló round-trip\n";
        return 2;
    }

    for (std::size_t i = 0;
         i < perturbed.size();
         ++i) {

        const auto a =
            std::bit_cast<std::uint64_t>(
                perturbed[i]
            );

        const auto b =
            std::bit_cast<std::uint64_t>(
                roundtrip[i]
            );

        if (a != b) {
            std::cerr
                << "Round-trip no exacto en "
                << i
                << '\n';
            return 2;
        }
    }

    DVec delta(perturbed.size());

    for (std::size_t i = 0;
         i < delta.size();
         ++i) {
        delta[i] =
            perturbed[i] - reference[i];
    }

    const double relative_distance =
        norm2(delta) / norm2(reference);

    std::cout
        << std::setprecision(17)
        << "N=" << perturbed.size()
        << '\n'
        << "rho=" << rho
        << '\n'
        << "seed=" << seed
        << '\n'
        << "relative_distance="
        << relative_distance
        << '\n'
        << "BITWISE_ROUNDTRIP=PASS\n";

    return 0;
}
