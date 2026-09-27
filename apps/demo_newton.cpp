#include <numa/methods/newton.hpp>
#include <numa/problems/chemical_reactor.hpp>
#include <cstdio>

int main() {
    using namespace numa;
    using namespace numa::problems;

    ChemicalReactor f;
    Vec<2> x0 = {0.8, 0.4};
    auto r = newton<2>(f, x0, 1e-12, 20);

    std::printf("METODO DE NEWTON (2x2)\n");
    std::printf("x* = (%.17g, %.17g)\n", r.x[0], r.x[1]);
    std::printf("iteraciones = %zu, error = %.3e\n", r.iterations, r.er_abs);
    return 0;
}
