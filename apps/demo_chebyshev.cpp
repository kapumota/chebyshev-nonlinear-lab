#include <numa/methods/chebyshev.hpp>
#include <numa/problems/chandrasekhar.hpp>
#include <cstdio>

int main() {
    using namespace numa;
    using namespace numa::problems;

    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;

    auto r = chebyshev<8>(f, x0, 1e-12, 20);

    std::printf("METODO ITERATIVO DE CHEBYSHEV (8x8)\n");
    for (std::size_t i = 0; i < 8; ++i)
        std::printf("x[%zu] = %.17g\n", i, r.x[i]);
    std::printf("iteraciones = %zu, error = %.3e\n", r.iterations, r.er_abs);
    return 0;
}
