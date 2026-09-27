#include <numa/methods/kantorovich.hpp>
#include <numa/problems/chandrasekhar.hpp>
#include <cstdio>

int main() {
    using namespace numa;
    using namespace numa::problems;

    Chandrasekhar<8> f;
    Vec<8> x0{};
    for (auto& v : x0) v = 1.0;

    const double L = 0.530394219033;
    const double a = 0.265197109517;

    auto r = kantorovich<8>(f, x0, L, a, 1e-12, 20);

    std::printf("NEWTON-KANTOROVICH (8x8)\n");
    for (std::size_t i = 0; i < 8; ++i)
        std::printf("x[%zu] = %.17g\n", i, r.x[i]);
    std::printf("iteraciones = %zu, er_abs = %.3e, er_teo = %.3e\n",
                r.iterations, r.er_abs, r.er_teorico);
    return 0;
}
