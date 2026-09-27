#pragma once
#include "../core/matrix.hpp"
#include "../core/vector.hpp"

namespace numa {

template <std::size_t N>
void gauss_seidel(Matrix<N, N>& D, Vec<N>& d) {
    for (std::size_t k = 0; k + 1 < N; ++k) {
        for (std::size_t i = k + 1; i < N; ++i) {
            const double m = D[i][k] / D[k][k];
            D[i][k] = 0.0;
            for (std::size_t j = k + 1; j < N; ++j)
                D[i][j] -= m * D[k][j];
            d[i] -= m * d[k];
        }
    }
}

} // namespace numa
