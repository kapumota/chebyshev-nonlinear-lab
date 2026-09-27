#pragma once
#include "../core/matrix.hpp"
#include "../core/vector.hpp"

namespace numa {

template <std::size_t N>
Vec<N> tri_sup(const Matrix<N, N>& D, const Vec<N>& d) {
    Vec<N> x{};
    x[N - 1] = d[N - 1] / D[N - 1][N - 1];
    for (std::size_t kk = N - 1; kk-- > 0; ) {
        double s = 0.0;
        for (std::size_t j = kk + 1; j < N; ++j)
            s += D[kk][j] * x[j];
        x[kk] = (d[kk] - s) / D[kk][kk];
    }
    return x;
}

} // namespace numa
