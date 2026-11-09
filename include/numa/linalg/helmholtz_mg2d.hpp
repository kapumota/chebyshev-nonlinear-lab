#pragma once
// Multigrid geométrico 2D de ciclo V para alpha I - nu Laplaciano_h en el cuadrado unitario,
// con Dirichlet nulo, stencil de 5 puntos y n^2 nodos interiores. Suavizador Gauss-Seidel
// rojo-negro, restricción de 9 puntos e interpolación bilineal. Se engrosa mientras n sea
// impar, con jerarquía completa si n = 2^m - 1. No es seguro entre hilos.

#include "../operators/operator.hpp"
#include <memory>
#include <vector>

namespace numa::op {

class HelmholtzMG2D {
public:
    HelmholtzMG2D(std::size_t n, double alpha, double nu, int pre = 2, int post = 2)
        : alpha_(alpha), pre_(pre), post_(post)
    {
        std::size_t m = n;
        double s = nu * (static_cast<double>(n) + 1.0) * (static_cast<double>(n) + 1.0);
        while (true) {
            Level L;
            L.n = m;
            L.s = s;
            L.x.assign(m * m, 0.0);
            L.b.assign(m * m, 0.0);
            L.r.assign(m * m, 0.0);
            levels_.push_back(std::move(L));
            if (m % 2 == 0 || m <= 1) break;
            m = (m - 1) / 2;
            s /= 4.0;
        }
    }

    std::size_t levels() const { return levels_.size(); }

    void operator()(const DVec& r, DVec& z) const {
        auto& L0 = levels_[0];
        L0.b = r;
        std::fill(L0.x.begin(), L0.x.end(), 0.0);
        cycle(0);
        z = L0.x;
    }

private:
    struct Level {
        std::size_t n = 0;
        double s = 0.0;   // nu / h^2
        DVec x, b, r;
    };

    void smooth(Level& L, int order) const {
        const std::size_t n = L.n;
        const double den = alpha_ + 4.0 * L.s;
        for (int pass = 0; pass < 2; ++pass) {
            const std::size_t color = (order == 0) ? pass : 1 - pass;
#pragma omp parallel for schedule(static) if (n * n > 4096)
            for (long ii = 0; ii < static_cast<long>(n); ++ii) {
                const std::size_t i = static_cast<std::size_t>(ii);
                for (std::size_t j = (i + color) & 1; j < n; j += 2) {
                    const std::size_t id = i * n + j;
                    double nb = 0.0;
                    if (i > 0)     nb += L.x[id - n];
                    if (i + 1 < n) nb += L.x[id + n];
                    if (j > 0)     nb += L.x[id - 1];
                    if (j + 1 < n) nb += L.x[id + 1];
                    L.x[id] = (L.b[id] + L.s * nb) / den;
                }
            }
        }
    }

    void residual(Level& L) const {
        const std::size_t n = L.n;
#pragma omp parallel for schedule(static) if (n * n > 4096)
        for (long ii = 0; ii < static_cast<long>(n); ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            for (std::size_t j = 0; j < n; ++j) {
                const std::size_t id = i * n + j;
                double nb = 0.0;
                if (i > 0)     nb += L.x[id - n];
                if (i + 1 < n) nb += L.x[id + n];
                if (j > 0)     nb += L.x[id - 1];
                if (j + 1 < n) nb += L.x[id + 1];
                L.r[id] = L.b[id] - (alpha_ * L.x[id] + L.s * (4.0 * L.x[id] - nb));
            }
        }
    }

    static void restrict_to(const Level& F, Level& C) {
        const std::size_t nf = F.n, nc = C.n;
        static const double w[3] = {0.25, 0.5, 0.25};
#pragma omp parallel for schedule(static) if (nc * nc > 4096)
        for (long II = 0; II < static_cast<long>(nc); ++II)
            for (std::size_t J = 0; J < nc; ++J) {
                const std::size_t I = static_cast<std::size_t>(II);
                double s = 0.0;
                for (int a = 0; a < 3; ++a) {
                    const std::size_t fi = 2 * I + a;
                    if (fi >= nf) continue;
                    for (int b = 0; b < 3; ++b) {
                        const std::size_t fj = 2 * J + b;
                        if (fj >= nf) continue;
                        s += w[a] * w[b] * F.r[fi * nf + fj];
                    }
                }
                C.b[I * nc + J] = s;
            }
    }

    static void terms(std::size_t f, std::size_t nc, std::size_t c[2], double w[2], int& cnt) {
        cnt = 0;
        if (f % 2 == 1) { c[cnt] = (f - 1) / 2; w[cnt] = 1.0; ++cnt; }
        else {
            const std::size_t r = f / 2;
            if (r >= 1) { c[cnt] = r - 1; w[cnt] = 0.5; ++cnt; }
            if (r < nc) { c[cnt] = r;     w[cnt] = 0.5; ++cnt; }
        }
    }

    static void prolong_add(const Level& C, Level& F) {
        const std::size_t nf = F.n, nc = C.n;
#pragma omp parallel for schedule(static) if (nf * nf > 4096)
        for (long ii = 0; ii < static_cast<long>(nf); ++ii) {
            const std::size_t i = static_cast<std::size_t>(ii);
            std::size_t ci[2]; double wi[2]; int ni;
            terms(i, nc, ci, wi, ni);
            for (std::size_t j = 0; j < nf; ++j) {
                std::size_t cj[2]; double wj[2]; int nj;
                terms(j, nc, cj, wj, nj);
                double s = 0.0;
                for (int a = 0; a < ni; ++a)
                    for (int b = 0; b < nj; ++b) s += wi[a] * wj[b] * C.x[ci[a] * nc + cj[b]];
                F.x[i * nf + j] += s;
            }
        }
    }

    void cycle(std::size_t l) const {
        Level& L = levels_[l];
        if (l + 1 == levels_.size()) {
            if (L.n == 1) L.x[0] = L.b[0] / (alpha_ + 4.0 * L.s);
            else for (int s = 0; s < 30; ++s) smooth(L, s % 2);
            return;
        }
        for (int s = 0; s < pre_; ++s) smooth(L, 0);
        residual(L);
        Level& C = levels_[l + 1];
        restrict_to(L, C);
        std::fill(C.x.begin(), C.x.end(), 0.0);
        cycle(l + 1);
        prolong_add(C, L);
        for (int s = 0; s < post_; ++s) smooth(L, 1);
    }

    double alpha_;
    int pre_, post_;
    mutable std::vector<Level> levels_;
};

} // namespace numa::op
