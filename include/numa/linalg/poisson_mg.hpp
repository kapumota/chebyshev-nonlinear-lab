#pragma once
// Multigrid geométrico de un ciclo V para el operador -Laplaciano_h en el cubo
// unitario con condición de Dirichlet nula, stencil de 7 puntos y n^3 nodos
// interiores. Se usa como precondicionador de Bratu3D, donde el término
// exponencial es de orden menor que el Laplaciano.
//
// Suavizador: Gauss-Seidel rojo-negro.
// Restricción: promedio ponderado de 27 puntos.
// Prolongación: interpolación trilineal.
// Se engrosa mientras n sea impar; n = 2^m - 1 produce la jerarquía completa.
// El objeto guarda áreas de trabajo y no es seguro entre hilos.

#include "../operators/operator.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numa::op {

class PoissonMG {
public:
    explicit PoissonMG(std::size_t n, int pre = 2, int post = 2)
        : pre_(pre), post_(post)
    {
        if (n == 0) {
            throw std::invalid_argument(
                "PoissonMG requiere un tamaño de malla positivo.");
        }
        if (pre < 0 || post < 0) {
            throw std::invalid_argument(
                "PoissonMG requiere números no negativos de suavizados.");
        }

        const std::size_t max_size = std::numeric_limits<std::size_t>::max();
        if (n > max_size / n || n * n > max_size / n) {
            throw std::overflow_error(
                "PoissonMG no puede representar n^3 en size_t.");
        }

        std::size_t m = n;
        double h2inv =
            (static_cast<double>(n) + 1.0) *
            (static_cast<double>(n) + 1.0);

        while (true) {
            Level level;
            level.n = m;
            level.h2inv = h2inv;

            const std::size_t N = m * m * m;
            level.x.assign(N, 0.0);
            level.b.assign(N, 0.0);
            level.r.assign(N, 0.0);
            levels_.push_back(std::move(level));

            if (m % 2 == 0 || m <= 1) break;
            m = (m - 1) / 2;
            h2inv /= 4.0;
        }
    }

    std::size_t levels() const { return levels_.size(); }

    // z = aproximación de A^{-1} r con un ciclo V.
    void operator()(const DVec& r, DVec& z) const {
        auto& finest = levels_[0];
        if (r.size() != finest.x.size()) {
            throw std::invalid_argument(
                "PoissonMG recibió un residual con dimensión incompatible.");
        }

        finest.b = r;
        std::fill(finest.x.begin(), finest.x.end(), 0.0);
        cycle(0);
        z = finest.x;
    }

private:
    struct Level {
        std::size_t n = 0;
        double h2inv = 0.0;
        DVec x;
        DVec b;
        DVec r;
    };

    static std::size_t index(
        std::size_t n, std::size_t i, std::size_t j, std::size_t k)
    {
        return (i * n + j) * n + k;
    }

    // Una barrida de Gauss-Seidel por colores.
    // order = 0 usa rojo y luego negro; order = 1 invierte el orden.
    static void smooth(Level& level, int order) {
        const std::size_t n = level.n;
        const double inv = 1.0 / level.h2inv;

        for (int pass = 0; pass < 2; ++pass) {
            const std::size_t color =
                static_cast<std::size_t>(order == 0 ? pass : 1 - pass);

            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    for (std::size_t k = (i + j + color) & 1;
                         k < n; k += 2) {
                        const std::size_t id = index(n, i, j, k);
                        double sum = level.b[id] * inv;

                        if (i > 0) sum += level.x[id - n * n];
                        if (i + 1 < n) sum += level.x[id + n * n];
                        if (j > 0) sum += level.x[id - n];
                        if (j + 1 < n) sum += level.x[id + n];
                        if (k > 0) sum += level.x[id - 1];
                        if (k + 1 < n) sum += level.x[id + 1];

                        level.x[id] = sum / 6.0;
                    }
                }
            }
        }
    }

    static void compute_residual(Level& level) {
        const std::size_t n = level.n;

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t k = 0; k < n; ++k) {
                    const std::size_t id = index(n, i, j, k);
                    double sum = 6.0 * level.x[id];

                    if (i > 0) sum -= level.x[id - n * n];
                    if (i + 1 < n) sum -= level.x[id + n * n];
                    if (j > 0) sum -= level.x[id - n];
                    if (j + 1 < n) sum -= level.x[id + n];
                    if (k > 0) sum -= level.x[id - 1];
                    if (k + 1 < n) sum -= level.x[id + 1];

                    level.r[id] =
                        level.b[id] - level.h2inv * sum;
                }
            }
        }
    }

    // El nodo grueso I corresponde al nodo fino 2I+1.
    static void restrict_to(const Level& fine, Level& coarse) {
        const std::size_t nf = fine.n;
        const std::size_t nc = coarse.n;
        static const double weights[3] = {0.25, 0.5, 0.25};

        for (std::size_t I = 0; I < nc; ++I) {
            for (std::size_t J = 0; J < nc; ++J) {
                for (std::size_t K = 0; K < nc; ++K) {
                    double sum = 0.0;

                    for (int a = 0; a < 3; ++a) {
                        const std::size_t fi =
                            2 * I + static_cast<std::size_t>(a);
                        if (fi >= nf) continue;

                        for (int b = 0; b < 3; ++b) {
                            const std::size_t fj =
                                2 * J + static_cast<std::size_t>(b);
                            if (fj >= nf) continue;

                            for (int c = 0; c < 3; ++c) {
                                const std::size_t fk =
                                    2 * K + static_cast<std::size_t>(c);
                                if (fk >= nf) continue;

                                sum +=
                                    weights[a] * weights[b] * weights[c] *
                                    fine.r[index(nf, fi, fj, fk)];
                            }
                        }
                    }

                    coarse.b[index(nc, I, J, K)] = sum;
                }
            }
        }
    }

    // Devuelve los términos de interpolación 1D para un nodo fino.
    static void interpolation_terms(
        std::size_t fine_index,
        std::size_t coarse_size,
        std::size_t coarse_indices[2],
        double weights[2],
        int& count)
    {
        count = 0;

        if (fine_index % 2 == 1) {
            coarse_indices[count] = (fine_index - 1) / 2;
            weights[count] = 1.0;
            ++count;
            return;
        }

        const std::size_t right = fine_index / 2;

        if (right >= 1) {
            coarse_indices[count] = right - 1;
            weights[count] = 0.5;
            ++count;
        }

        if (right < coarse_size) {
            coarse_indices[count] = right;
            weights[count] = 0.5;
            ++count;
        }
    }

    static void prolong_add(const Level& coarse, Level& fine) {
        const std::size_t nf = fine.n;
        const std::size_t nc = coarse.n;

        for (std::size_t i = 0; i < nf; ++i) {
            std::size_t ci[2];
            double wi[2];
            int ni = 0;
            interpolation_terms(i, nc, ci, wi, ni);

            for (std::size_t j = 0; j < nf; ++j) {
                std::size_t cj[2];
                double wj[2];
                int nj = 0;
                interpolation_terms(j, nc, cj, wj, nj);

                for (std::size_t k = 0; k < nf; ++k) {
                    std::size_t ck[2];
                    double wk[2];
                    int nk = 0;
                    interpolation_terms(k, nc, ck, wk, nk);

                    double sum = 0.0;
                    for (int a = 0; a < ni; ++a) {
                        for (int b = 0; b < nj; ++b) {
                            for (int c = 0; c < nk; ++c) {
                                sum +=
                                    wi[a] * wj[b] * wk[c] *
                                    coarse.x[index(
                                        nc, ci[a], cj[b], ck[c])];
                            }
                        }
                    }

                    fine.x[index(nf, i, j, k)] += sum;
                }
            }
        }
    }

    void cycle(std::size_t level_index) const {
        Level& level = levels_[level_index];
        const bool coarsest = level_index + 1 == levels_.size();

        if (coarsest) {
            if (level.n == 1) {
                level.x[0] =
                    level.b[0] / (6.0 * level.h2inv);
            } else {
                for (int sweep = 0; sweep < 30; ++sweep) {
                    smooth(level, sweep % 2);
                }
            }
            return;
        }

        for (int sweep = 0; sweep < pre_; ++sweep) {
            smooth(level, 0);
        }

        compute_residual(level);

        Level& coarse = levels_[level_index + 1];
        restrict_to(level, coarse);
        std::fill(coarse.x.begin(), coarse.x.end(), 0.0);

        cycle(level_index + 1);
        prolong_add(coarse, level);

        for (int sweep = 0; sweep < post_; ++sweep) {
            smooth(level, 1);
        }
    }

    int pre_;
    int post_;
    mutable std::vector<Level> levels_;
};

// Proveedor para problemas cuyo operador principal es el Laplaciano.
// No depende del punto x.
struct PoissonMGProvider {
    std::shared_ptr<PoissonMG> mg;

    explicit PoissonMGProvider(std::size_t n)
        : mg(std::make_shared<PoissonMG>(n)) {}

    struct Bound {
        const PoissonMG* mg = nullptr;

        void operator()(const DVec& r, DVec& z) const {
            (*mg)(r, z);
        }
    };

    Bound make(const DVec&) const {
        return Bound{mg.get()};
    }
};

} // namespace numa::op
