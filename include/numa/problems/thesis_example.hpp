#pragma once

#include "../functions/function.hpp"
#include <cmath>
#include <limits>

namespace numa::problems {

struct ThesisExample final : IFunction<2> {
    Vec<2> eval(const Vec<2>& x) const override {
        const double X = x[0];
        const double Y = x[1];

        return {
            X * X + Y * Y - 2.0,
            std::exp(X - 1.0) + Y * Y * Y - 2.0
        };
    }

    Matrix<2, 2> jacobian(const Vec<2>& x) const override {
        const double X = x[0];
        const double Y = x[1];

        return {{
            {2.0 * X, 2.0 * Y},
            {std::exp(X - 1.0), 3.0 * Y * Y}
        }};
    }

    Matrix<2, 4> hessian(const Vec<2>& x) const override {
        const double X = x[0];
        const double Y = x[1];

        return {{
            {2.0, 0.0, 0.0, 2.0},
            {std::exp(X - 1.0), 0.0, 0.0, 6.0 * Y}
        }};
    }

    double lipschitz_k() const override {
        return std::numeric_limits<double>::infinity();
    }
};

} // namespace numa::problems
