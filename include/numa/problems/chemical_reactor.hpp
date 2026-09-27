#pragma once
#include "../functions/function.hpp"
#include <cmath>

namespace numa::problems {

struct ChemicalReactor final : IFunction<2> {
    Vec<2> eval(const Vec<2>& x) const override {
        const double X = x[0], Y = x[1];
        return {
            1.6*X*X + 3.6*X*Y - 7.8*X - 2.6*Y + 5.2,
            3.1*X*X + 0.9*Y*Y - 6.2*X + 6.2*Y
        };
    }
    Matrix<2,2> jacobian(const Vec<2>& x) const override {
        const double X = x[0], Y = x[1];
        return {{
            { 3.2*X + 3.6*Y - 7.8, 3.6*X - 2.6 },
            { 6.2*X - 6.2,         1.8*Y + 6.2 }
        }};
    }
    Matrix<2,4> hessian(const Vec<2>&) const override {
        return {{
            { 3.2, 3.6, 3.6, 0.0 },
            { 6.2, 0.0, 0.0, 1.8 }
        }};
    }
    double lipschitz_k() const override { return 6.8; }
};

} // namespace numa::problems
