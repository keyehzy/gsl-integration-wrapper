#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <numbers>

#include "integ.h"

std::vector<double> linspace(double start, double end, int n) {
  std::vector<double> x(n);
  double step = (end - start) / static_cast<double>(n);
  for (int i = 0; i < n; i++) {
    x[i] = start + static_cast<double>(i) * step;
  }
  return x;
}

void quantum_metric_1d() {
  for (double m : linspace(-0.5, 0.5, 50)) {
    auto f = make_function([&](double k) { return 0.25 / (k * k + m * m); });
    std::cout << m << " " << cquad(f, -1000.0, 1000.0) << '\n';
  }
}

void quantum_metric_2d() {
  for (double m : linspace(-0.5, 0.5, 50)) {
    auto x = make_function([&](double kx) {
      auto y = make_function([&](double ky) {
        return 0.25 / ((kx * kx + ky * ky) + m * m) -
               0.25 * kx * kx / gsl_pow_2((kx * kx + ky * ky) + m * m);
      });
      return cquad(y, -1000.0, 1000.0);
    });
    std::cout << m << " " << cquad(x, -1000.0, 1000.0) << '\n';
  }
}

void quantum_metric_3d() {
  for (double m : linspace(-0.5, 0.5, 50)) {
    auto x = make_function([&](double kx) {
      auto y = make_function([&](double ky) {
        auto z = make_function([&](double kz) {
          return 0.25 / ((kx * kx + ky * ky + kz * kz) + m * m) -
                 0.25 * kx * kx /
                     gsl_pow_2((kx * kx + ky * ky + kz * kz) + m * m);
        });
        return cquad(z, -1000.0, 1000.0);
      });
      return cquad(y, -1000.0, 1000.0);
    });
    std::cout << m << " " << cquad(x, -1000.0, 1000.0) << '\n';
  }
}

int main(int argc, char **argv) {
  quantum_metric_3d();
  return 0;
}
