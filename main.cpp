#include <cstdio>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <array>
#include <numeric>
#include <algorithm>

#include "integ.h"

template <typename F>
class AdaptativeMesh {
private:
     struct Interval {
          double start, end;
          double mid () const { return 0.5 * (start + end); }
          double length() const { return end - start; }
     };

     struct ReusableNodes {
          double left, mid, right;
     };

public:
     AdaptativeMesh(GSLFunction<F> f) : m_function(f) {}

     std::vector<double> const& points() const {
          return m_points;
     }

     std::vector<double> mesh(double start, double end) {
          Interval segment{start, end};
          ReusableNodes renodes {
               m_function(scaled_node(-1.0, segment)),
               m_function(scaled_node(0.0, segment)),
               m_function(scaled_node(1.0, segment))
          };
          return mesh_interval(segment, renodes);          
     }

private:
     std::vector<double> mesh_interval(Interval segment, ReusableNodes renodes,
                                       double epsabs = 1e-8,
                                       double epsrel = 1e-8);

     double scaled_node(double x, Interval segment) {
          return 0.5 * (x + 1.0) * (segment.end - segment.start) + segment.start;
     }

     double scaled_weight(double w, Interval segment) {
          return 0.5 * w * (segment.end - segment.start);
     }

     static constexpr std::array<double, 5> rule3 = {
          1./3., 0., 4./3., 0., 1./3.
     };

     static constexpr std::array<double, 5> rule5 = {
          1./3., 4./3., 2./3., 4./3., 1./3.
     };

     GSLFunction<F> m_function;
     std::vector<double> m_points;
};

template <typename F>
std::vector<double> AdaptativeMesh<F>::mesh_interval(Interval segment,
                                                     ReusableNodes renodes,
                                                     double epsabs,
                                                     double epsrel) {
     double mid_left = m_function(scaled_node(-0.5, segment));
     double mid_right = m_function(scaled_node(0.5, segment));

     std::vector<double> nodes{renodes.left, mid_left, renodes.mid, mid_right, renodes.right};

     double i0 = std::inner_product(nodes.begin(), nodes.end(), rule3.begin(), 0);
     double i1 = 0.5 * std::inner_product(nodes.begin(), nodes.end(), rule5.begin(), 0);

     double epsabs_estimate = fabs(i1 - i0);
     double epsrel_estimate = fabs(1.0 - i0 / i1);

     if (epsabs_estimate < epsabs && epsrel_estimate < epsrel) {
          std::transform(nodes.begin(), nodes.end(), nodes.begin(), [&] (double n) {
               return scaled_node(n, segment);
          });
          return nodes;
     }

     double segment_mid = 0.5 * (segment.start + segment.end);
     std::vector<double> left_mesh = mesh_interval(Interval{segment.start, segment_mid}, ReusableNodes{ renodes.left, mid_left, renodes.mid });
     std::vector<double> right_mesh = mesh_interval(Interval{segment_mid, segment.end}, ReusableNodes{ renodes.mid, mid_right, renodes.right });
     std::vector<double> merge_mesh;
     std::merge(left_mesh.begin(), left_mesh.end(), right_mesh.begin(), right_mesh.end(), std::back_inserter(merge_mesh));
     return merge_mesh;
}

int main(int argc, char** argv)
{
     auto f = make_function([&](double x) {
          return x;
     });
     AdaptativeMesh adaptative(f);
     std::vector<double> mesh = adaptative.mesh(0.0, 1.0);
     for (double x : mesh) {
          std::cout << x << std::endl;
     }
}

// std::vector<double> linspace(double start, double end, int n) {
//   std::vector<double> x(n);
//   double step = (end - start) / static_cast<double>(n);
//   for (int i = 0; i < n; i++) {
//     x[i] = start + static_cast<double>(i) * step;
//   }
//   return x;
// }

// double quantum_metric_1d(double m) {
//   auto f = make_function_inv([&](double k) {
//     return 0.25 / (k * k + m * m);
//   });
//   return cquadi(f);
// }

// double quantum_metric_2d(double m) {
//   auto x = make_function_inv([&](double kx) {
//     auto y = make_function_inv([&](double ky) {
//       return 0.25 / (kx * kx + ky * ky + m * m) - 0.25 * kx * kx / pow(kx * kx + ky * ky + m * m, 2);
//     });
//     return cquadi(y);
//   });
//   return cquadi(x);
// }

// double quantum_metric_3d(double m) {
//   auto x = make_function_inv([&](double kx) {
//     auto y = make_function_inv([&](double ky) {
//       auto z = make_function_inv([&](double kz) {
//         return 0.25 / (kx * kx + ky * ky + kz * kz + m * m) -
//           0.25 * kx * kx /
//           pow(kx * kx + ky * ky + kz * kz + m * m, 2);
//       });
//       return cquadi(z);
//     });
//     return cquadi(y);
//   });
//   return cquadi(x);
// }

// int main(int argc, char **argv) {
//   for (double m : linspace(-0.5, 0.5 + 1e-4, 100)) {
//     std::printf("%lf %lf\n", m, quantum_metric_2d(m));
//   }
//   return 0;
// }
