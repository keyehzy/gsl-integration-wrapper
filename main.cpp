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
     AdaptativeMesh(GSLFunction<F> f, double start, double end,
                    double epsabs = 1e-8, double epsrel = 1e-8, int maxiter = 10)
          : m_function(f), m_epsabs(epsabs), m_epsrel(epsrel),
            m_maxiter(maxiter) {
          Interval segment{start, end};
          ReusableNodes renodes {
               m_function(scaled_node(-1.0, segment)),
               m_function(scaled_node(0.0, segment)),
               m_function(scaled_node(1.0, segment))
          };
          mesh_interval(segment, renodes);
     }

     double point(int i) const {
          return m_points[i];
     }

     double value(int i) const {
          return m_values[i];
     }

     std::vector<double> const& points() const {
          return m_points;
     }

     std::vector<double> const& values() const {
          return m_values;
     }

     int size() const {
          return m_size;
     }

private:
     void mesh_interval(Interval segment, ReusableNodes renodes, int depth = 0);

     double scaled_node(double x, Interval segment) {
          return 0.5 * (x + 1.0) * (segment.end - segment.start) + segment.start;
     }

     double scaled_weight(double w, Interval segment) {
          return 0.5 * w * (segment.end - segment.start);
     }

     static constexpr std::array<double, 5> pts = { 
          -1.0, -0.5, 0.0, 0.5, 1.0
     };

     static constexpr std::array<double, 5> rule3 = {
          1.0/3.0, 0.0, 4.0/3.0, 0.0, 1.0/3.0
     };

     static constexpr std::array<double, 5> rule5 = {
          1.0/3.0, 4.0/3.0, 2.0/3.0, 4.0/3.0, 1.0/3.0
     };

     GSLFunction<F> m_function;
     double m_epsabs;
     double m_epsrel;
     int m_maxiter;
     int m_size = 0;
     std::vector<double> m_values;
     std::vector<double> m_points;
};

template <typename F>
void AdaptativeMesh<F>::mesh_interval(Interval segment, ReusableNodes renodes,
                                      int depth) {
     // std::printf("start: %lf end: %lf depth: %d\n", segment.start, segment.end, depth);
     
     if (depth > m_maxiter) {
          std::fprintf(stderr, "max interations, bailing out\n");
          return;
     }

     double mid_left = m_function(scaled_node(-0.5, segment));
     double mid_right = m_function(scaled_node(0.5, segment));

     std::array<double, 5> nodes{renodes.left, mid_left, renodes.mid, mid_right, renodes.right};

     double i0 = 0, i1 = 0;
     
     for (int i = 0; i < 5; i++) {
          i0 += scaled_weight(rule3[i], segment) * nodes[i];
          i1 += 0.5 * scaled_weight(rule5[i], segment) * nodes[i];
     }

     double epsabs_estimate = fabs(i1 - i0);
     double epsrel_estimate = fabs(1.0 - i0 / i1);

     if (epsabs_estimate < m_epsabs && epsrel_estimate < m_epsrel) {
          std::fprintf(stderr, "tolerance for error reached, bailing out\n");
          for (int i = 0; i < 5; i++) {
               m_points.push_back(scaled_node(pts[i], segment));
               m_values.push_back(nodes[i]);
               m_size += 1;
          }
     } else {
          mesh_interval(Interval{segment.start, segment.mid()},
                        ReusableNodes{renodes.left, mid_left, renodes.mid},
                        depth + 1);
          mesh_interval(Interval{segment.mid(), segment.end},
                        ReusableNodes{renodes.mid, mid_right, renodes.right},
                        depth + 1);
     }
}


double quantum_metric_1d(double m) {
  auto f = make_function_inv([&](double k) {
    return 0.25 / (k * k + m * m);
  });
  return cquadi(f);
}

double quantum_metric_2d(double m) {
  auto x = make_function_inv([&](double kx) {
    auto y = make_function_inv([&](double ky) {
      return 0.25 / (kx * kx + ky * ky + m * m) - 0.25 * kx * kx / pow(kx * kx + ky * ky + m * m, 2);
    });
    return cquadi(y);
  });
  return cquadi(x);
}

double quantum_metric_3d(double m) {
  auto x = make_function_inv([&](double kx) {
    auto y = make_function_inv([&](double ky) {
      auto z = make_function_inv([&](double kz) {
        return 0.25 / (kx * kx + ky * ky + kz * kz + m * m) -
          0.25 * kx * kx /
          pow(kx * kx + ky * ky + kz * kz + m * m, 2);
      });
      return cquadi(z);
    });
    return cquadi(y);
  });
  return cquadi(x);
}

std::vector<double> linspace(double start, double end, int n) {
  std::vector<double> x(n);
  double step = (end - start) / static_cast<double>(n);
  for (int i = 0; i < n; i++) {
    x[i] = start + static_cast<double>(i) * step;
  }
  return x;
}

int main(int argc, char** argv)
{
     double m = 0.5;


     auto f = ([&](double kx, double ky, double kz) {
          return 0.25 / (kx * kx + ky * ky + kz * kz + m * m) - 0.25 * kx * kx / pow(kx * kx + ky * ky + kz * kz + m * m, 2);
     });

     // auto f = ([&](double kx, double ky) {
     //      return 0.25 / (kx * kx + ky * ky + m * m) - 0.25 * kx * kx / pow(kx * kx + ky * ky + m * m, 2);
     // });

     // auto f = make_function_inv([&](double k) {
     //      return 0.25 / (k * k + m * m);
     // });

     for (double kx : linspace(-M_PI, M_PI, 100)) {
          for (double ky : linspace(-M_PI, M_PI, 100)) {
               std::printf("%.20lf %.20lf %.20lf\n", kx, ky, f(kx, ky, m));
          }
     }

     return 0;
}

// int main(int argc, char** argv)
// {
//      // Acquire mesh
//      auto g = make_function_inv([&](double m) {
//           return quantum_metric_3d(m);
//      });

//      AdaptativeMesh mesh(g, -0.5, 0.5, 1e-8, 1e-8, 25);

//      for (int i = 0; i < mesh.size(); i++) {
//           std::printf("%.20lf %.20lf\n", mesh.point(i), mesh.value(i));
//      }
// }
