#pragma once

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <vector>
#include <cassert>

enum class GSLFunctionDomain {
     Euclidean,
     Stereographic,
};

template <typename F>
class GSLFunction : public gsl_function {
public:
     GSLFunction(const F &f, GSLFunctionDomain domain)
          :  gsl_function{.function = invoke, .params = this },
             m_function(f),
             m_domain(domain) {}

     operator gsl_function*() {
          return this;
     }

     double operator()(double x) const {
          return m_function(x);
     }

     GSLFunctionDomain domain() const { return m_domain; }

private:
     static double invoke(double x, void *params) {
          GSLFunction* function = static_cast<GSLFunction*>(params);
          switch(function->domain()) {
          case GSLFunctionDomain::Euclidean:
               return function->m_function(x);
          case GSLFunctionDomain::Stereographic:
               return function->m_function(1/x - 1);
          }
     }

     const F m_function;
     GSLFunctionDomain m_domain;
};

template <typename F>
GSLFunction<F> make_function(const F& func) {
     return GSLFunction<F>(func, GSLFunctionDomain::Euclidean);
}

template <typename F>
GSLFunction<F> make_function_inv(const F& func) {
     return GSLFunction<F>(func, GSLFunctionDomain::Stereographic);
}

class GSLIntegrationWorkspace {
public:
     GSLIntegrationWorkspace(std::size_t size = 1000)
          : m_size(size),
            m_workspace(gsl_integration_workspace_alloc(size)) {};

     ~GSLIntegrationWorkspace() {
          gsl_integration_workspace_free(m_workspace);
     }

     operator gsl_integration_workspace*() const {
          return m_workspace;
     }

     std::size_t size() const {
          return m_size;
     }

private:
     std::size_t m_size;
     gsl_integration_workspace *m_workspace;
};

enum class QAG_Order {
     K15 = 1,
     K21 = 2,
     K31 = 3,
     K41 = 4,
     K51 = 5,
     K61 = 6,
};

template <typename F>
double qag(GSLFunction<F> f, double a, double b,
           QAG_Order order = QAG_Order::K21,
           double epsabs = 1e-8, double epsrel = 1e-8) {
     double result, abserr;
     GSLIntegrationWorkspace w;
     gsl_integration_qag(f, a, b, epsabs = 1e-8, epsrel = 1e-8, w.size(),
                         static_cast<int>(order), w, &result,
                         &abserr);
     return result;
}

template <typename F>
double qags(GSLFunction<F> f, double a, double b, double epsabs = 1e-8,
            double epsrel = 1e-8) {
     double result, abserr;
     GSLIntegrationWorkspace w;
     gsl_integration_qags(f, a, b, epsabs = 1e-8, epsrel = 1e-8, w.size(), w,
                          &result, &abserr);
     return result;
}

template <typename F>
double qagp(GSLFunction<F> f, std::vector<double> pts, double epsabs = 1e-8,
            double epsrel = 1e-8) {
     double result, abserr;
     GSLIntegrationWorkspace w;
     gsl_integration_qagp(f, pts.data(), pts.size(), epsabs = 1e-8, epsrel = 1e-8, w.size(),
                          w, &result, &abserr);
     return result;
}

template <typename F>
double qagi(GSLFunction<F> f, double epsabs = 1e-8, double epsrel = 1e-8) {
     double result, abserr;
     GSLIntegrationWorkspace w;
     gsl_integration_qagi(f, epsabs = 1e-8, epsrel = 1e-8, w.size(), w, &result, &abserr);
     return result;
}

class GSLIntegrationQuadWorkspace {
public:
     GSLIntegrationQuadWorkspace(std::size_t size = 1000)
          : m_size(size),
            m_workspace(gsl_integration_cquad_workspace_alloc(size)) {};

     ~GSLIntegrationQuadWorkspace() {
          gsl_integration_cquad_workspace_free(m_workspace);
     }

     operator gsl_integration_cquad_workspace*() const {
          return m_workspace;
     }

     std::size_t size() const {
          return m_size;
     }

private:
     std::size_t m_size;
     gsl_integration_cquad_workspace *m_workspace;
};

template <typename F>
double cquad(GSLFunction<F> f, double a, double b,
           double epsabs = 1e-8, double epsrel = 1e-8) {
     double result, abserr;
     std::size_t nevals;
     GSLIntegrationQuadWorkspace w;
     gsl_integration_cquad(f, a, b, epsabs = 1e-8, epsrel = 1e-8,
                         w, &result, &abserr, &nevals);
     return result;
}

template <typename F>
double cquadi(GSLFunction<F> f, double epsabs = 1e-8, double epsrel = 1e-8) {
     assert(f.domain() == GSLFunctionDomain::Stereographic);
     return cquad(f, 0, 1, epsabs, epsrel);
}
