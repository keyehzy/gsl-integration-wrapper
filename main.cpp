#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <utility>

#include "integ.h"

int main(int argc, char** argv) {
     auto f = make_function([&] (double x) {
          return exp(-x*x); 
     });

     double result = qag(f, 0, 1);
     std::cout << result << std::endl;
     return 0;
}
