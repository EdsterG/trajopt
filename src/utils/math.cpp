#include "utils/math.hpp"
#include "stdlib.h"
#include <cmath>

namespace util {

float randf() {
  return (float)rand()/(float)RAND_MAX;
}

// This implementation of mod returns sign of the divisor
double mod(double a, double n) {
    return fmod(fmod(a, n) + n, n);
}

double min_angle_diff(double x, double y) {
    return mod(x - y + M_PI, 2*M_PI) - M_PI;
}

}
