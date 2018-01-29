#include "utils/math.hpp"
#include "stdlib.h"

namespace util {

float randf() {
  return (float)rand()/(float)RAND_MAX;
}

}
