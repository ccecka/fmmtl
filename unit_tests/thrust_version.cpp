#include "fmmtl/config.hpp"

#if defined(FMMTL_USE_THRUST)
#include <thrust/version.h>
#endif

#include <iostream>

int main(void)
{
#if defined(FMMTL_USE_THRUST)
  int major = THRUST_MAJOR_VERSION;
  int minor = THRUST_MINOR_VERSION;
  std::cout << "Thrust v" << major << "." << minor << std::endl;
#else
  std::cout << "Thrust is not being used." << std::endl;
#endif

  return 0;
}
