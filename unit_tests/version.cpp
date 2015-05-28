#include "fmmtl/config.hpp"

#include <iostream>

int main(void)
{
#if defined(THRUST_MAJOR_VERSION) & defined(THRUST_MINOR_VERSION)
  int major = THRUST_MAJOR_VERSION;
  int minor = THRUST_MINOR_VERSION;
  std::cout << "Using Thrust v" << major << "." << minor << std::endl;
#else
  std::cout << "Thrust is not being used." << std::endl;
#endif
  return 0;
}
