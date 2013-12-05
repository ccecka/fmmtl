#include "fmmtl/config.hpp"

#include <iostream>

int main(void)
{
  std::cout << "Using Boost v"
            << BOOST_VERSION / 100000     << "."  // major version
            << BOOST_VERSION / 100 % 1000 << "."  // minior version
            << BOOST_VERSION % 100                // patch level
            << std::endl;

#if defined(THRUST_MAJOR_VERSION) & defined(THRUST_MINOR_VERSION)
  int major = THRUST_MAJOR_VERSION;
  int minor = THRUST_MINOR_VERSION;
  std::cout << "Using Thrust v" << major << "." << minor << std::endl;
#else
  std::cout << "Thrust is not being used." << std::endl;
#endif
  return 0;
}
