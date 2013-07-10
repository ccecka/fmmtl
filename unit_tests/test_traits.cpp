
#include "KernelTraits.hpp"

#include "UnitKernel.kern"

#include <iostream>

int main() {
  std::cout << ExpansionTraits<UnitExpansion>() << std::endl;
}
