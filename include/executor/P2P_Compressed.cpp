#include "P2P_Compressed.hpp"

#include <iostream>

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(const Kernel& K)
    : K_(K) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(
    const Kernel& K,
    std::vector<std::pair<unsigned,unsigned> >&,
    std::vector<unsigned>&,
    std::vector<std::pair<unsigned,unsigned> >&,
    std::vector<typename Kernel::source_type>&,
    std::vector<typename Kernel::target_type>&)
    : K_(K) {
}

template <typename Kernel>
void
P2P_Compressed<Kernel>::execute(const typename Kernel::charge_type*,
                                typename Kernel::result_type*) {

  std::cout << "Launching CPU P2P Kernel!\n";
}
