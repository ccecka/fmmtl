#include "P2P_Compressed.hpp"

#include <iostream>

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed() {
}

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(
    std::vector<std::pair<unsigned,unsigned> >&,
    std::vector<unsigned>&,
    std::vector<std::pair<unsigned,unsigned> >&,
    const std::vector<typename Kernel::source_type>&,
    const std::vector<typename Kernel::target_type>&) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::~P2P_Compressed() {
  // TODO: delete
}

template <typename Kernel>
void P2P_Compressed<Kernel>::execute(
    const Kernel&,
    const std::vector<typename Kernel::charge_type>&,
    std::vector<typename Kernel::result_type>&) {
}
