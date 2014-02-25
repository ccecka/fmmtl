#pragma once

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
    const std::vector<source_type>&,
    const std::vector<target_type>&) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::~P2P_Compressed() {
  // TODO: delete
}

template <typename Kernel>
void P2P_Compressed<Kernel>::execute(
    const Kernel&,
    const std::vector<charge_type>&,
    std::vector<result_type>&) {
  std::cerr << "ERROR: Calling unimplemented P2P_compressed CPU" << std::endl;
}

template <typename Kernel>
void P2P_Compressed<Kernel>::execute(
    const Kernel&,
    const std::vector<source_type>&,
    const std::vector<charge_type>&,
    const std::vector<target_type>&,
    std::vector<result_type>&) {
  std::cerr << "ERROR: Calling unimplemented P2P_compressed CPU" << std::endl;
}
