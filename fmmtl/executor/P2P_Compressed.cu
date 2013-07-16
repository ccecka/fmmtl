#include "P2P_Compressed.hpp"

#include <thrust/device_malloc.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include <cstdio>

template <typename value_type, typename Container>
inline value_type* gpu_copy(const Container& c) {
  typedef thrust::device_vector<typename Container::value_type> d_vector;
  d_vector* d = new d_vector(c.begin(), c.end());
  return reinterpret_cast<value_type*>(thrust::raw_pointer_cast(d->data()));
}


template <typename Kernel>
__global__ void dumb_p2p(const Kernel K,
                         const unsigned* target_ranges,  // pair<unit,uint>
                         const unsigned* target_ptrs,
                         const unsigned* source_ranges,  // pair<uint,uint>
                         const typename Kernel::source_type* source,
                         const typename Kernel::charge_type* charge,
                         const typename Kernel::target_type* target,
                         typename Kernel::result_type* result) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  unsigned t_first = target_ranges[2*blockIdx.x+0] + threadIdx.x;
  unsigned t_last  = target_ranges[2*blockIdx.x+1];

  unsigned sr_first = target_ptrs[blockIdx.x+0];
  unsigned sr_last  = target_ptrs[blockIdx.x+1];

  // For each target until the last
  for (; t_first < t_last; t_first += blockDim.x) {
    const target_type t = target[t_first];
    result_type r = result_type();

    // For each source range
    for (; sr_first != sr_last; ++sr_first) {
      unsigned s_first = source_ranges[2*sr_first+0];
      unsigned s_last  = source_ranges[2*sr_first+1];
      // TODO: Load in shared memory and reuse

      // For each source in the source range
      for (; s_first < s_last; ++s_first) {
        const source_type s = source[s_first];
        const charge_type c = charge[s_first];
        r += K(t,s) * c;
      }
    }

    // Assign the result
    result[t_first] = r;
  }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

struct Data {
  unsigned num_sources;
  unsigned num_targets;
  unsigned num_blocks;

  Data(unsigned s, unsigned t, unsigned b)
      : num_sources(s), num_targets(t), num_blocks(b) {
  }
};

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(const Kernel& K)
    : K_(K) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(
    const Kernel& K,
    std::vector<std::pair<unsigned,unsigned> >& target_ranges,
    std::vector<unsigned>& target_ptrs,
    std::vector<std::pair<unsigned,unsigned> >& source_ranges,
    std::vector<typename Kernel::source_type>& sources,
    std::vector<typename Kernel::target_type>& targets)
    : K_(K),
      data(new Data(sources.size(), targets.size(), target_ranges.size())),
      d_target_ranges(gpu_copy<unsigned>(target_ranges)),
      d_target_ptrs(gpu_copy<unsigned>(target_ptrs)),
      d_source_ranges(gpu_copy<unsigned>(source_ranges)),
      d_sources(gpu_copy<typename Kernel::source_type>(sources)),
      d_targets(gpu_copy<typename Kernel::target_type>(targets)) {
}

template <typename Kernel>
void P2P_Compressed<Kernel>::execute(
    const typename Kernel::charge_type* charges,
    typename Kernel::result_type* results) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  std::cout << "Launching GPU P2P Kernel\n";

  unsigned num_sources = reinterpret_cast<Data*>(data)->num_sources;
  unsigned num_targets = reinterpret_cast<Data*>(data)->num_targets;

  thrust::device_vector<charge_type> d_charges(charges, charges + num_sources);
  thrust::device_vector<result_type> d_results(num_targets);

  std::cout << "Host to Device\n";

  unsigned num_blocks = reinterpret_cast<Data*>(data)->num_blocks;

  // Launch kernel <<<grid_size, block_size>>>
  dumb_p2p<<<num_blocks,256>>>(K_,
                               d_target_ranges,
                               d_target_ptrs,
                               d_source_ranges,
                               d_sources,
                               thrust::raw_pointer_cast(d_charges.data()),
                               d_targets,
                               thrust::raw_pointer_cast(d_results.data()));

  gpuErrchk( cudaPeekAtLastError() );
  gpuErrchk( cudaDeviceSynchronize() );

  std::cout << "GPU Kernel Done\n";

  // Copy results back and add
  thrust::host_vector<typename Kernel::result_type> h_results = d_results;

  std::cout << "Device to Host\n";

  for (unsigned k = 0; k < h_results.size(); ++k) {
    results[k] += h_results[k];
  }
}
