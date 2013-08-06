#include "P2P_Compressed.hpp"

#include <thrust/device_malloc.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <cstdio>

template <typename value_type, typename Container>
inline value_type* gpu_copy(const Container& c) {
  typedef thrust::device_vector<typename Container::value_type> d_vector;
  d_vector* d = new d_vector(c.begin(), c.end());
  return reinterpret_cast<value_type*>(thrust::raw_pointer_cast(d->data()));
}

inline void gpuAssert(cudaError_t code, char* file, int line, bool abort=true)
{
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n",cudaGetErrorString(code),file,line);
    if (abort) exit(code);
  }
}
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }


template <unsigned BLOCKSIZE,
          typename Kernel,
          typename RandomAccessIterator1,  // pair<uint,uint>
          typename RandomAccessIterator2,  // Chained uint
          typename RandomAccessIterator3>  // pair<uint,uint>
__global__ void
blocked_p2p(const Kernel K,  // The Kernel to apply
            // BlockIdx -> pair<uint,uint> target range
            RandomAccessIterator1 target_range,
            // BlockIdx,BlockIdx+1 -> pair<uint,uint> into source_range
            RandomAccessIterator2 source_range_ptr,
            // Idx -> pair<uint,uint> source range
            RandomAccessIterator3 source_range,
            const typename Kernel::source_type* source,
            const typename Kernel::charge_type* charge,
            const typename Kernel::target_type* target,
            typename Kernel::result_type* result) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  // Get the target range this block is responsible for
  thrust::pair<unsigned,unsigned> t_range = target_range[blockIdx.x];
  unsigned t_first = t_range.first;
  unsigned t_last  = t_range.second;

  // Get the range of source ranges this block is responsible for
  RandomAccessIterator3 sr_first = source_range + source_range_ptr[blockIdx.x+0];
  RandomAccessIterator3 sr_last  = source_range + source_range_ptr[blockIdx.x+1];

  // Parallel for each target until the last
  for (t_first += threadIdx.x; t_first < t_last; t_first += BLOCKSIZE) {
    const target_type t = target[t_first];
    result_type r = result_type();

    // For each source range
    for (; sr_first < sr_last; ++sr_first) {
      thrust::pair<unsigned,unsigned> s_range = *sr_first;
      unsigned s_first = s_range.first;
      unsigned s_last  = s_range.second;
      // TODO: Load in shared memory and reuse?

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


struct Data {
  unsigned num_sources;
  unsigned num_targets;
  unsigned num_target_blocks;

  Data(unsigned s, unsigned t, unsigned b)
      : num_sources(s), num_targets(t), num_target_blocks(b) {
  }
};

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed()
    : data(0) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::P2P_Compressed(
    std::vector<std::pair<unsigned,unsigned> >& target_ranges,
    std::vector<unsigned>& source_range_ptrs,
    std::vector<std::pair<unsigned,unsigned> >& source_ranges,
    std::vector<typename Kernel::source_type>& sources,
    std::vector<typename Kernel::target_type>& targets)
    : data(new Data(sources.size(), targets.size(), target_ranges.size())),
      target_ranges(gpu_copy<unsigned>(target_ranges)),
      source_range_ptrs(gpu_copy<unsigned>(source_range_ptrs)),
      source_ranges(gpu_copy<unsigned>(source_ranges)),
      sources(gpu_copy<typename Kernel::source_type>(sources)),
      targets(gpu_copy<typename Kernel::target_type>(targets)) {
}

template <typename Kernel>
P2P_Compressed<Kernel>::~P2P_Compressed() {
  // TODO: delete
}

template <typename Kernel>
void P2P_Compressed<Kernel>::execute(
    const Kernel& K,
    const std::vector<typename Kernel::charge_type>& charges,
    std::vector<typename Kernel::result_type>& results) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  thrust::device_vector<charge_type> d_charges(charges);
  thrust::device_vector<result_type> d_results(results.size());

  const unsigned threads_per_block = 256;
  const unsigned num_blocks = reinterpret_cast<Data*>(data)->num_target_blocks;

  // Launch kernel <<<grid_size, block_size>>>
  blocked_p2p<threads_per_block><<<num_blocks,threads_per_block>>>(
      K,
      reinterpret_cast<thrust::pair<unsigned,unsigned>*>(target_ranges),
      source_range_ptrs,
      reinterpret_cast<thrust::pair<unsigned,unsigned>*>(source_ranges),
      sources,
      thrust::raw_pointer_cast(d_charges.data()),
      targets,
      thrust::raw_pointer_cast(d_results.data()));

  //gpuErrchk( cudaPeekAtLastError() );
  //gpuErrchk( cudaDeviceSynchronize() );

  // Copy results back and add
  thrust::host_vector<result_type> h_results = d_results;

  for (unsigned k = 0; k < h_results.size(); ++k)
    results[k] += h_results[k];
}


struct target_range_maker
    : public thrust::unary_function<unsigned,
                                    thrust::pair<unsigned,unsigned> > {
  __device__
  thrust::pair<unsigned,unsigned> operator()(unsigned target_block) const {
    unsigned start_block = target_block * 256;
    return thrust::make_pair(start_block, start_block + 256);
  }
};

template <typename Kernel>
void
P2P_Compressed<Kernel>::execute(const Kernel& K,
                                const std::vector<source_type>& s,
                                const std::vector<charge_type>& c,
                                const std::vector<target_type>& t,
                                std::vector<result_type>& r) {
  typedef Kernel kernel_type;
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  thrust::device_vector<source_type> d_sources(s);
  thrust::device_vector<charge_type> d_charges(c);
  thrust::device_vector<target_type> d_targets(t);
  thrust::device_vector<result_type> d_results(r);

  const unsigned threads_per_block = 256;
  const unsigned num_blocks = (t.size() + threads_per_block - 1) / threads_per_block;

  // Launch kernel <<<grid_size, block_size>>>
  blocked_p2p<threads_per_block><<<num_blocks, threads_per_block>>>(
      K,
      thrust::make_transform_iterator(thrust::make_counting_iterator(0),
                                      target_range_maker()),
      thrust::make_constant_iterator(0),
      thrust::make_constant_iterator(thrust::make_pair(0,s.size())),
      thrust::raw_pointer_cast(d_sources.data()),
      thrust::raw_pointer_cast(d_charges.data()),
      thrust::raw_pointer_cast(d_targets.data()),
      thrust::raw_pointer_cast(d_results.data()));

  // Copy results back and add
  thrust::host_vector<result_type> h_results = d_results;

  for (unsigned k = 0; k < h_results.size(); ++k)
    r[k] += h_results[k];
}

