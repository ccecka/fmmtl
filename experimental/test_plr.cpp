
#include "PLR.hpp"





int main(int argc, char** argv) {
  unsigned N = 1 << 14;     // rows
  unsigned K = 1;           // cols of the rhs
  //unsigned leaf_size = 32;  // maximum size of the tree leaves

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    }
    if (strcmp(argv[i],"-K") == 0) {
      K = atoi(argv[++i]);
    }
  }

  //flens::verbose::ClosureLog::start("hodlr_log.txt");


  // Define types from FLENS
  using namespace flens;
  //using ZeroBased      = IndexOptions<int, 0>;
  using T              = double;
  //using T              = std::complex<double>;
  using MatrixType     = GeMatrix<FullStorage<T, ColMajor> >;
  using VectorType     = DenseVector<Array<T> >;
  using IndexType      = typename MatrixType::IndexType;
  using IndexVector    = DenseVector<Array<IndexType> >;
  const Underscore<IndexType> _;

  // Initialize the sources/targets as random values
  std::vector<double> source = fmmtl::random_n(N);
  //for (auto& s : source) s *= s;
  //std::sort(source.begin(), source.end());

  //auto polynomial = [](double x) { return 1 / std::sqrt(1e-2 + x*x); };

  // Create a test matrix
  MatrixType A(N,N);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = i+1; j <= N; ++j)
      //A(i,j) = 0;
      //A(i,j) = 1;
      A(i,j) = std::exp(-norm_2_sq(source[i-1] - source[j-1]));
      //A(i,j) = polynomial(norm_2(source[i-1] - source[j-1]));
      //A(i,j) = std::abs(source[i-1] - source[j-1]);
      //A(i,j) = std::exp(std::complex<double>(0,i*j*6.28318530718/N));
      //A(i,j) = std::exp(-norm_2_sq(std::sin(6.28*(source[i-1]-source[j-1]))));
  A.diag(0) = 1;

  A.lower() = transpose(A.upper());

  // Initialize a random RHS,   A*X = B
  MatrixType B(N,K);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= K; ++j)
      B(i,j) = fmmtl::random<T>::get();

  MatrixType testR, exactR;
  { ScopeClock timer("Direct Matvec: ");
    exactR = A * B;
  }


  //auto H = plr(A);

  auto H = geplr<1,1>('C',
                      A.data(), A.numRows(), A.numCols(), A.leadingDimension(),
                      source.data(), source.data(),
                      10, 1e-10, 2);

  { ScopeClock timer("PLR Matvec:    ");
    testR = H * B;
  }


  MatrixType ResR = exactR - testR;
  std::cout << "Matvec rel norm_F = " << norm_f(ResR)/norm_f(exactR) << std::endl;

  return 0;
}
