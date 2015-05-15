//#define FLENS_DEBUG_CLOSURES

#include "HODLR.hpp"


int main(int argc, char** argv) {

  unsigned N = 1 << 12;     // rows
  unsigned K = 1;           // cols of the rhs
  unsigned leaf_size = 32;  // maximum size of the tree leaves

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
  using T              = double; // std::complex<double>;
  using MatrixType     = GeMatrix<FullStorage<T, ColMajor> >;
  using VectorType     = DenseVector<Array<T> >;
  using IndexType      = typename MatrixType::IndexType;
  using IndexVector    = DenseVector<Array<IndexType> >;
  const Underscore<IndexType> _;

  // Initialize the sources/targets as random values
  std::vector<double> source = fmmtl::random_n(N);
  //for (auto& s : source) s *= s;
  std::sort(source.begin(), source.end());

  auto polynomial = [](double x) { return 1 / std::sqrt(1e-2 + x*x); };

  // Create a test matrix
  MatrixType A(N,N);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = i+1; j <= N; ++j)
      //A(i,j) = 0;
      //A(i,j) = 1;
      //A(i,j) = std::exp(-norm_2_sq(source[i-1] - source[j-1]));
      A(i,j) = polynomial(norm_2(source[i-1] - source[j-1]));
      //A(i,j) = std::abs(source[i-1] - source[j-1]);
      //A(i,j) = std::exp(std::complex<double>(0,i*j*6.28318530718/N));
      //A(i,j) = std::exp(-norm_2_sq(std::sin(6.28*(source[i-1]-source[j-1]))));
  A.diag(0) = 2;

  A.lower() = transpose(A.upper());

  // Initialize a random RHS,   A*X = B
  MatrixType B(N,K);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= K; ++j)
      B(i,j) = fmmtl::random<T>::get();

#if 1
  // Ge MATRIX
  {
  std::cout << "N = " << N << std::endl;
  std::cout << "K = " << K << std::endl;

  Clock timer;
  auto H = gehodlr('C', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;

  std::cout << "HODLR Compression: " << H.compression() << std::endl;


  // Matmats
  MatrixType testR, exactR;
  { ScopeClock timer("HODLR  Matvec: ");
    testR = H * B;
  }

  { ScopeClock timer(std::string("Direct Matvec: "));
    exactR = A * B;
  }
  MatrixType ResR = exactR - testR;
  std::cout << "Matvec rel norm_F = " << norm_f(ResR)/norm_f(exactR) << std::endl;

  // Solves
  IndexVector ipiv;
  MatrixType testX = B;
  { ScopeClock timer("HODLR  Solve1: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve2: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve3: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  std::pair<T,int> det;
  { ScopeClock timer("HODLR Det: ");
    det = flens::lapack::det(H);
  }
  std::cout << "|H| = " << det.first << "e" << det.second << std::endl;

  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer("Direct Solve1: ");
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer("Direct Solve2: ");
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }

  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }
#endif


#if 0
  // Sy MATRIX
  {
  std::cout << "N = " << N << std::endl;
  std::cout << "K = " << K << std::endl;

  Clock timer;
  auto H = syhodlr('C', 'U', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;

  std::cout << "HODLR Compression: " << H.compression() << std::endl;

  // Matmats
  MatrixType testR, exactR;
  { ScopeClock timer("HODLR  Matvec: ");
    testR = H * B;
  }

  { ScopeClock timer(std::string("Direct Matvec: "));
    exactR = A * B;
  }
  MatrixType ResR = exactR - testR;
  std::cout << "Matvec rel norm_F = " << norm_f(ResR)/norm_f(exactR) << std::endl;

  // Solves
  IndexVector ipiv;
  MatrixType testX = B;
  { ScopeClock timer("HODLR  Solve1: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve2: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve3: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  std::pair<T,int> det;
  { ScopeClock timer("HODLR Det: ");
    det = flens::lapack::det(H);
  }
  std::cout << "|H| = " << det.first << "e" << det.second << std::endl;

  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer("Direct Solve1: ");
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer("Direct Solve2: ");
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }

  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }
#endif

#if 0
  // He MATRIX
  {
  std::cout << "N = " << N << std::endl;
  std::cout << "K = " << K << std::endl;

  Clock timer;
  auto H = hehodlr('C', 'U', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;

  std::cout << "HODLR Compression: " << H.compression() << std::endl;

  // Matmats
  MatrixType testR, exactR;
  { ScopeClock timer("HODLR  Matvec: ");
    testR = H * B;
  }

  { ScopeClock timer(std::string("Direct Matvec: "));
    exactR = A * B;
  }
  MatrixType ResR = exactR - testR;
  std::cout << "Matvec rel norm_F = " << norm_f(ResR)/norm_f(exactR) << std::endl;

  // Solves
  IndexVector ipiv;
  MatrixType testX = B;
  { ScopeClock timer("HODLR  Solve1: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve2: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve3: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  std::pair<T,int> det;
  { ScopeClock timer("HODLR Det: ");
    det = flens::lapack::det(H);
  }
  std::cout << "|H| = " << det.first << "e" << det.second << std::endl;


  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer("Direct Solve1: ");
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer("Direct Solve2: ");
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }

  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }
#endif

  //flens::verbose::ClosureLog::stop();
}
