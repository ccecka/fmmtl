#include "HODLR.hpp"


int main(int argc, char** argv) {

  unsigned N = 1 << 12;     // rows
  unsigned leaf_size = 64;  // maximum size of the tree leaves

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    }
  }

  unsigned M = N;

  // Define types from FLENS
  using namespace flens;
  //using ZeroBased      = IndexOptions<int, 0>;
  using MatrixType     = GeMatrix<FullStorage<double, ColMajor> >;
  using VectorType     = DenseVector<Array<double> >;
  using IndexType      = typename MatrixType::IndexType;
  using IndexVector    = DenseVector<Array<IndexType> >;
  const Underscore<IndexType> _;

  // Initialize the sources/targets as random values
  std::vector<double> source = fmmtl::random_n(M);
  for (auto& s : source) s *= s;
  std::sort(source.begin(), source.end());

  // Create a test matrix
  MatrixType A(N,M);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= M; ++j)
      //A(i,j) = 1;
      //A(i,j) = std::exp(-norm_2_sq(p_targets[i-1] - p_sources[j-1]));
      A(i,j) = std::exp(-norm_2_sq(std::sin(6.28*(source[i-1]-source[j-1]))));
  A.diag(0) = 2;

  // Initialize a random RHS,   A*X = B
  unsigned Nb = N;
  MatrixType B(N,Nb);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= Nb; ++j)
      B(i,j) = fmmtl::random<double>::get();


  // Ge MATRIX
  {
  Clock timer;
  auto H = gehodlr('C', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;

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
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }
  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }

  A.lower() = transpose(A.upper());

  // Sy MATRIX
  {
  Clock timer;
  auto H = syhodlr('C', 'U', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;


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
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }

  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }

  // He MATRIX
  {
  Clock timer;
  auto H = hehodlr('C', 'U', A.data(), A.numRows(), A.leadingDimension(), leaf_size);
  std::cout << "HODLR Construction: " << timer.seconds() << std::endl;


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
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::sv(H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }
  testX = B;
  { ScopeClock timer("HODLR  Solve: ");
    flens::lapack::trs(NoTrans, H, ipiv, testX);
  }

  MatrixType Acpy = A;
  MatrixType exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::sv(Acpy, ipiv, exactX);
  }
  exactX = B;
  { ScopeClock timer(std::string("Direct Solve: "));
    flens::lapack::trs(NoTrans, Acpy, ipiv, exactX);
  }

  MatrixType ResX = exactX - testX;
  std::cout << "Solve  rel norm_F = " << norm_f(ResX)/norm_f(exactX) << std::endl;
  }

}
