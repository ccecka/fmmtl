#pragma once
/** @file HelmholtzQuad.hpp
 * @brief Define the Quadratures to use in the Helmholtz kernel
 */

#include "Bessel.hpp"
#include "FFT.hpp"

#include <iostream>
#include <cmath>
#include <complex>
#include <numeric>

// Ceil/Round/Floor to nearest multiple of K
inline int ceil (double x, unsigned K) { return K*((int)ceil(x/K));  }
inline int round(double x, unsigned K) { return K*((int)round(x/K)); }
inline int floor(double x, unsigned K) { return K*((int)floor(x/K)); }

typedef Vec<3,double> Vec3;
typedef std::complex<double> complex;

constexpr complex CI(0,1);
constexpr double  PI = M_PI;

// Overloaded complex output
inline std::ostream& operator<<(std::ostream& os, complex a)
{
  //ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
  //ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

  int oldp = os.precision(6);

  os << std::real(a);
  if( std::imag(a) != 0 ) {
    if( std::imag(a) < 0 )
      os << " - " << -std::imag(a) << "*i";
    else
      os << " + " << std::imag(a) << "*i";
  }

  //os.setf(olda,ios::adjustfield);
  //os.setf(oldf,ios::floatfield);
  os.precision(oldp);
  os << "";

  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& a)
{
  std::ios::fmtflags olda = os.setf(std::ios::right,std::ios::adjustfield);
  std::ios::fmtflags oldf = os.setf(std::ios::scientific,std::ios::floatfield);

  int oldp = os.precision(8);

  int N = a.size();
  for( int k = 0; k < N; ++k ) {
    os << k << "\t" << a[k] << "\n";
  }

  os.setf(olda,std::ios::adjustfield);
  os.setf(oldf,std::ios::floatfield);
  os.precision(oldp);
  os << "";

  return os;
}




inline complex gegenbauer_series(int ell, double kappa,
                                 const Vec3& r, const Vec3& r0)
{
  double kr0 = kappa * norm(r0);
  double kr  = kappa * norm(r);
  double rr0 = r0.dot(r) / (norm(r0) * norm(r));

  complex G = 0;
  for (int n = 0; n <= ell; ++n) {
    G += double(2*n+1)*bessel_h(n,kr0)*bessel_j(n,kr)*legendre_P(n,rr0);
    ++n;
    if (n > ell) break;
    G -= double(2*n+1)*bessel_h(n,kr0)*bessel_j(n,kr)*legendre_P(n,rr0);
  }

  return CI * kappa * G;
}

// Get the Gegenbauer truncation, L, via a number of possible methods
inline static int gegenbauer_truncation(double kappa,
                                        double norm_r, double norm_r0,
                                        double eps,
                                        int method = 4)
{
  int ell = 0;

  switch (method) {
    case 1: {// Old Log Method
      ell = (int) (kappa*norm_r - log10(eps)*log(PI + kappa*norm_r));
    } break;
    case 2: { // EBF Formula (Chew)
      double knormr = kappa * norm_r;
      ell = ceil(knormr + 1.8*pow(-log10(eps), 2.0/3.0) * pow(knormr, 1.0/3.0));
    } break;
    case 3: { // Direct method (Collino) with worst case cos = -1
      // Compute the EBF_L as an upperbound
      double knorm_r  = kappa*norm_r;
      double knorm_r0 = kappa*norm_r0;

      // A slightly more accurate EBF ell to use as a starting point
      ell = gegenbauer_truncation(kappa, norm_r, norm_r0, 1e-2 * eps, 2);

      double eM = knorm_r*knorm_r0 / std::abs(norm_r0-norm_r);
      double eP = knorm_r*knorm_r0 / std::abs(norm_r0+norm_r);

      complex hl = bessel_h(ell, knorm_r0), hlp1 = bessel_h(ell+1, knorm_r0);
      double  jl = bessel_j(ell, knorm_r ), jlp1 = bessel_j(ell+1, knorm_r );

      // Decrease ell until we hit the cutoff
      double error = std::max(eM*std::abs(hlp1*jl-hl*jlp1), eP*std::abs(hlp1*jl+hl*jlp1));
      //cout << "Trunc Error: " << error << "    ell: " << ell << endl;
      while (error < eps && ell > 1) {
        --ell;
        hlp1 = hl;  hl = bessel_h( ell, knorm_r0 );
        jlp1 = jl;  jl = bessel_j( ell, knorm_r  );
        error = std::max(eM*std::abs(hlp1*jl-hl*jlp1), eP*std::abs(hlp1*jl+hl*jlp1));
        //cerr << "Trunc Error: " << error << "    ell: " << ell << endl;
      }
      // Went one step over
      ++ell;

      // Check this against the true Gegenbauer error
      // If it agrees, we're good
      // If it doesn't agree, then |r| is too small and we should
      //          do an ultradirect check (using the true Geg error)

    } break;
    case 4: { // Ultradirect method - check against gegenbauer series itself
      // Use the EBF ell as a starting point
      ell = gegenbauer_truncation(kappa, norm_r, norm_r0, eps, 2);

      double max_error = 1e100;
      while (max_error > eps) {
        ++ell;
        max_error = 0;

        // For 20 directions of r0
        int N = 20;
        for (int n = 0; n <= N; ++n) {
          Vec3 r0 = Vec3(std::sin((PI*n)/N),
                         0,
                         std::cos((PI*n)/N)) * norm_r0;
          Vec3 r  = Vec3(0,0,1) * norm_r;

          double R = norm(r+r0);
          complex Ie = exp(CI * kappa* R) / R;

          complex GS = gegenbauer_series(ell, kappa, r, r0);

          //cout << Ie << "\t" << GS << "\t" << abs(Ie-GS) << endl;

          max_error = std::max(max_error, std::abs(Ie - GS));
        }

        //cout << ell << ": " << max_error << endl;
      }

      //cout << "Truncation - EBF: " << ell_ebf << "  UD: " << ell << endl;
      //ell = min(ell,ell_ebf);

    } break;
  }
  return ell;
}

// A class to quickly evaluate the transfer function T_{ell,r_0}(s)
// for many possible directions s
class Transfer_Function_Eval
{
  Vec3 r0hat;
  std::vector<complex> H;   // Coefficients of the series
  std::vector<double> P;
 public:
  // Constructor
  Transfer_Function_Eval(double kappa, int ell, const Vec3& r0)
      : r0hat(r0/norm(r0)), H(ell+1), P(ell+1) {
    double knormr0 = kappa * norm(r0);
    complex i_k = (CI*kappa)/(4*PI);
    for (unsigned n = 0; n < H.size(); ++n) {
      H[n] = i_k * double(2*n+1) * bessel_h(n, knormr0);
      i_k *= CI;
    }
  }

  inline complex operator()(double sdotr0) {
    gsl_sf_legendre_Pl_array(P.size()-1, sdotr0, P.data());
    return std::inner_product(H.begin(), H.end(), P.begin(), complex(0,0));
  }

  inline complex operator()(const Vec3& s) {
    assert(std::abs(norm(s)-1) < 1e-15);
    return operator()(r0hat.dot(s));
  }
};


/* Computes a low-pass |sin(phi)|
 *
 * nF: number of desired frequencies of abs(sin)
 * nR: number of desired real space points of abs(sin), nR >= 2*nF+1
 */
inline std::vector<complex> FabsSin( int nF, int nR )
{
  assert(nR >= 2*nF+1);
  std::vector<complex> abssin(nR,0);

  // Compute all the Fourier coefficients
  abssin[0] = 2.0/PI;
  for( int k = 2; k <= nF; k += 2 )
    abssin[nR - k] = abssin[k] = 2.0/(PI*(1 - k*k));

  fft(abssin.data(), nR);
  return abssin;
}


// Computes the low-pass modified transfer function matrix N_rows x N_phi
// ell:  Gegenbauer truncation
// kappa: Wavenumber
// N_phi: the number of quadrature points in phi
// N_rows: the number of rows of the quadrature N_theta = 2*N_rows - 2
//
// Returns: T[n][m] nth row, mth phi value of the low-pass modified trans fn
inline std::vector<std::vector<complex> > mod_transfer_fn(int ell,
                                                          double kappa,
                                                          int N_rows,
                                                          int N_phi,
                                                          const Vec3& r0)
{
  assert(!(N_phi & 1)); // N_phi is even

  Vec3 r0hat = r0 / norm(r0);

  // Construct a Transfer_Function_Eval T_{ell,r_0}
  Transfer_Function_Eval Tfn(kappa, ell, r0);

  int N_theta = 2*N_rows - 2;           // # points used for ModTransFn
  int N_t = 2*ell + 1;                    // # points needed for TransFn

  // Get the low-pass |sin(theta)|
  int N_sinF = (N_theta/2-1) + ell;       // # freqs needed in |sin|
  //int N_conv = 2*(ell+N_sinF) + 1;
  int N_conv = N_theta + 2*ell - 1;       // size of convolution with |sin|
  std::vector<complex> abssin = FabsSin(N_sinF, N_conv);

  // The N_rows x N_phi transfer function
  std::vector<std::vector<complex>> T(N_rows, std::vector<complex>(N_phi,0));

  // Temp storage for t_phi(theta)
  std::vector<complex> t(N_conv);

  // FFTs for t(theta)    (Optimized)
  fftw_complex* a = reinterpret_cast<fftw_complex*>(&t[0]);
  fftw_plan tFFT1 = fftw_plan_dft_1d(N_t,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan tIFFT1 = fftw_plan_dft_1d(N_conv,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_plan tFFT2 = fftw_plan_dft_1d(N_conv,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan tIFFT2 = fftw_plan_dft_1d(N_theta,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);

  const double scale = 0.5 / (N_t*N_conv);    // 1/2 and interp/anterp coeffs

  // For each phi [0,PI)
  for (int m = 0; m < N_phi/2; ++m) {
    double phi = m * (2*PI)/N_phi;
    double cosP = std::cos(phi);
    double sinP = std::sin(phi);
    double xydotr0 = cosP*r0hat[0] + sinP*r0hat[1];

    // Sample the transfer function at 2ell+1 points in theta [0,2PI)
    for (int n = 0; n < N_t; ++n) {
      double theta = n * (2*PI)/N_t;
      double sinT = std::sin(theta);
      double cosT = std::cos(theta);

      //Vec3 s( sintheta*cosphi, sintheta*sinphi, costheta );
      //t[n] = scale * transfer_function(ell, kappa, r0, s);
      double sdotr0 = sinT*xydotr0 + cosT*r0hat[2];
      t[n] = scale * Tfn(sdotr0);    // Optimized
    }

    // Interpolate the transfer function to N_conv
    fftw_execute(tFFT1);
    f_cut(t.data(), N_t, N_conv);
    fftw_execute(tIFFT1);

    // Multiply by low-pass |sin| to get the modified transfer function
    for (int k = 0; k < N_conv; ++k)
      t[k] *= abssin[k];

    // Anterpolate to N_theta
    fftw_execute(tFFT2);
    f_cut(t.data(), N_conv, N_theta);
    t[N_theta/2] = 0;                    // Set the oddball freq to zero
    fftw_execute(tIFFT2);

    // Explicit convolution to check the answer
    //std::vector<complex> exact = slowConv(t, ell, N_theta/2-1);

    // Unwrap into the N_rows x N_phi matrix
    T[0][m] = T[0][N_phi/2+m] = t[0];       // North Poles
    for (int n = 1; n < N_rows; ++n) {
      T[n][m]         = t[n];
      T[n][N_phi/2+m] = t[N_theta-n];
    }
  }

  fftw_destroy_plan(tFFT1);
  fftw_destroy_plan(tIFFT1);
  fftw_destroy_plan(tFFT2);
  fftw_destroy_plan(tIFFT2);

  return T;
}






inline int get_N_theta(int ell, double kappa,
                       double norm_r, double norm_r0,
                       double eps)
{
  int N_rows_max = ell + 30;
  int N_theta_max = 2*N_rows_max - 2;

  // EBF METHOD
  /*
    int maxL_ebf = ell_ebf(kappa, norm_r, eps);
    int N_theta_ebf = 2*maxL_ebf + 2;
    cerr << "N_theta_ebf = " << N_theta_ebf << endl;
    return N_theta_ebf;
  */

  // Another Heuristic
  //return 2 * L + 2;

  // Get the transfer matrix for phi = 0, theta = [0,2PI)
  Vec3 r0(0,0,norm_r0);
  std::vector<std::vector<complex>> T = mod_transfer_fn(ell,kappa,N_rows_max,2,r0);

  // Copy real-space values into a single std::vector t(theta)
  std::vector<complex> t(N_theta_max);
  t[0] = T[0][0];                        // North Pole
  for (int k = 1; k < N_rows_max; ++k) {
    t[k] = T[k][0];
    t[N_theta_max-k] = T[k][1];
  }

  // Get the Fourier coefficients
  ifft(t.data(), N_theta_max);

  // Get the absolute value of the Fourier coefficients
  std::vector<double> absT( N_theta_max/2-1 );
  for (int k = 0; k < N_theta_max/2-1; ++k) {
    absT[k] = std::abs(t[k]);
    // Make sure this is correct and symmetric
    //cout << k << ": " << abs( abs(t[k])-abs(t[N_theta_max-k]) ) << endl;
    //assert( k == 0 || abs( abs(t[k])-abs(t[N_theta_max-k]) ) < 1e-12 );
  }

  // Compute the bessel functions
  std::vector<double> absJ(N_theta_max + 1);
  for (int k = 0; k <= N_theta_max; ++k) {
    absJ[k] = std::abs(bessel_J(k, kappa*norm_r));
    //std::cerr << "|J_" << k << "(" << kappa*norm_r << ")| = " << absJ[k] << std::endl;
  }

  // Search for N_theta
  for (int N_theta = std::max(2*(ell-10),0); N_theta < N_theta_max; N_theta+=2) {

    // Get the truncation error sum_{n >= N_theta/2} + sum_{n <= -N_theta/2}
    double error_trunc = 0;
    for( int n = N_theta/2; n < N_theta_max/2-1; ++n )
      error_trunc += absT[n] * absJ[n];

    // Get the aliasing error
    double error_alias = absT[0] * absJ[N_theta];
    for( int n = 1; n <= N_theta/2-1; ++n )
      error_alias += absT[n] * absJ[N_theta - n];

    double error = 4*PI*PI * 2*(error_alias + error_trunc);
    //std::cerr << N_theta << ": " << error << std::endl;

    if (error < eps) return N_theta;
  }

  std::cerr << "Error: N_theta convergence failure" << std::endl;
  exit(1);
}




struct S2Quadrature {
  //! Helmholtz wavenumber
  double kappa;
  //! Level in the tree
  int level;
  //! Gegenbauer truncation
  int ell;
  //! Maximum number of points in a Quadrature row (constant theta)
  unsigned max_N_phi;
  //! Row offsets. kth row: Indices [row_offset[k], row_offset[k+1]) into q
  std::vector<unsigned> row_offset;

  struct S2QuadraturePoint {
    double theta, phi, w;
    Vec<3,double> s;
    S2QuadraturePoint(double theta_ = 0, double phi_ = 0, double w_ = 0)
        : theta(theta_), phi(phi_), w(w_),
          s(std::sin(theta)*std::cos(phi),
            std::sin(theta)*std::sin(phi),
            std::cos(theta)) {}
  };

  //! Quadrature points
  std::vector<S2QuadraturePoint> q;

  S2Quadrature(double _kappa, int _level,
               double box_size, double alpha, double eps)
      : kappa(_kappa), level(_level) {
    // Smallest transfer vector magnitude (worst case)
    double norm_r0 = 2 * box_size;
    // Largest translation vector magnitude (worst case)
    double norm_r = alpha * sqrt(3) * box_size;

    ell = gegenbauer_truncation(kappa, norm_r, norm_r0, eps);

    std::cerr << "norm_r0 = " << norm_r0 << std::endl;
    std::cerr << "norm_r = " << norm_r << std::endl;
    std::cerr << "ell = " << ell << std::endl;

    /* Determine the number of quadrature rows (N_rows = N_theta/2 + 1) */
    int N_theta = get_N_theta(ell, kappa, norm_r, norm_r0, eps);
    int N_rows = N_theta/2 + 1;
    int N_phiT = 2*ell + 2;
    row_offset.resize(N_rows+1, 0);
    row_offset[0] = 0;

    Vec3 r0(norm_r0, 0, 0);
    std::vector<std::vector<complex>> TL = mod_transfer_fn(ell,kappa,N_rows,N_phiT,r0);

    // Prune N_phi(theta) if possible
    for (int n = 0; n < N_rows; ++n) {
      double theta = n*PI/(N_rows-1);
      int N_phi = 0;

      if (n == 0 || n == N_rows-1) {   // The Poles
        N_phi = 1;
      } else if(n <= N_rows/2) {      // Positive z row
        // Determine N_phi for this row
        double knormrs = kappa * norm_r * std::sin(theta);

        // EBF METHOD (This is pretty good...)
        /*
          int maxL_ebf = ceil(knormrs + 1.8*pow(-log10(eps), 2.0/3.0)
          *pow(knormrs, 1.0/3.0));
          int N_phi_ebf = round(2*maxL_ebf + 1, 4);
          int N_phi = N_phi_ebf;
        */
        // Another Heuristic
        //int N_phi = round(2*L+1,4);

        // Compute the transfer function coefficients
        ifft(TL[n].data(), TL[n].size());

        std::vector<double> absT(ell+1);
        for (int m = 0; m <= ell; ++m) {
          absT[m] = std::abs(TL[n][m]);
          // Make sure this is correct and symmetric
          //assert
        }

        // Precompute Bessel J
        int max_ell = ell + 20;
        std::vector<double> absJ(ceil(2*max_ell+1, 4) + 1);
        for( int k = 0; k < (int) absJ.size(); ++k ) {
          absJ[k] = std::abs(bessel_J(k, knormrs));
        }

        // Compute N_phi by searching for error
        int N_phi_max = ceil(2*max_ell+8, 4);
        for (N_phi = 4; N_phi <= N_phi_max; N_phi += 4) {
          //for( N_phi = ceil(2*maxL+1, 4); N_phi > 4; N_phi -= 4 ) {

          // Get the truncation error sum_{|k| >= N_phi/2}
          double error_trunc = 0;
          for (int m = N_phi/2; m <= ell; ++m) {
            error_trunc += absT[m] * absJ[m];
          }

          // Get the aliasing error
          double error_alias = absT[0] * absJ[N_phi];
          for (int m = 1; m <= N_phi/2-1; ++m) {
            error_alias += absT[m] * absJ[N_phi - m];
          }

          double error = 4*PI*PI * 2*( error_alias + error_trunc );
          //cerr << N_phi << ": " << error << "\t" << error_alias << "\t" << error_trunc << endl;

          if (error < eps) break;
          //if( error > eps ) { N_phi += 4; break; }
        }

        if (N_phi >= N_phi_max)
          std::cerr << "WARNING: N_phi ceiling hit" << std::endl;
      } else {                         // Negative z row
        // Copy the points from positive z row
        N_phi = n_phi(N_rows-1-n);
      }

      // Save the N_phi for this row
      assert(N_phi >= 1);
      row_offset[n+1] = row_offset[n] + N_phi;
    }

    // Construct the quadrature points
    q.resize(row_offset.back());
    // Poles
    q.front() = S2QuadraturePoint(0, 0, 4*PI*PI/n_theta());
    q.back()  = S2QuadraturePoint(PI, 0, 4*PI*PI/n_theta());
    max_N_phi = 1;
    // Other Rows
    for (unsigned n = 1; n < n_row()-1; ++n) {
      max_N_phi = std::max(max_N_phi, n_phi(n));
      double theta = n*PI/(n_row()-1);
      for (unsigned k = 0; k < n_phi(n); ++k) {
        double phi = k * 2*PI/n_phi(n);
        double w = 2 * 2*PI*2*PI/(n_theta()*n_phi(n));
        q[row_offset[n] + k] = S2QuadraturePoint(theta, phi, w);
      }
    }
  }
  Vec<3,double>& operator[](int k) {
    return q[k].s;
  }
  const Vec<3,double>& operator[](int k) const {
    return q[k].s;
  }
  unsigned size() const {
    return q.size();
  }
  /** The number of points in theta [0,2pi) (polar angle) */
  unsigned n_theta() const {
    return 2*n_row() - 2; // Don't double-count the poles
  }
  /** The number of rows (in theta [0,pi])*/
  unsigned n_row() const {
    return row_offset.size() - 1;
  }
  /** The maximum number of points in phi [0,2pi) (azimuthal angle) */
  unsigned n_phi() const {
    return max_N_phi;
  }
  /** The number of points (in phi) in row k */
  unsigned n_phi(int k) const {
    return row_offset[k+1] - row_offset[k];
  }
  double weight(int k) const {
    return q[k].w;
  }
  unsigned row_begin(int k) const {
    return row_offset[k];
  }
  unsigned row_end(int k) const {
    return row_offset[k+1];
  }
  friend std::ostream& operator<<(std::ostream& os, const S2Quadrature& q) {
    os << "(" << q.ell
       << "," << q.n_row()
       << "," << q.size() << ")" << std::endl;
    for (unsigned k = 0; k < q.n_row(); ++k) {
      os << "\t" << q.n_phi(k) << std::endl;
    }
    return os;
  }
};


struct S2Function {
  S2Quadrature* quad;
  std::vector<complex> c;
  S2Function()
      : quad(nullptr) {
  }
  S2Function(S2Quadrature* _quad)
      : quad(_quad), c(quad->size()) {
  }
  unsigned size() const {
    return c.size();
  }
  complex& operator[](int k) {
    return c[k];
  }
  const complex& operator[](int k) const {
    return c[k];
  }
  std::vector<complex>::iterator row_begin(int k) {
    return c.begin() + quad->row_begin(k);
  }
  std::vector<complex>::const_iterator row_begin(int k) const {
    return c.begin() + quad->row_begin(k);
  }
  std::vector<complex>::iterator row_end(int k) {
    return c.begin() + quad->row_end(k);
  }
  std::vector<complex>::const_iterator row_end(int k) const {
    return c.begin() + quad->row_end(k);
  }
  friend std::ostream& operator<<(std::ostream& os, const S2Function& A) {
    return os << A.c;
  }
};



struct Transfer_Function : public S2Function {
  Transfer_Function(S2Quadrature* quad, const Vec3& r0)
      : S2Function(quad)
  {
    //cerr << "r0 = " << r0 << endl;

    double kappa = quad->kappa;
    unsigned ell = quad->ell;
    unsigned N_rows = quad->n_row();
    unsigned N_phi = 2*ell + 2;

    // Get the low-pass modified transfer matrix N_rows x N_phi
    std::vector<std::vector<complex>> T = mod_transfer_fn(ell,kappa,N_rows,N_phi,r0);

    // maxRowSize can be larger than Nphi -> shouldn't interpolate in place
    std::vector<complex> t(std::max(quad->n_phi(),N_phi));

    // Anterpolate the rows for an optimized quadrature
    for (unsigned n = 0; n < N_rows; ++n) {
      int N_phi_n = quad->n_phi(n);

      // Anterpolate to N_phi_k
      std::copy(T[n].begin(), T[n].begin() + N_phi, t.begin());
      fftinterp(t.data(), N_phi, N_phi_n);

      // Copy into the NFunction
      std::copy(t.begin(), t.begin() + N_phi_n, row_begin(n));
    }
  }
};
