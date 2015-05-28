#pragma once
/** Simple C-library interfaces for PLR matrices
 */

#include <complex.h>

// An opaque type that we'll use as a handle
struct ZPLR2D_Handle;
typedef struct ZPLR2D_Handle ZPLR2D_Handle;


/** C-interface to the complex general matrix PLR decomposition
 *   over two-dimensional domains
 *
 * @param[in] order   'R': Row-major ordering. 'C': Column-major ordering
 * @param[in] data    Array containing the coefficients of the matrix.
 * @param[in] n       Number of rows of the matrix.
 * @param[in] m       Number of columns of the matrix
 * @param[in] lda     Leading dimension of the matrix, >= max(1,m).
 *                If order == 'R', data[i*lda+j] is the (i,j)-th element
 *                If order == 'C', data[i+lda*j] is the (i,j)-th element
 * @param[in] trgs    Coordinate-major points corresponding to rows of the matrix
 * @param[in] srcs    Coordinate-major points corresponding to cols of the matrix
 * @param[in] max_rank The maximum rank of a block in the PLR structure
 * @param[in] eps_tol  The maximum err tolerance of a block in the PLR structure
 * @param[in] init_depth Optimization parameter to start the decomposition tests
 *                       at a lower level of the tree hierarchy.
 *
 * @pre size(data) >= (lda, n)
 * @pre size(targets) == 2*n
 * @pre size(sources) == 2*m
 */
ZPLR2D_Handle*
zplr2d(char order, double _Complex* data, int n, int m, int lda,
       const double* trgs, const double* srcs,
       unsigned max_rank, double eps_tol, unsigned init_depth = 0);

/** y = beta*y + alpha*M*x */
void
zplr2dmv(char trans, const ZPLR2D_Handle* M, int n, int m,
         double _Complex alpha, const double _Complex* x, int incx,
         double _Complex beta, double _Complex* y, int incy);

/** Destruct a ZPLR2D Matrix and free its memory */
void
free_zplr2d(ZPLR2D_Handle* plr);
