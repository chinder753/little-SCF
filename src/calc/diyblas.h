#ifndef DIYBLAS_H
#define DIYBLAS_H

// BLAS
typedef enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142} CBLAS_SIDE;
void cblas_daxpy(const int n, const double alpha, const double *x, const int incx, double *y, const int incy);
void cblas_dcopy(const int n, const double *x, const int incx, double *y, const int incy);
double cblas_ddot(const int n, const double *x, const int incx, const double *y, const int incy);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_dswap(const int n, double *x, const int incx, double *y, const int incy);
double cblas_dnrm2 (const int N, const double *X, const int incX);
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const int M, const int N, const int K,
                 const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N, const double alpha, const double *Ap,
                 const double *X, const int incX, const double beta, double *Y, const int incY);

// LAPACKE
#define lapack_int     int
#define LAPACK_COL_MAJOR               102
lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w );
lapack_int LAPACKE_dspgvd( int matrix_layout, lapack_int itype, char jobz,
                           char uplo, lapack_int n, double* ap, double* bp,
                           double* w, double* z, lapack_int ldz );

void spprintf(const int n, const char str[], const double * const a);

void geprintf(const int n, const char str[], const double * const a);

int sp2sy(const int n
          , const double * const x
          , double * const y
          );

int sy2sp(const int n
          , const double * const x
          , double * const y
          );

int tran(const int n
        , const double * const x
        , double * const y
        );

#endif
