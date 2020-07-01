/*---------------------------------------------------------------------*/
/*! \file

\brief build options and definition of fortran functions

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef COMPILER_DEFINITIONS_H
#define COMPILER_DEFINITIONS_H


/*----------------------------------------------------------------------*
 | special definitions for special compilers.....                       |
 *----------------------------------------------------------------------*/
/* append underslashs, if necessary. Important for linking to fortran routines! */

#undef CCA_APPEND_U

/* append underslash for gnu's linux compiler gcc and g77 */
/* refer to src/fortran for the respective routines */

#define CCA_APPEND_U (1)

#ifdef CCA_APPEND_U
/* required for lapack access. Do not remove! */
#define dgesv dgesv_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define dgetrs dgetrs_
#define dggev dggev_
#define dsyev dsyev_
#define dsyevd dsyevd_
#define dsygv dsygv_
#define dsytrf dsytrf_
#define dsytri dsytri_
#define dsytrs dsytrs_
#define dveczero dveczero_
/* required in shell8 element. Do not remove! */
#define fortranpow fortranpow_
#define s8jacb s8jacb_
#define dhgeqz dhgeqz_
#define dgghrd dgghrd_
#define dgeqp3 dgeqp3_
#define dggbal dggbal_

#endif

/* fortran routines from the lapack package, used in
 * src/linalg/linalg_utils_densematrix_determinant.cpp,
 * src/linalg/linalg_utils_densematrix_eigen.cpp and src/linalg/linalg_utils_densematrix_svd.cpp */
#ifdef __cplusplus
extern "C"
{
#endif

  void dsytrf(
      char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
  void dgetrf(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
  void dgetri(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dggev(const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda,
      double *b, const int *ldb, double *alphar, double *alphai, double *beta, double *vl,
      const int *ldvl, double *vr, const int *ldvr, double *work, const int *lwork, int *info);
  void dhgeqz(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, double *h, int *ldh,
      double *t, int *ldt, double *alphar, double *alphai, double *beta, double *q, int *ldq,
      double *z, int *ldz, double *work, int *lwork, int *info);
  void dgghrd(char *compq, char *compz, int *n, int *ilo, int *ihi, double *a, int *lda, double *b,
      int *ldb, double *q, int *ldq, double *z, int *lzd, int *info);
  void dgeqp3(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork,
      int *info);
  void dggbal(const char *job, const int *n, double *A, const int *lda, double *B, const int *ldb,
      int *ilo, int *ihi, double *lscale, double *rscale, double *work, int *info);

#ifdef __cplusplus
}

#endif

#endif /* COMPILER_DEFINITIONS_H */
