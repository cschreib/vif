#ifndef VIF_MATH_LAPACK_HPP
#define VIF_MATH_LAPACK_HPP

// LAPACK functions imported from fortran library
// ----------------------------------------------

namespace lapack {
    extern "C" void dgetrf_(int* n, int* m, double* a, int* lda, int* ipiv, int* info);
    extern "C" void dgetri_(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork,
        int* info);
    extern "C" void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,
        int* ldb, int* info);
    extern "C" void dgesdd_(char* jobz, int* m, int* n, double* a, int* lda, double* s, double* u,
        int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* iwork, int* info);

    extern "C" void dsytrf_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work,
        int* lwork, int* info);
    extern "C" void dsytri_(char* uplo, int* n, double* a, int* lda, int* ipiv, double* work,
        int* info);
    extern "C" void dsysv_(char* uplo, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b,
        int* ldb, double* work, int* lwork, int* info);
    extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
        double* work, int* lwork, int* info);

    extern "C" void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
    extern "C" void dpotri_(char* uplo, int* n, double* a, int* lda, int* info);
    extern "C" void dpotrs_(char* uplo, int* n, int* nrhs, const double* a, int* lda, double* b,
        int* ldb, int* info);

    extern "C" void dtrtri_(char* uplo, char* diag, int* n, double* a, int* lda, int* info);
}

#endif
