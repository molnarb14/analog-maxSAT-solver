#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil-g++.h"
#define NR_END 1

void nrerror(const char error_text[]) // Numerical Recipes standard error handler
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n", error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh) // allocate a float vector with subscript range v[nl..nh]
{
    float *v;
    v = new float[nh + NR_END]();
    if (v == nullptr) nrerror("allocation failure in vector()");
    return v;
}

int *ivector(long nl, long nh) // allocate an int vector with subscript range v[nl..nh]
{
    int *v;
    v = new int[nh + NR_END]();
    if (v == nullptr) nrerror("allocation failure in ivector()");
    return v;
}


long long *llvector(long nl, long nh) // allocate an int vector with subscript range v[nl..nh]
{
	long long *v;
    v = new long long[nh + NR_END]();
	if (v == nullptr) nrerror("allocation failure in ivector()");
	return v;
}

unsigned char *cvector(long nl, long nh) // allocate an unsigned char vector with subscript range v[nl..nh]
{
    unsigned char *v;
    v = new unsigned char[nh + NR_END]();
    if (v == nullptr) nrerror("allocation failure in cvector()");
    return v;
}

unsigned long *lvector(long nl, long nh) // allocate an unsigned long vector with subscript range v[nl..nh]
{
    unsigned long *v;
    v = new unsigned long[nh + NR_END]();
    if (v == nullptr) nrerror("allocation failure in lvector()");
    return v;
}

double *dvector(long nl, long nh) // allocate a double vector with subscript range v[nl..nh]
{
    double *v;
    v = new double[nh + NR_END]();
    if (v == nullptr) nrerror("allocation failure in dvector()");
    return v;
}

float **matrix(long nrl, long nrh, long ncl, long nch) // allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
{
    float **m;
    m = new float*[nrh + NR_END]();
    if(m == nullptr) nrerror("allocation failure 1 in matrix()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        m[i] = new float[nch + NR_END]();
        if(m[i] == nullptr) nrerror("allocation failure 2 in matrix()");
    }
    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch) // allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
{
    double **m;
    m = new double*[nrh + NR_END]();
    if(m == nullptr) nrerror("allocation failure 1 in dmatrix()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        m[i] = new double[nch + NR_END]();
        if(m[i] == nullptr) nrerror("allocation failure 2 in dmatrix()");
    }
    return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch) // allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
{
    int **m;
    m = new int*[nrh + NR_END]();
    if(m == nullptr) nrerror("allocation failure 1 in imatrix()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        m[i] = new int[nch + NR_END]();
        if(m[i] == nullptr) nrerror("allocation failure 2 in imatrix()");
    }
    return m;
}

long long **llmatrix(long nrl, long nrh, long ncl, long nch) // allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]
{
    long long **m;
    m = new long long*[nrh + NR_END]();
    if(m == nullptr) nrerror("allocation failure 1 in llmatrix()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        m[i] = new long long[nch + NR_END]();
        if(m[i] == nullptr) nrerror("allocation failure 2 in llmatrix()");
    }
    return m;
}

//TODO: needs to be transformed into C++ form
/*float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl) // point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch]
{
    long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
    float **m;
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;
    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
    return m;
}*/

//TODO: need to be transformed into C++ form
/*float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch) // allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1 and ncol=nch-ncl+1. The routine should be called with the address &a[0][0] as the first argument.
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;
    m[nrl]=a-ncl;
    for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
    return m;
}*/

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) // allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
{
    float ***t;
    t = new float**[nrh + NR_END]();
    if (t == nullptr) nrerror("allocation failure 1 in f3tensor()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        t[i] = new float*[nch + NR_END]();
        if(t[i] == nullptr) nrerror("allocation failure 2 in f3tensor()");
    }
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        for(int j = 0; i < nch + NR_END; ++j)
        {
            t[i][j] = new float[ndh]();
            if(t[i][j] == nullptr) nrerror("allocation failure 3 in f3tensor()");
        }
    }
    return t;
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) // allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]
{
    int ***t;
    t = new int**[nrh + NR_END]();
    if (t == nullptr) nrerror("allocation failure 1 in f3tensor()");
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        t[i] = new int*[nch + NR_END]();
        if(t[i] == nullptr) nrerror("allocation failure 2 in f3tensor()");
    }
    for (int i = 0; i < nrh + NR_END; ++i)
    {
        for(int j = 0; i < nch + NR_END; ++j)
        {
            t[i][j] = new int[ndh]();
            if(t[i][j] == nullptr) nrerror("allocation failure 3 in f3tensor()");
        }
    }
    return t;
}

void free_vector(float *v, long nl, long nh) // free a float vector allocated with vector()
{
    delete[] v;
}

void free_ivector(int *v, long nl, long nh) // free an int vector allocated with ivector()
{
    delete[] v;
}

void free_llvector(long long *v, long nl, long nh) // free an int vector allocated with ivector()
{
	delete[] v;
}

void free_cvector(unsigned char *v, long nl, long nh) // free an unsigned char vector allocated with cvector()
{
    delete[] v;
}

void free_lvector(unsigned long *v, long nl, long nh) // free an unsigned long vector allocated with lvector()
{
    delete[] v;
}

void free_dvector(double *v, long nl, long nh) // free a double vector allocated with dvector()
{
    delete[] v;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) // free a float matrix allocated by matrix()
{
    for(int i = 0; i < nrh; ++i)
        delete[] m[i];
    delete[] m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) // free a double matrix allocated by dmatrix()
{
    for(int i = 0; i < nrh; ++i)
        delete[] m[i];
    delete[] m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) // free an int matrix allocated by imatrix()
{
    for(int i = 0; i < nrh; ++i)
        delete[] m[i];
    delete[] m;
}

void free_llmatrix(long long **m, long nrl, long nrh, long ncl, long nch) // free an int matrix allocated by imatrix()
{
    for(int i = 0; i < nrh; ++i)
        delete[] m[i];
    delete[] m;
}

//TODO: rewrite in C++
/*void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch) // free a submatrix allocated by submatrix()
{
    free((FREE_ARG) (b+nrl-NR_END));
}*/

//TODO: rewrite in C++
/*void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch) // free a matrix allocated by convert_matrix()
{
    free((FREE_ARG) (b+nrl-NR_END));
}*/

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh) // free a float f3tensor allocated by f3tensor()
{
    for(int i = 0; i < nrh; ++i)
    {
        for(int j = 0; i < nch; j++)
        {
            delete[] t[i][j];
        }
    }
    for(int i = 0; i < nrh; ++i)
    {
        delete[] t[i];
    }
    delete[] t;
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh) // free a float f3tensor allocated by f3tensor()
{
    for(int i = 0; i < nrh; ++i)
    {
        for(int j = 0; i < nch; j++)
        {
            delete[] t[i][j];
        }
    }
    for(int i = 0; i < nrh; ++i)
    {
        delete[] t[i];
    }
    delete[] t;
}
