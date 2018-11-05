//v8.2.2BST
#include <math.h>
#include "nrutil.h"

#define SWAP(a, b) { swap = (a); (a) = (b); (b) = swap; }

void funcs(double x, double a[], double *y, double dyda[], int na, double fconst)
{
    //f(x) = fconst + a[1] * x^a[2]
    
    double r = pow(x, a[2]);
    
    *y = fconst + a[1] * r;
    
    dyda[1] = r;
    dyda[2] = a[1] * r * log(x);
}

int gaussj(double **a, int n, double **b, int m)
{
    int *indxc, *indxr, *ipiv;
    int i, icol, irow, j, k, l, ll;
    double big, dum, pivinv, swap;
    
    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for(j = 1; j <= n; j++)
    {
        ipiv[j]=0;
    }
    
    for(i = 1; i <= n; i++)
    {
        big = 0.0;
        for(j = 1; j <= n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (k = 1;k <= n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if(fabs(a[j][k]) >= big)
                        {
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
            
        if(irow != icol)
        {
            for(l = 1; l <= n; l++)
            {
                SWAP(a[irow][l], a[icol][l])
            }
            for(l = 1; l <= m; l++)
            {
                SWAP(b[irow][l], b[icol][l])
            }
        }
        
        indxr[i] = irow;
        indxc[i] = icol;
        
        if(a[icol][icol] == 0.0)
        {
            //nrerror("gaussj: Singular Matrix");
            return -1;
        }
        
        pivinv = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        
        for(l = 1; l <= n; l++)
        {
            a[icol][l] *= pivinv;
        }
        
        for(l = 1; l <= m; l++)
        {
            b[icol][l] *= pivinv;
        }
        
        for(ll = 1; ll <= n; ll++)
        {
            if (ll != icol)
            {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for(l = 1; l <= n; l++)
                {
                    a[ll][l] -= a[icol][l] * dum;
                }
                for(l = 1; l <= m; l++)
                {
                    b[ll][l] -= b[icol][l] * dum;
                }
            }
        }
    }
    
    for(l = n; l >= 1; l--)
    {
        if(indxr[l] != indxc[l])
        {
            for(k = 1;k <= n; k++)
            {
                SWAP(a[k][indxr[l]], a[k][indxc[l]])
            }
        }
    }
    
    free_ivector(ipiv, 1, n);
    free_ivector(indxr, 1, n);
    free_ivector(indxc, 1, n);
    
    return 0;
}

void covsrt(double **covar, int ma, int ia[], int mfit)
{
    int i, j, k;
    double swap;
    for (i = mfit + 1; i <= ma; i++)
    {
        for (j = 1; j <= i; j++)
        {
            covar[i][j] = covar[j][i] = 0.0;
        }
    }
    
    k = mfit;
    
    for(j = ma; j >= 1; j--)
    {
        if(ia[j])
        {
            for(i = 1; i <= ma; i++)
            {
                SWAP(covar[i][k], covar[i][j])
            }
            for (i = 1; i <= ma; i++)
            {
                SWAP(covar[k][i], covar[j][i])
            }
            k--;
        }
    }
}

//void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[], int ma, double **alpha, double beta[], double *chisq, void (*funcs)(double, double [], double *, double [], int))
void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[], int ma, double **alpha, double beta[], double *chisq, double fconst)
{
    int i, j, k, l, m, mfit = 0;
    double ymod, wt, sig2i, dy, *dyda;
    
    dyda = dvector(1, ma);
    
    for(j = 1; j <= ma; j++)
    {
        if(ia[j])
        {
            mfit++;
        }
    }
    
    for(j = 1; j <= mfit; j++)
    {
        for(k = 1; k <= j; k++)
        {
            alpha[j][k] = 0.0;
        }
        beta[j] = 0.0;
    }
    
    *chisq = 0.0;
    for(i = 1; i <= ndata; i++)
    {
        (*funcs)(x[i], a, &ymod, dyda, ma, fconst);
        
        sig2i = 1.0/(sig[i] * sig[i]);
        dy = y[i] - ymod;

        for(j = 0, l = 1; l <= ma; l++)
        {
            if(ia[l])
            {
                wt = dyda[l] * sig2i;
                for(j++, k = 0, m = 1; m <= l; m++)
                {
                    if(ia[m])
                    {
                        alpha[j][++k] += wt * dyda[m];
                    }
                    beta[j] += dy * wt;
                }
            }
        }
        *chisq += dy * dy * sig2i;
    }
    
    for(j = 2; j <= mfit; j++)
    {
        for(k = 1; k < j; k++)
        {
            alpha[k][j] = alpha[j][k];
        }
    }
    free_dvector(dyda,1,ma);
}

//void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[], int ma, double **covar, double **alpha, double *chisq, void (*funcs)(double, double [], double *, double [], int), double *alamda)
int mrqmin(double x[], double y[], double sig[], int ndata, double *a, int ia[], int ma, double **covar, double **alpha, double *chisq, double *alamda, double fconst)
{
    int j, k, l;
    static int mfit;
    static double ochisq, *atry, *beta, *da, **oneda;
    int error;
    
    if (*alamda < 0.0)
    {
        atry = dvector(1, ma);
        beta = dvector(1, ma);
        da = dvector(1, ma);
        for(mfit = 0, j = 1; j <=ma ; j++)
        {
            if(ia[j])
            {
                mfit++;
            }
        }
        
        oneda = dmatrix(1, mfit, 1, 1);
        *alamda = 0.001;
        
        mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, fconst);
        
        ochisq = (*chisq);
        
        for(j = 1; j <= ma; j++)
        {
            atry[j] = a[j];
        }
    }
    
    for(j = 1; j <= mfit; j++)
    {
        for(k = 1; k <= mfit; k++)
        {
            covar[j][k] = alpha[j][k];
        }
        covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
        oneda[j][1] = beta[j];
    }
    
    error = gaussj(covar, mfit, oneda, 1);
    if(error == -1)
    {
        return -1;
    }
    
    for(j = 1;j <= mfit; j++)
    {
        da[j] = oneda[j][1];
    }
    
    if(*alamda == 0.0)
    {
        covsrt(covar, ma, ia, mfit);
        covsrt(alpha, ma, ia, mfit);
        
        free_dmatrix(oneda, 1, mfit, 1, 1);
        free_dvector(da, 1, ma);
        free_dvector(beta, 1, ma);
        free_dvector(atry, 1, ma);
        return 0;
    }
    
    for (j = 0, l = 1; l <= ma; l++)
    {
        if(ia[l])
        {
            atry[l] = a[l] + da[++j];
        }
    }

    mrqcof(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, fconst);

    if(*chisq < ochisq)
    {
        *alamda *= 0.1;
        ochisq = (*chisq);
        for(j = 1; j <= mfit; j++)
        {
            for(k = 1; k <= mfit; k++)
            {
                alpha[j][k]=covar[j][k];
            }
            beta[j] = da[j];
        }
        for(l = 1; l <= ma; l++)
        {
            a[l] = atry[l];
        }
    }
    else
    {
        *alamda *= 10.0;
        *chisq = ochisq;
    }
    
    return 0;
}
