#include "maxSAT.h"
#include <math.h>
//#include <nrutil.h>
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

class sanda
{
public:
    double *s, *a, *cl, **cli;
    ksat *kSAT;
    int *K, Kmax, N, M, step, problem;
    double t;
    double A, B, tmax;
    int sc;
    double b;

    sanda(int kmax, int *kk, int nn, double alpha)
    {
        N = nn;					//number of variables
        M = (int)(nn * alpha);		//number of clauses
        K = kk;					// K-SAT
        Kmax = kmax;
        s = dvector(1, N);		//s variables
        a = dvector(1, M);		//auxiliary  variables
        cl = dvector(1, M);		//cl[m] is the value of clause m, the equivalent of K_m in the paper
        cli = dmatrix(1, M, 1, Kmax);
        b = 0.05; //0.5; //alpha;
    }

    void init(ksat *sat, int prob)			//initializing the s and a variables, with s=0
    {
        int i, j;
        kSAT = sat;
        problem = prob;
        for(i = 1; i <= N; i++) s[i] = 0.0;
        for (i = 1; i <= M; i++) a[i] = 1.0;
        for(i = 1; i <= M; i++)
        {
            cl[i] = 0.0;
            for(j = 1; j <= K[i]; j++)
            {
                cli[i][j] = 0.0;
            }
        }
    }


    //initializing variables s and a, with |s[i]| having   random values between smin and smax
    //normally smin=0 and smax=1, we just used this for special measurements, which are not included in this paper
    void init_srand(ksat *sat, int prob, double smin, double smax)
    {
        int i, j;
        kSAT = sat;
        problem = prob;
        for(i = 1; i <= N; i++)
        {
            s[i] = (rand_int(0, 1) * 2 - 1); //* rand_real(smin, smax);
            //printf("%lf\t",s[i]);
        }

        for(i = 1; i <= M; i++)
        {
            a[i] = 1.0; //rand_real(0, 1);
        }

        for(i = 1; i <= M; i++)
        {
            cl[i] = 0.0;
            for(j = 1; j <= K[i]; j++)
            {
                cli[i][j] = 0.0;
            }
        }

    }

    //initializing variables s and a, with s[i] having   random values in a small box placed in a random point
    //again this was not used for any of the results presented in our paper
    void init_srand_randombox(ksat *sat, int prob, double *boxmiddle, double box)   //if ptof=1 then print to file
    {
        int i, j;
        kSAT = sat;

        problem = prob;
        step = 0;
        t = 0.0;


        for(i = 1; i <= N; i++)
        {
            s[i] = boxmiddle[i] + rand_real(-box, box);
            //printf("%lf\t",s[i]);
        }

        for (i = 1; i <= M; i++) a[i] = 1.0;

        for(i = 1; i <= M; i++)
        {
            cl[i] = 0.0;
            for(j = 1; j <= K[i]; j++)
            {
                cli[i][j] = 0.0;
            }
        }
    }


    //calculates the clauses
//    void clauses() //CNN specific
//    {
//        int i, j, m, l;
//        double *minus, *plus;
//        double most, pow2k, pow2km1;
//        minus = dvector(1, N);
//        plus = dvector(1, N);
//        minus[0] = 0;
//        plus[0] = 0;
//        for(i = 1; i <= N; i++)				//calculating these ahead, saves computation
//        {
//            minus[i] = 1 - function(s[i]);
//            plus[i] = 1 + function(s[i]);
//        }
//
//        pow2k = 1.0/pow(2, K);			//normalizing factor in K_m (see paper)
//        pow2km1=1.0/pow(2,K-1);
//
//        for(m = 1; m <= M; m++)
//        {
//            cl[m] = pow2k;
//            for (j = 1; j <= K; j++)
//            {
//                cli[m][j] = pow2km1;
//            }
//            for (j = 1; j <= K; j++)
//            {
//                i = abs(kSAT->var[m][j]);
//                if(kSAT->var[m][j] > 0)
//                {
//                    most = minus[i];
//                }
//                else
//                {
//                    most = plus[i];
//                }
//                cl[m] *= most;
//                for(l = 1; l <= K; l++)
//                {
//                    if(i != abs(kSAT->var[m][l])) cli[m][l] *= most;
//                }
//            }
//        }
//
//        free_dvector(plus, 1, N);
//        free_dvector(minus, 1, N);
//    }

    void clauses() //ORIG specific
    {
        int i, j, l, m;
        double *minus, *plus;
        double most, *pow2k, *pow2km1;
        pow2k = dvector(1, M);
        pow2km1 = dvector(1, M);
        minus = dvector(1, N);
        plus = dvector(1, N);

        minus[0] = 0;
        plus[0] = 0;

        for(i = 1; i <= N; i++)				//calculating these ahead, saves computation
        {
            minus[i] = 1 - s[i];
            plus[i] = 1 + s[i];
        }

        for(int i = 1; i <= M; i++)
        {
            pow2k[i] = 1.0/pow(2, K[i]);			//normalizing factor in K_m (see paper)
            pow2km1[i] = 1.0/pow(2, K[i]);
        }

        for(m = 1; m <= M; m++)
        {
            cl[m] = pow2k[m]; //printf("%lf ",cl[m]);
            for(j = 1; j <= K[m]; j++)
            {
                cli[m][j] = pow2km1[m];
            }
            for(j = 1; j <= K[m]; j++)
            {
                i = abs(kSAT->var[m][j]);
                if(kSAT->var[m][j] > 0)
                {
                    most = minus[i];
                }
                else
                {
                    most = plus[i];
                }
                cl[m] *= most;
                for(l = 1; l <= K[m]; l++)
                {
                    if(i != abs(kSAT->var[m][l]))
                    {
                        cli[m][l] *= most;
                    }
                }
            }
        }

        free_dvector(pow2k, 1, M);
        free_dvector(pow2km1, 1, M);
        free_dvector(plus, 1, N);
        free_dvector(minus, 1, N);
    }

    //some useful functions for monitoring the behavior of the dynamics
    double energy_V()
    {
        double v = 0, avgam = 0; //0
        int m, i;
        clauses();
        for(m = 1; m <= M; m++)
        {
            v += a[m] * cl[m] * cl[m]; //v+=gunction(a[m])*cl[m]*cl[m];
            avgam += a[m];
        }
        avgam /= N;

        for(i = 1; i <= N; i++)
        {
            v += avgam * b * cos(s[i] * M_PI/2.0) * cos(s[i] * M_PI/2.0);
        }
        //printf("V=%lf\n",v);
        return(v);
    }

    double energy_E()
    {
        double v=0;
        int i, m;
        clauses();
        for(m=1; m<=M; m++)  	v+=cl[m]*cl[m];
//        for(i = 1; i <= N; i++) v += b * cos(s[i] * M_PI/2.0) * cos(s[i] * M_PI/2.0);
        //printf("V=%lf\n",v);
        return(v);
    }

    double radius()
    {
        double R = 0, Rn = 0;
        int i;
        for(i = 1; i <= N; i++)
        {
            R += s[i] * s[i];//R+=function(s[i])*function(s[i]);
        }
        //printf("%lf\t%lf\n",R,Rn);getchar();
        return R;
    }

    double maximal_s()
    {
        double smax=0;
        int i;
        for(i=1; i<=N; i++)  if (fabs(s[i])>smax) smax=fabs(s[i]);
        printf("%lf\n",smax);
        getchar();
        return(smax);
    }

    double minimal_s()
    {
        double smin=1;
        int i;
        for(i=1; i<=N; i++)  if (fabs(s[i])<smin) smin=fabs(s[i]);
        printf("%lf\n",smin);
        getchar();
        return(smin);
    }

    int countUnsatisfiedClauses(double ss[])   //check if a given real array ss (with values between -1 and 1) is a solution
    {
        int m, i, j, correct, satisfied = 0;

        for(m = 1; m <= M; m++)
        {
            correct = 0;
            for(j = 1; j <= K[m]; j++)
            {
                i = abs(kSAT->var[m][j]);
                if(ss[i] * kSAT->var[m][j] > 0.0)
                {
                    correct++;
                }
            }
            if(correct >= 1) satisfied++;
        }
        return M - satisfied;
    }


    //equations of CNN differential equations
    double function(double x)
    {
        return 0.5 * (fabs(x + 1) - fabs(x - 1));
    }

    double gunction(double x)
    {
        return 0.5 * (fabs(x) - fabs(x - 1) + 1);
    }

    double noise()
    {
        double ns;
        double min = -0.005;
        double max = 0.005;
        ns = rand_real(min, max);
        //printf("*debugInfo: %lf\n", ns);//getchar();
        return ns;
    }

    //calculating the derivatives for the RK
//    void derivs(double x,double y[],double dydx[]) //CNN specific
//    {
//        int i,m,sign,j;
//        for(i=1; i<=N; i++)  s[i]=y[i];			//array y includes both s and a variables
//        for(m=1; m<=M; m++)  a[m]=y[N+m];
//
//        //printf("*debugInfo: starting derivs calcualtion\n");
//
//        for(i=1; i<=N; i++)
//        {
//            dydx[i] = -1.0 * s[i] + A * function(s[i]);// + noise(); //!!!!!!!!!!!!!!!!
//        }
//        //printf("*debugInfo: dxdy[i] processed\n");
//        for(m=1; m<=M; m++)
//        {
//            dydx[N+m] = -1.0 * a[m] + B * gunction(a[m]) + 1 - K;// + noise(); //!!!!!!!!!!!!!!!!
//        }
//        //printf("*debugInfo: dxdy[N+m] processed\n");
//
//        for(m=1; m<=M; m++)
//        {
//            for(j=1; j<=K; j++)
//            {
//                i=abs(kSAT->var[m][j]);
//
//                if (kSAT->var[m][j]>0)
//                {
//                    sign=1;
//                }
//                else
//                {
//                    sign=-1;
//                }
//                //printf("*debugInfo: sign defined\n");
//                dydx[i] += sign * gunction(a[m]) ; //!!!!!!!!!!!!!!!!!
//                //printf("*debugInfo: dxdy[%d] refreshed\n", i);
//                dydx[N+m] +=  -1.0 * sign * function(s[i]);    //!!!!!!!!!!!
//                //printf("*debugInfo: dxdy[N+%d] refreshed\n", m);
//            }
//        }
//
//        //printf("*debugInfo: derivs calculated\n");
//    }

    void derivs(double x, double y[], double dydx[]) //ORIG specific
    {
        int i, m, sign, j;
        double avgam = 0;
        for(i = 1; i <= N; i++)
        {
            s[i] = y[i];			//array y includes both s and a variables
        }
        for(m = 1; m <= M; m++)
        {
            a[m] = y[N + m];
            avgam += a[m];
        }
        avgam /= N;

        clauses();

        for(i = 1; i <= N; i++)
        {
            dydx[i] = avgam * b * sin(s[i] * M_PI) * M_PI/2; //0
        }

        for(m = 1; m <= M; m++)
        {
            for(j = 1; j <= K[m]; j++)
            {
                i = abs(kSAT->var[m][j]);
                if(kSAT->var[m][j] > 0) sign = 1;
                else sign = -1;
                dydx[i] += 2 * a[m] * sign * cl[m] * cli[m][j];
            }
        }

        for(m = 1; m <= M; m++)
        {
            dydx[N + m] = cl[m] * a[m];		//derivatives of the auxiliary variables
        }
    }


    void print_sanda()
    {
        int i;
        printf("%d\t%lf\n\nS:\t", step, t);
        for(i = 1; i <= N; i++) printf("%lf\t%lf\t ", s[i], function(s[i]));
        printf("\n\nA:\t");
        for(i = 1; i <= M; i++) printf("%lf\t%lf\t ", a[i], gunction(a[i]));
        printf("%lf\n", b);
    }

    void printtofile_ER()
    {
        double E0, E, R;
        char name[250];
        FILE *f;
//        sprintf(name,"../../results/N%d-tng2/pt/ER/P%d-ER-N%da%.2lf-A%.3lf-B%.3lf.dat", N, problem, N, (double)(M)/N, A, B);
        sprintf(name, "../../results/maxSAT/SP2/ER/P%d-ER-N%da%.2lf.dat", problem, N, (double)(M)/N);
        if (step==1) f=fopen(name,"w");
        else f=fopen(name,"a");
        E0 = energy_E();
        //E = energy_V();
        E = countUnsatisfiedClauses(s);
        R = radius();
        //fprintf(f,"%d\t%lf\t%lf\t%lf\t%lf\n",step,t,E0,E,R);
//        fprintf(f,"%d\t%lf\t%lf\n",step,t,R);
        fprintf(f, "%d\t%lf\t%lf\t%lf\t%lf\n", step, t, E0, E, R);
        fclose(f);
    }

    void printtofile_sanda()
    {
        int i;
        char name[150];//, name2[150];
        FILE *file;//, *gfile;
        sprintf(name,"../../results/N%d-tng2/pt/sanda/P%d-sanda-N%da%.2lf-tmax%lf-A%.3lf-B%.3lf.dat", N, problem,N, (double)(M)/N, tmax, A, B);
//        sprintf(name2,"../../results/pt/fandg/P%d-fandg-N%da%.2lf-tmax%lf-A%.3lf-B%.3lf.dat",problem,N,(double)(M)/N, tmax, A, B);

        if (step==1)
        {
            file=fopen(name,"w");
//            gfile = fopen(name2, "w");
        }
        else
        {
            file=fopen(name,"a");
//            gfile = fopen(name2, "a");
        }

        fprintf(file,"%d\t%lf\tS:",step,t);
        for(i=1; i<=N; i++) fprintf(file,"\t%lf",s[i]);
        fprintf(file,"\tA:");
        for(i=1; i<=M; i++) fprintf(file,"\t%lf",a[i]);
        fprintf(file,"\n");

//        fprintf(gfile,"%d\t%lf\tF(Si):",step,t);
//        for(i=1; i<=N; i++) fprintf(gfile,"\t%lf",function(s[i]));
//        fprintf(gfile,"\tG(Am):");
//        for(i=1; i<=M; i++) fprintf(gfile,"\t%lf",gunction(a[i]));
//        fprintf(gfile,"\n");

        fclose(file);
//        fclose(gfile);
    }

    void printtofile_sanda(FILE *file)
    {
        int i;

        fprintf(file, "%d\t%lf\tS:", step, t);
        for(i = 1; i <= N; i++) fprintf(file, "\t%lf", s[i]);
        fprintf(file, "\tA:");
        for(i = 1; i <= M; i++) fprintf(file, "\t%lf", a[i]);
        fprintf(file, "\n");
    }

    //cycle detection utils
    void MoveElementsForward(double *vector, int size)
    {
        for(int i = 1; i < size; i++)
        {
            vector[i - 1] = vector[i];
        }
    }

    double getMax(double *vector, int size)
    {
        double max = vector[0];
        for(int i = 1; i < size; i++)
        {
            if(vector[i] > max) max = vector[i];
        }
        return max;
    }

    double getMin(double *vector, int size)
    {
        double min = vector[0];
        for(int i = 1; i < size; i++)
        {
            if(vector[i] < min) min = vector[i];
        }
        return min;
    }

    void normalize(double *source, double *destination, int size)
    {
        for(int i = 0; i < size; i++)
        {
            destination[i] = source[i]/source[size - 1];
            //printf("%lf\t%lf\t%lf\n", source[i], source[size - 1], destination[i]); //getchar();
        }
    }

    void normalizeWithSum(double *source, double *destination, int size)
    {
        double sum_temp = 0;
        int i;
        for(i = 0; i < size; i++)
        {
            sum_temp += source[i];
        }
        sum_temp = sum_temp*1.0/(size - 1);
        for(i = 0; i < size; i++)
        {
            destination[i] = source[i]/sum_temp;
        }
    }

    double sum(double *vector, int size)
    {
        double s = 0;
        for(int i = 0; i < size; i++)
        {
            s += vector[i];
        }
        return s;
    }

    double sumWithDt(double *vector, double dT, int size)
    {
        double sum = 0;
        for(int i = 0; i < size; i++)
        {
            sum += ((vector[i] - dT)*1.0/dT)*((vector[i] - dT)*1.0/dT);
        }
        return sum;
    }

    double getAverage(double *vector, int size)
    {
        //printf("**debugInfo: size=%d\n", size); //getchar();
        double avg = 0;
        for(int i = 0; i < size; i++)
        {
            avg += vector[i];
        }
        return avg/size;
    }

    double averageDt(double *source, int size)
    {
        double sum_temp = 0;
        int i;
        for(i = 0; i < size; i++)
        {
            sum_temp += source[i];
        }
        return sum_temp*1.0/size;
    }

    double getAverage(double *vector, int size, bool print)
    {
        printf("**debugInfo: size=%d\n", size); //getchar();
        double avg = 0;
        for(int i = 0; i < size; i++)
        {
            avg += vector[i];
            if(print) printf("[%d] = %lf\t", i, vector[i]);
            if(i%50 == 0) getchar();
        }
        printf("\n%lf\t%lf\n", avg, avg/size); getchar();
        return avg/size;
    }
    //cycle detection utils - end

    //adaptive Runge-Kutta taken from numerical recipies
    void rkck(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[])
    {
        int i;
        static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2;
        static double b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2;
        static double b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0;
        static double b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0;
        static double b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0;
        static double c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0;
        static double dc5 = -277.00/14336.0;
        double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,dc6=c6-0.25;
        double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
        ak2 = dvector(1, n);
        ak3 = dvector(1, n);
        ak4 = dvector(1, n);
        ak5 = dvector(1, n);
        ak6 = dvector(1, n);
        ytemp = dvector(1, n);

        for(i=1; i<=n; i++) //first step
            ytemp[i]=y[i]+b21*h*dydx[i];
        derivs(x+a2*h,ytemp,ak2); //second step
        for (i=1; i<=n; i++)
            ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
        derivs(x+a3*h,ytemp,ak3); //Third step.
        for (i=1; i<=n; i++)
            ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
        derivs(x+a4*h,ytemp,ak4); //Fourth step.
        for (i=1; i<=n; i++)
            ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
        derivs(x+a5*h,ytemp,ak5); //Fifth step.
        for (i=1; i<=n; i++)
            ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
        derivs(x+a6*h,ytemp,ak6); //Sixth step.
        for (i=1; i<=n; i++) //Accumulate increments with proper weights.
            yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
        for (i=1; i<=n; i++)
            yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
        //Estimate error as difference between fourth and  fth order methods.
        free_dvector(ytemp, 1, n);
        free_dvector(ak6, 1, n);
        free_dvector(ak5, 1, n);
        free_dvector(ak4, 1, n);
        free_dvector(ak3, 1, n);
        free_dvector(ak2, 1, n);
    }

    int rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,double yscal[], double *hdid, double *hnext)
    {
        int i, j;
        double errmax, h, htemp, xnew, *yerr, *ytemp;
        yerr = dvector(1, n);
        ytemp = dvector(1, n);
        h = htry;
        for(;;)
        {
            rkck(y, dydx, n, *x, h, ytemp, yerr);
            errmax = 0.0;
            for(i = 1; i <= n; i++) errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
            errmax /= eps;
            
            if(errmax <= 1.0) break;
            htemp = SAFETY * h * pow(errmax, PSHRNK);
            //Truncation error too large, reduce stepsize.
            h = (h >= 0.0 ? FMAX(htemp, 0.1 * h) : FMIN(htemp, 0.1 * h));
            //No more than a factor of 10.
            xnew = (*x) + h;
            if(xnew == *x) return 0; //nrerror("stepsize underflow in rkqs");
        }
        if(errmax > ERRCON) *hnext = SAFETY * h * pow(errmax, PGROW);
        else *hnext = 5.0 * h; //No more than a factor of 5 increase.
        *x += (*hdid = h);
        for(i = 1; i <= n; i++) y[i] = ytemp[i];
        free_dvector(ytemp, 1, n);
        free_dvector(yerr, 1, n);
        return 1;
    }

    int rungekutta_maxt(double ystart[], int nvar, double x1, double x2, double eps_original, double h1, double hmin, int *nok, int *nbad, double *deltat, double MAXT, double aa, double bb, int smallWindowSize, double cycleThreshold, double timeCycleThreshold, double tbin, double **minE, int prob, int &binindex, int *lastIndex)
    {
        A = aa;
        B = bb;
        sc = 0;
        tmax = MAXT;

        int nstp, i, m;
        double x, hnext, hdid, h, eps;
        double *yscal, *y, *dydx;
        int sign_changed, nullavaltozas, underflow;
        double *s_previous;

        int En, En_min = M;
        double TINY;
        double sum_fsi, sum_gam, sum_si, sum_am, sum_derivs;
        TINY = 1.0e-20;
        s_previous = dvector(1, N);
        yscal = dvector(1, nvar);
        y = dvector(1, nvar);
        dydx = dvector(1, nvar);
        x = x1;
        t = x1;
        step = 0;
        h = SIGN(h1, x2 - x1);
        *nok = 0;
        *nbad = 0;
        double l, R0; //cycleThreshold

        double dT_small, dT_big, sum_max_small, sum_min_small, sum_max_big, sum_min_big;
        double r, r_old;
        bool cycle_small, cycle_big, cycle_time;

        for(i = 1; i <= nvar; i++) y[i] = ystart[i];
        eps = eps_original;
        //E0=energy_E();
        //E=energy_V();
        nullavaltozas = 0;
        t = 0;
        nstp = 0;
        l = 2;
        r = 0;
        r_old = 0;
        R0 = 0;
        sum_max_small = 0;
        sum_min_small = 0;
        sum_max_big = 0;
        sum_min_big = 0;

        cycle_small = false;
        cycle_big = false;
        cycle_time = false;

        dT_small = smallWindowSize;
        dT_big = 1.5 * smallWindowSize;

        int prev_binindex;
        int *bin_flag = ivector(0, (int)(MAXT/tbin) + MAXT);

        for(int i = 0; i < (int)(MAXT/tbin) + MAXT; i++)
        {
            bin_flag[i] = 0;
        }

        do
        {
            //printf("**inside Runge-Kutta\tt = %lf\tlastIndex[%d] = %d\n", t, prob, lastIndex[prob]); getchar();
            nstp++;
            prev_binindex = binindex;
            //printf("%lf\t%d\n",t,nstp);
            for(i = 1; i <= N; i++)
            {
                s_previous[i] = y[i];
            }
            derivs(x, y, dydx);   //also updates s and a

            //printtofile_sanda();//getchar();
            //print_sanda(); //getchar();
            //printtofile_ER();
            //printf("itt...");
            for(i = 1; i <= nvar; i++)	//Scaling used to monitor accuracy. This general-purpose choice can be modiif need be.
            {
                yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
            }

            if((x + h - x2) * (x + h - x1) > 0.0)
            {
                h = x2 - x; //If stepsize can overshoot, decrease.
            }

            underflow = rkqs(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext);
            if(underflow == 0)
            {
                return 0;
            }

            t = x;
            step = nstp;

            if (hdid == h)
            {
                ++(*nok);
            }
            else
            {
                ++(*nbad);
            }

            for(i=1; i<=N; i++)  s[i] = y[i];
            for(i=1; i<=M; i++)  a[i] = y[N + i];

            //There are different possible startegies for detecting if we found the solution
            //none of them really affects the analog time, and teh results obtained are almost perfectly identical
            //the first strategy we use here in this routine, is that we check for solution whenever one of the s variables changed its sign (we enter a different N-quadrant)
            //in the next routine we use an energy threshold

            //detect if any of the s variables changed its sign
            sign_changed = 0;
            sum_fsi = 0;
            sum_gam = 0;
            sum_si = 0;
            sum_am = 0;
            sum_derivs = 0;

            for(i = 1; i <= N; i++)
            {
                if(s[i] * s_previous[i] <= 0)
                {
                    sign_changed = 1;
                }
            }

            /*for(i = 1; i <= N; i++)
            {
                sum_fsi += function(s[i]) * function(s[i]);
                sum_si += s[i];
            }

            for(m = 1; m <= M; m++)
            {
                sum_gam += gunction(a[m]);
                sum_am += a[m];
            }

            for(i = 1; i <= N + M; i++)
            {
                sum_derivs += dydx[i] * dydx[i];
            }*/

            //Energy approximation for hard-SAT v. 2.0
            En = countUnsatisfiedClauses(s);

            if(En_min > En)
            {
                En_min = En;
                /*if(minE[En][prob] == -1)
                {*/
                    minE[En][prob] = t;
                    lastIndex[prob] = En;
                /*}
                else
                {
                    lastIndex[prob] = En_min;
                }*/
            }

            //end - energy approx.

            if(fabs(hnext) <= hmin)
            {
                return 0;/*nrerror("Step size too small in odeint");*/
            }
            h = hnext;
        }
        while(t < MAXT);

        free_dvector(s_previous, 1, N);
        free_dvector(yscal, 1, nvar);
        free_dvector(y, 1, nvar);
        free_dvector(dydx, 1, nvar);
        free_ivector(bin_flag, 0, (int)(MAXT/tbin) + MAXT);

        return(nstp);
    }

    ~sanda()
    {
        free_dvector(s, 1, N);
        free_dvector(a, 1, M);
        free_dvector(cl, 1, M);
        free_dmatrix(cli, 1, M, 1, Kmax);
    }
};
