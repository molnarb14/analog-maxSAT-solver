////functions:
//	init_random(unsigned int the_seed)
//	double rand_real(double r_min,double r_max)
//	int rand_int(int r_min,int r_max)
#include <math.h>

#define PI 3.141592654
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum_ACT=-18276459;

float ran2(long *idum)
// call it with a negative seed
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7; j>=0; j--)
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 +=IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

void init_random(long the_seed)
{
    idum_ACT=-18276459+the_seed;
}

void init_random_time()
{
    idum_ACT=-18276459+time(0);
}


double rand_real(double r_min,double r_max)
{
    return (ran2(&idum_ACT)*(r_max-r_min)+r_min);
}

long long rand_int(long long r_min,long long r_max)
{
    return ((long long)(ran2(&idum_ACT)*(r_max+1-r_min)+r_min));
}

/*void randomize(double x[],unsigned long N){
	double *x_old;
	short int *gotya;
	unsigned long i,j,n;

	x_old=new double[N+1];
	gotya=new short int[N+1];

	for(i=1;i<=N;i++) {
		x_old[i]=x[i];
		gotya[i]=0;
	}

	for(n=0;n<N;n++){
		j=rand_int(1,N-n);
		x[n+1]=x_old[j];
		x_old[j]=x_old[N-n];
		x_old[N-n]=x[n+1];
	}
	delete[] x_old;x_old=NULL;
}*/

float expdev()
{
    // p(result)=
//Returns an exponentially distributed, positive, random deviate of unit mean,
// using ran2(&idum_ACT) as thesource of uniformdeviates. ititialize with
//init_randrom(negative seed)
    float dum;
    do dum=ran2(&idum_ACT);
    while (dum == 0.0);
    return -log(dum);
}

float gasdev()
{
//Returns anormally distributed deviate with zero mean and unit variance,
// using ran2(&idum_ACT) as the source of uniformdeviates.
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
//if (*idum < 0) iset=0; //Reinitialize.
    if (iset == 0)   //We don’t have an extra deviate handy, so
    {
        do
        {
            v1=2.0*ran2(&idum_ACT)-1.0; //pick two uniform numbers in the square extending from -1 to +1 in each direction,
            v2=2.0*ran2(&idum_ACT)-1.0;
            rsq=v1*v1+v2*v2; //see if they are in the unit circle,
        }
        while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
        fac=sqrt(-2.0*log(rsq)/rsq);      //Now make the Box-Muller transformation to get two normal deviates.
        //Return one and save the other for next time.
        gset=v1*fac;
        iset=1;			  //Set flag.
        return v2*fac;
    }
    else  						//We have an extra deviate handy,
    {
        iset=0;					//so unset the flag,
        return gset;				//and return it.
    }
}

double rand_Gauss(double mean_value, double szigma_mean)
{
    return(gasdev()*szigma_mean+mean_value);
}
/*
float Poisson_random(float xm){
//Returns as af loating-point number an integer value that is a random deviate
//drawn from a Poisson distribution of mean xm, using ran2(&idum_ACT)
//as a source of uniform random deviates.

 static float sq,alxm,g,oldm=(-1.0); //oldm is a flag for whether xm has changed
									 //since last call.
 float em,t,y;
 if (xm < 12.0) {						//Usedirect method.
	 if (xm != oldm) {
		oldm=xm;
		g=exp(-xm);					// If xm is new, compute the exponential.
	 }
	 em = -1;
	 t=1.0;
	 do {				//Instead of adding exponential deviates it is equivalent to
						//multiply uniform deviates. Wenever actually have to take the
						//log, merely compare to the pre-computed exponential.
		 ++em;
		 t *= ran2(&idum_ACT);
		 }
	 while (t > g);
	 }
	 else {				//Use rejection method.
		if (xm != oldm) { 		//If xm has changed since the last call,
								//then precompute some functions that occur below.
				oldm=xm;
				sq=sqrt(2.0*xm);
				alxm=log(xm);
				g=xm*alxm-gammln(xm+1.0); //The function gammln is the natural log of
											//the gamma function, asgivenin§6.1.
				}
		do {
			do {			//y is a deviate from a Lorentzian comparison function.
				y=tan(PI*ran2(&idum_ACT));
				em=sq*y+xm; //em is y, shifted and scaled.
			}
			while (em < 0.0); //Reject if in regime of zero probability.
			em=floor(em);	//The trick for integer-valued distributions.
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
							//The ratio of the desired distribution to the
							//comparison function; we accept or reject by comparing it
							// to another uniform deviate. The factor 0.9 is chosen
							//so that t never exceeds 1.
		}
		while (ran2(&idum_ACT) > t);
	}
	return em;
}
*/
