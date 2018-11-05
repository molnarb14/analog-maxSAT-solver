#include "rand_fg.h"
#include <math.h>
#include "nrutil-g++.cpp"

class ksat{

	public:
        unsigned int Kmax, N, M; //, *K_gen, K_max;	//number of variables (N) and number of clause (M) in a K-SAT instance
        int *K;
		double alpha; //ratio alpha = M/N
		int **var; //the matrix encoding the SAT instance. Each row gives the variables of a given clause: with positive sign if it's in normal form, and negative sign if it's negated


        ksat(unsigned int nn, unsigned int mm, double alf)
        {
            N = nn;
            M = mm;
            K = ivector(1, M);
            var = imatrix(1, M, 1, Kmax);
        }
    
		ksat(unsigned int kmax, int *kk, unsigned int nn, double alf)
        {
			Kmax = kmax; //number of variables in a clause   K-SAT
            K = kk;
			N = nn;	//number of variables
			alpha = alf; //ratio alpha=M/N
            M = (int)(alpha * N); //number of clauses
			var = imatrix(1, M, 1, Kmax); //the matrix encoding the SAT instance. Each row gives the variables of a given clause: with positive sign if it's in normal form, and negative sign if it's negated
        }
    
//        ksat(unsigned int *kk, unsigned int nn, unsigned int mm, int kmax)
//        {
//            K_gen = kk;					//number of variables in a clause   K-SAT
//            K_max = kmax;
//            N = nn;					//number of variables
//            M = mm;                 //number of clauses
//            alpha = M * 1.0/N;				//ratio alpha=M/N
////            M = (int)(alpha * N);		//number of clauses
//            var = imatrix(1, M, 1, K_max);   //the matrix encoding the SAT instance. Each row gives the variables of a given clause: with positive sign if it's in normal form, and negative sign if it's negated
//        }


		//initializing a random K-SAT instance based on a given random-seed. The seeds we use were previously randomly selected numbers, which we solved with another program (minisat)
		//and separated them in two different lists: satisfiable and unsatisfiable. This is how we can check when the solution is not found for SATISFIABLE instances.
		void initksat(int rseed)
        {
            int i, j, m, p, already_in_clause;
			init_random(rseed); //initialize random generator with the given seed
			for(i = 1; i <= M; i++)
            {
				for(j = 1; j <= Kmax; j++)
                {
                    var[i][j] = 0;
                }
            }
            
			for(m = 1; m <= M; m++)
            {
				for(i = 1; i <= Kmax; i++)
                {
                    p = (int)(rand_int(1, N));				//randomly choose a variable
					already_in_clause = 0;
					for(j = 1; j <= i - 1; j++)
                    {
                        if(abs(var[m][j]) == p) already_in_clause=1;     //check if it's already in the clause
                    }
                    
					if (already_in_clause == 0)					//if it's not in the clause yet, choose between normal and negated form, which is encoded as the sign of var[m][i]
                    {
                        var[m][i] = (int)((2 * rand_int(0, 1) - 1) * p);
                    }
					else
                    {
                        i--;
                    }
                }
            }
        }

		void printksat()	//print out the k-SAT instance encoded in matrix var
        {
            int i, k;
			printf("N=%d M=%d\n", N, M);
			for(i = 1; i <= M; i++)
			{
//			    printf("%d: ", i);
			    for(k = 1; k <= K[i]; k++)
			    {
			        printf("%d\t", var[i][k]);
			    }
			    printf("\n");
            }
        }

        void readfromfile(char nev[])	//print to file the k-SAT instance encoded in matrix var
        {
            int m, k;
            int dummy;
            FILE *f;
            f = fopen(nev, "r");
            for(m = 1; m <= M; m++)
            {
                for(k = 1; k <= K[m]; k++)
                {
                    if(k == 3)
                    {
                        dummy = fscanf(f, "%d%d", &var[m][k], &dummy);  //az m-edik mondat k-adik valtozoja
                    }
                    else
                    {
                        dummy = fscanf(f, "%d", &var[m][k]);  //az m-edik mondat k-adik valtozoja
                    }
                }
            }
            fclose(f);
        }
    
		void printtofile(char nev[])	//print to file the k-SAT instance encoded in matrix var
			{int i,m,k;
			FILE *f;
			f=fopen(nev,"w");
			fprintf(f,"p cnf %d %d\n",N,M);
			for(m=1;m<=M;m++)
				{for(k=1;k<=K[m];k++)
					{i=var[m][k];  //az m-edik mondat k-adik valtozoja
					fprintf(f,"%d ",i);
					}
				fprintf(f,"0\n");
				}
			fclose(f);
			}


		int checksolution_binary(int ss[])   //check if a given binary (0,1) array ss is a solution
			{int m,j,i,correct,sbin;
			//for(i=1;i<=N;i++)  printf("%d\t",ss[i]);
			//printf("\n");
			for (m=1;m<=M;m++)
				{correct=0;
				for (j=1;j<=K[m];j++)
					{i=abs(var[m][j]);			//var[m][j] can be positive or negative, depending on the normal or negated form, so first we get the index of the variable
					if (ss[i]==0)  sbin=-1;
					if (ss[i]==1) sbin=1;
					if (sbin*var[m][j]>0) correct++;    //it is TRUE
					}
				if (correct==0) return(0);	 //this clause is FALSE => it is not a solution
				}
			return(1);
			}

		int checksolution_real(double ss[])   //check if a given real array ss (with values between -1 and 1) is a solution
			{int m,j,i,correct;
			//for(i=1;i<=N;i++)  printf("%lf\t",ss[i]);
			//printf("\n");//getchar();
			for (m=1;m<=M;m++)
				{correct=0;
				for (j=1;j<=K[m];j++)
					{i=abs(var[m][j]);
					if (ss[i]*var[m][j]>0.0) correct++;
					//	printf("%d\t%lf\t%d\n",var[m][j],ss[i],correct);
					}
				//printf("\n");
				//getchar();
				if (correct==0) return(0);	//this clause is FALSE => it is not a solution
				}
			return(1);
			}

/// the following functions were only used when generating the two-dimensional maps shown in Fig4. they are not part of the algorithm solving the k-SAT instance

		//encodes the solotion array as one number (we use it when we do the two-dimensional maps on Fig.4 and we need to save the solution obtained in each point of the map)
		long long int solutionnumber(double ss[])
			{int m,j,i,sbin;
			//for(i=1;i<=N;i++)  printf("%lf\t",ss[i]);
			//printf("\n");//getchar();
			long long int sol_number=0;
			for (i=1;i<=N;i++)
				{if (ss[i]>=0.0)  sbin=1;
				else sbin=0;
				sol_number=sol_number*2+sbin;
				}
			return(sol_number);
			}

		//given a solution number gives back the s array
		void inverse_solutionnumber(long long int sol_number, double ss[])
			{int m,j,i;
			//printf("sol_number=%lld\n",sol_number);
			for (i=N;i>=1;i--)
				{
				ss[i]=(sol_number%2)*2.0-1.0;
				sol_number=sol_number/2;
				}
			//for(i=1;i<=N;i++)  printf("%lf\t",ss[i]);
			//printf("\n");getchar();
			}

		//we want to find out if the two solutions are in the same solutioncluster. At very small alpha or higher N, this can take very long time, but we used this only for Figure 4, it is not a function normally used when solving a k-SAT instance
		int checkcluster(long long int solution1, long long int solution2)
				{int *X;
				long long int **sol,thissolution;
				int layer,m,i,jj,kk,already_included;
				double *s;
				//printf("solution%lld ,solution2=%lld\n",solution1,solution2);getchar();
				X=ivector(0,N);
				sol=llmatrix(0,N,1,N*N*N*N*N);   //this stores the solutions in different layers (layer 1 are solutions being one-flip away from solution1, in layer 2 they are two-flips away and so on)
				s=dvector(1,N);
				for(i=1;i<=N;i++)  X[i]=0;			//number of solutions in a given layer
				X[0]=1;
				sol[0][1]=solution1;
				for(layer=1;layer<=N;layer++)      //layer 1 are solutions being one-flip away from solution1, in layer 2 they are two-flips away and so on
					{//printf("layer=%d\n",layer);//getchar();
					for(m=1;m<=X[layer-1];m++)		//we build the next layer, looking at the neighbors of solutions in the previous layer
						{//printf("m=%d\n",m);
						inverse_solutionnumber(sol[layer-1][m],s);	//take the solution number and obtain the s array itself
						for(i=1;i<=N;i++)
							{s[i]=-s[i];					//check if flipping s[i] still satisfies the SAT instance
							if (checksolution_real(s)==1)		//if yes
								{thissolution=solutionnumber(s);
								//than check if it's already included in previous layers
								already_included=0;
									for(jj=1;jj<=layer;jj++)
										for(kk=1;kk<=X[jj];kk++)  if (thissolution==sol[jj][kk]) already_included=1;
								if (already_included==0)			// if it's not included
									{X[layer]++;					// increase the number of solutions in this layer and store it in the sol array
									sol[layer][X[layer]]=thissolution;
									//printf("Thissolution=%lld\tsolution2=%lld\n",thissolution,solution2);
									if (thissolution==solution2)	//if this is the solution2,  what we are looking for the function returns 1
										{free_ivector(X,0,N);
										free_llmatrix(sol,0,N,1,N*N*N*N*N);
										free_dvector(s,1,N);
										//printf("found X[%d]=%d\n",layer,X[layer]);getchar();
										return(1);
										}
									}
								}
							s[i]=-s[i];	    //we flip it back so we can check the other "one-flip neighbors"
							}
						}
					//printf("X[%d]=%d\n",layer,X[layer]);getchar();
					if (X[layer]==0)		//we found the whole solution cluster and solution2 is not part of it
						{free_ivector(X,0,N);
						free_llmatrix(sol,0,N,1,N*N*N*N*N);
						free_dvector(s,1,N);
						return(0);
						}
					}

                    return 0;
                }

		//the function used for generating the maps in Fig4
		void twodmap(char *filename,char *filename2,int grid1,int grid2,double tmax)
			{int i,j,m,colors,p,satisfied,steps,*color,included,*cluster;
			 int clusters,solutions;
			long long int sol_number,*lista;
			double tt,dt,s1,s2,*ss;
			FILE *f,*g;
                int dummy;
			lista=llvector(1,100000);
			for(i=1;i<=100000;i++)  lista[i]=-1;

			color=ivector(1,100000);
			cluster=ivector(1,10000);

			f=fopen(filename,"r");					//this file contains for each of hte points on the grid the result given by the algorihtm
			g=fopen(filename2,"w");

			colors=1;				//colors are needed only when making figures, it has no algorithmic importance (the fact that it start from 1 is only for my convinience when using xmgrace to generate the figures)
			solutions=0;			//this is the number of different solutions found so far
			clusters=0;				//this is the number of solution-clusters found so far
			for(i=-grid1;i<=grid1;i++)				//we go through each point of the grid
				for(j=-grid2;j<=grid2;j++)
					{//printf("%d\t%d\n",i,j);
					dummy = fscanf(f,"%lf%lf%lld%d%d%lf%lf",&s1,&s2,&sol_number,&satisfied,&steps,&tt,&dt);   //tt is the analog time,   dt is the running time of the code solving the given instance, we don't use it
					//printf("sol_number=%lld\n",sol_number);
					included=0;
					for(m=1;m<=solutions;m++)			//check if this solution is alerady included in our list.
						{
						if (sol_number==lista[m])	//if yes, than we already know its cluster and color
							{
							//fprintf(g,"%lf\t%lf\t%d\t%d\n",s1,s2,color[m]%50,2+(int)(tt/tmax*90));
							fprintf(g,"%lf\t%lf\t%d\t%d\t%d\t%lf\t%d\n",s1,s2,color[m],color[m]%50,cluster[m]+1,tt,2+(int)(tt/tmax*90));  //this might not seem to be too logical, it only helped me to prepare the figures easier in xmgrace (for example, color 1 is white in xmgrace, that's why I add 1, etc...)
							//printf("%lf\t%lf\t%d\t%d\t%lf\t%d\n",s1,s2,color[m],cluster[m]+1,tt,2+(int)(tt/tmax*90));
							included=1;
							}
						}
					if (included==0)		//if not, we increase the number of solutions and we include it in our list and give a new color
							{solutions++;
							//printf("solutions=%d\n",solutions);getchar();
							lista[solutions] =sol_number;
							colors++;
							color[solutions]=colors;
							//printf("solutions=%d\tclusters=%d\n",solutions,clusters);	//getchar();
							p=0;
							for(m=1;m<=solutions-1;m++)    //we check if it's in the same cluster with any of the previous solutions
								{
								if (p==0) {p=checkcluster(lista[m],sol_number);
											if (p==1)  {cluster[solutions]=cluster[m];}
										  }
								}
							if (p==0)  {clusters++;cluster[solutions]=clusters;}  //this is in a new cluster

							//fprintf(g,"%lf\t%lf\t%d\t%d\n",s1,s2,color[solutions]%50,2+(int)(tt/tmax*90));
							fprintf(g,"%lf\t%lf\t%d\t%d\t%d\t%lf\t%d\n",s1,s2,color[solutions],color[solutions]%50,cluster[solutions]+1,tt,2+(int)(tt/tmax*90));
							//printf("%lf\t%lf\t%d\t%d\t%lf\t%d\n",s1,s2,color[solutions],cluster[solutions]+1,tt,2+(int)(tt/tmax*90));

							}


					}

			fclose(g);
			fclose(f);

			/*ss=dvector(1,N);
			for(i=1;i<=solutions;i++)
				{printf("megoldas=%lld\tcluster=%d\n",lista[i],cluster[i]);
					inverse_solutionnumber(lista[i],ss);
					for(j=1;j<=N;j++)  printf("%.0lf ",ss[j]);
					printf("\n");
				}
			*/

			//printf("solutions=%d\tclusters=%d\n",solutions,clusters);
			free_ivector(color,1,100000);
			free_llvector(lista,1,100000);
			free_ivector(cluster,1,10000);
		}



		~ksat()
			{
                free_imatrix(var, 1, M, 1, Kmax);
                if(K != NULL)
                {
                    free_ivector(K, 1, M);
                }
			}

};

