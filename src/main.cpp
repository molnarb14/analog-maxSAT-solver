#include <iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../includes/s-a-variables.h"
#include "../includes/levenberg-marquardt.h"

using namespace std;

int N, M, K_max; //init
double alpha, tmax, A, B;
int Nprobb = 50000;
int Nprobb_initial;
int Nprobb_min = 1000;
int Nprobb_max = 200000;
int seed_max = 1000;
double tbin = 1; //0.01
double offset = 0.1; //offsetting energy levels at LM fit
int stat_threshold = 10000;
double p_exp_threshold = 0.000005;

string *svector(long nl, long nh)
{
    string *v;
    v = new string[nh + 1];
    if(!v)
    {
        nrerror("allocation failure in svector()");
    }
    return v;
}

void free_svector(string *v, long nl, long nh)
{
    delete[] v;
}

double *re_dvector(double *v, long nl, long nh)
{
    v = (double *)realloc(v, (size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if(!v)
    {
        nrerror("re-allocation failure in re_dvector()");
    }
    return v - nl + NR_END;
}

int *re_ivector(int *v, long nl, long nh)
{
    v = (int *)realloc(v, (size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
    if(!v)
    {
        nrerror("re-allocation failure in re_ivector()");
    }
    return v - nl + NR_END;
}

double **re_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

    m = (double **)realloc(m, (size_t)((nrow + NR_END) * sizeof(double*)));
    if(!m)
    {
        nrerror("reallocation failure 1 in matrix()");
    }

    m += NR_END;
    m -= nrl;

    m[nrl] = (double *)realloc(m[nrl], (size_t)((nrow * ncol + NR_END) * sizeof(double)));
    if(!m[nrl])
    {
        nrerror("reallocation failure 2 in matrix()");
    }

    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i = nrl + 1; i <= nrh; i++)
    {
        m[i] = m[i - 1] + ncol;
    }

    return m;
}

void checkFileOpening(FILE *file, char *path)
{
    if(!file)
    {
        printf("Error opening file %s...\n", path);
        exit(1);
    }
}

int findMax(int *vector)
{
    int max = vector[0];
    for(int i = 1; i <= M; i++)
    {
        if(vector[i] > max)
        {
            max = vector[i];
        }
    }
    return max;
}

void printkSAT2file(ksat k, int **matrix)	//print out the k-SAT instance encoded in matrix var
{
    char *read_test = (char *)malloc(256 * sizeof(char));
    strcpy(read_test, "read_test.dat");
    FILE *f = fopen(read_test, "w");
    checkFileOpening(f, read_test);
    fprintf(f, "N=%d M=%d\n", N, M);
    for(int i = 1; i <= M; i++)
    {
        for(int kk = 1; kk <= k.K[i]; kk++)
        {
            fprintf(f, "%d\t", matrix[i][kk]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    free(read_test);
}

int searchForOverallMin(int *vector)
{
    int min = vector[1];
    for(int i = 2; i <= Nprobb; i++)
    {
        if(vector[i] < min)
        {
            min = vector[i];
        }
    }
    return min;
}

void moveContainerElementsForward(int *container, int size)
{
    for(int i = 1; i < size; i++)
    {
        container[i] = container[i + 1];
    }
}

bool verifyChange(int *container, int size)
{
    bool change = false;
    int c = container[1];
    for(int i = 2; i <= size; i++)
    {
        if(container[i] != c)
        {
            change = true;
        }
    }
    return change;
}

void listContainer(int *container, int size)
{
    for(int i = 1; i <= size; i++)
    {
        printf("%d\t", container[i]);
    }
    printf("\n");
}

int pow10(int n)
{
    static int ten[10] = {
        1, 10, 100, 1000, 10000,
        100000, 1000000, 10000000, 100000000, 1000000000
    };

    return ten[n];
}

void createFolders(char **parts, char *prefix, int index, int elements, char *resultFolder, struct stat sb)
{
    if(index == elements)
    {
        return;
    }
    else
    {
        sprintf(resultFolder, "%s/%s", prefix, parts[index]);
        int dir_err = mkdir(resultFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (stat(resultFolder, &sb) != 0 && S_ISDIR(sb.st_mode))
        {
            if(-1 == dir_err)
            {
                printf("Error creating directory %s!\n", resultFolder);
                exit(1);
            }
        }

        createFolders(parts, resultFolder, index + 1, elements, resultFolder, sb);
    }
}

void createFolders(char *folderName, struct stat sb)
{
    int dir_err = mkdir(folderName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (stat(folderName, &sb) != 0 && S_ISDIR(sb.st_mode))
    {
        if(-1 == dir_err)
        {
            printf("Error creating directory %s!\n", folderName);
            exit(1);
        }
    }
}

void replaceCharacterInString(char *string, char source, char target)
{
    for(int i = 0; i < strlen(string); i++)
    {
        if(string[i] == source)
        {
            string[i] = target;
        }
    }
}

int main(int argc, char * const argv[])
{
    double smin = 0, smax = 1.0;
    int prob, prob_start, solved = 0, unsolved = 0, steps;
    double t1, t2, dt, dtalltime;
    double t_start, t_stop, t_run;
    double avgtime, avtime, totalTime;
    double deltat;
    char *filename = (char *)malloc(256 * sizeof(char));
    char *CFGfilename = (char *)malloc(256 * sizeof(char));
    char *filepath = (char *)malloc(256 * sizeof(char));
    char *logged = (char *)malloc(256 * sizeof(char));
    char *auxfilepath = (char *)malloc(256 * sizeof(char));
    int i, goodsteps, badsteps;
    double *ystart;
    double cycleThreshold = 0.0001;
    double timeCycleThreshold = 0.025;
    int probs;
    FILE *f;
    FILE *benchfile, *results;
    double **minE;
    double **koltsegek;
    int **p;
    int *p_emin_initial;
    int binindex = 0;
    int *bins;
    int *lastIndex;

    double hmin = 0.0000000000001;
    double h1 = 0.01;
    double maxerror = 0.001;
    int smallWindowSize = N * 20;

    int dummy;
    char dummy_c;
    bool stop = false;
    bool random = false;

    double *kappa, *p_nprobb;

    int ndata;
    int ma = 2; //no. of relevant parameters for LM: 2 for LM-2x, 3 for LM-orig, 3 for LM-6c
    double *x, *y, *sig, *a, **covar, **al_pha;
    double chisq = 0, alamdba = -1;
    int *ia;
    int fitted, reached;
    double a1, a2; //fitted parameters to be used at the end, calculating the results

    int lowest_e = M, noOfLowest = 0;
    int lowest_nprobb = 1;
    int first_e = M;
    double eta = 0;
    double predicted = 0;
    double f_constant = 1;
    double min_min;

    int *lowest_e_container, *fitted_container;
    int container_i = 1;
    int container_size = 5;
    bool believe = false;
    int believe_nprobb = Nprobb;
    bool believe_message = false, mind_change_message = false;
    int final;
    double p2exp;

    int *init, *seed_emin, seeds = 1;
    double *seed_time;
    string *benchmrks;
    char *benchmrksDummy = (char *)malloc(256 * sizeof(char));
    int going_under = 5;

    double fconst_good, chisq_good;

    int stop_code = 0;

    char *solutionFile = (char *)malloc(256 * sizeof(char));
    char *solutionsFolder = (char *)malloc(256 * sizeof(char));
    char *correspondingEnergyLevel = (char *)malloc(50 * sizeof(char));
    int iteration;
    char *iterationSuffix = (char *)malloc(50 * sizeof(char));
    FILE *solution;
    bool solutionExport = false;
    int error; //variable for checking ERROR: gaussj: Singular Matrix from Levenberg-Marquardt algorigm mrqmin/gaussj

    if(argc != 4)
    {
        printf("-------------------------------------------------------------------------\n");
		printf("AnalogSAT: Simulate a continuous-time dynamical system (CTDS) that solves\n");
		printf("a given SAT problem, or minimizes the number of violated clauses, after\n");
		printf("running a set of trajectories predicts the possible minimum number of\n");
		printf("unsatisfied clauses using the Levenberg-Marquardt non-linear curve\n");
		printf("fitting method. It also predicts the number of trajectories needed to\n");
		printf("reach the predicted minimum.\n");
		printf("-------------------------------------------------------------------------\n\n");
		printf("Usage %s <CFG filename> <tmax> <Nprobb>\n\n", argv[0]);
        printf("  CFG filename - config file filename without extension from cfg folder\n");
		printf("      CFG file must include the names of the problems in DIMACS cnf format\n");
		printf("      Use generateCFG.sh in the sh folder to generate CFG files for the\n");
		printf("      CNF files stored in the cnf folder.\n\n");
		printf("  tmax - maximum analog time to simulate\n");
		printf("      0 < tmax <= 150 (recommended 35 or 50)\n\n");
        printf("  Nprobb - initial number of trajectories\n");
		printf("      1000 <= Nprobb <= 2000000 (default: 50000)\n\n");
        printf("  ex.\n");
        printf("      %s 2016.ms_random.highgirth.3sat 50 100000\n", argv[0]);
        exit(1);
    }

    strcpy(filename, argv[1]);
    strcpy(CFGfilename, argv[1]);
    tmax = atoi(argv[2]);
    Nprobb = atoi(argv[3]);

    bool read = true;

    benchmrks = svector(1, seed_max);

    random = false;

    strcpy(filepath, "cfg/");
    strcat(filepath, filename);
    strcat(filepath, ".txt");

    f = fopen(filepath, "r");
    checkFileOpening(f, filepath);

    read = true;

    while(read)
    {
        if(!feof(f))
        {
            dummy = fscanf(f, "%s", benchmrksDummy);
            benchmrks[seeds].append(benchmrksDummy);
            seeds++;
        }
        else
        {
            read = false;
        }
    }

    char *pch;
    int elements = 0;
    char *today = (char *)malloc(256 * sizeof(char));
    int folderLayers = 10; //supposing no folders deeper then 10 layers
    char **folderPathParts = (char **)malloc(folderLayers * sizeof(char *));
    for(int i = 0; i < folderLayers; i++)
    {
        folderPathParts[i] = (char *)malloc(256 * sizeof(char));
    }
    strcpy(auxfilepath, filename);
    replaceCharacterInString(auxfilepath, '.', '/');

    pch = strtok(filename, ".");
    while(pch != NULL)
    {
        strcpy(folderPathParts[elements], pch);
        elements++;
        pch = strtok(NULL, ".");
    }

    char *resultFolder = (char *)malloc(256 * sizeof(char));
    char *logFolder = (char *)malloc(256 * sizeof(char));
    struct stat sb;
    char *prefix = (char *)malloc(256 * sizeof(char));
    strcpy(prefix, "results");
    createFolders(prefix, sb);
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    if(tm.tm_mon + 1 < 10)
    {
        if(tm.tm_mday < 10)
        {
            sprintf(today, "%d0%d0%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
        }
        else
        {
            sprintf(today, "%d0%d%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
        }
    }
    else
    {
        if(tm.tm_mday < 10)
        {
            sprintf(today, "%d%d0%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
        }
        else
        {
            sprintf(today, "%d%d%d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday);
        }
    }
    strcpy(folderPathParts[elements], today);
    elements++;
    //creating folders for results with corresponding hierarchy
    createFolders(folderPathParts, prefix, 0, elements, resultFolder, sb);

    bins = ivector(0, tmax);
    lastIndex = ivector(0, Nprobb_max); //needs to be reallocated if neccessary
    a = dvector(1, ma);
    ia = ivector(1, ma);
    covar = dmatrix(1, ma, 1, ma);
    al_pha = dmatrix(1, ma, 1, ma);

    lowest_e_container = ivector(1, container_size);
    fitted_container = ivector(1, container_size);

    seeds -= 1; //calculates +2 on read..., if only 1 file it is +1

    for(int mi = 1; mi <= seeds; mi++)
    {
        error = 0;

        strcpy(filepath, "cnf/");
        strcat(filepath, auxfilepath);
        strcat(filepath, "/");
        strcat(filepath, benchmrks[mi].c_str());
        strcat(filepath, ".cnf");

        benchfile = fopen(filepath, "r");
        checkFileOpening(benchfile, filepath);

        dummy = fscanf(benchfile, "%s%s%d%d", &dummy_c, &dummy_c, &N, &M);

        alpha = M * 1.0/N;

        fclose(benchfile);

        ksat k(N, M, alpha);

        i = 1;
        int kk = 1;
        int value;
        benchfile = fopen(filepath, "r");
        checkFileOpening(benchfile, filepath);

        dummy = fscanf(benchfile, "%s%s%d%d", &dummy_c, &dummy_c, &N, &M);

        while(i <= M)
        {
            dummy = fscanf(benchfile, "%d", &value);
            if(value != 0)
            {
                k.var[i][kk] = value;
                kk++;
            }
            else
            {
                k.K[i] = kk - 1;
                i++;
                kk = 1;
            }
            //if(feof(f)) //TODO: examine the need for this condition
            //{
            //    break;
            //}
        }
        fclose(benchfile);

        K_max = findMax(k.K);

        sanda sa(K_max, k.K, N, alpha);

        probs = Nprobb;
        avtime = 0;
        avgtime = 0;
        solved = unsolved = 0;

        kappa = dvector(0, M);
        p_nprobb = dvector(0, M);

        minE = dmatrix(0, M, 1, Nprobb_max); //needs to be reallocated if neccessary
        p = imatrix(0, (int)(tmax/tbin) + 1, 0, M);
        koltsegek = dmatrix(0, (int)(tmax/tbin) + 1, 0, M);
        p_emin_initial = ivector(0, M);

        x = dvector(1, M + 1);
        y = dvector(1, M + 1);
        sig = dvector(1, M + 1);

        ystart = dvector(1, k.N + k.M);

        for(i = 1; i <= 5; i++)
        {
            lowest_e_container[i] = fitted_container[i] = -1 * i;
        }

        for(i = 0; i <= Nprobb; i++)
        {
            for(int j = 0; j <= M; j++)
            {
                minE[j][i] = -1;
            }
            minE[M][i] = 0;
            lastIndex[i] = 0;
        }

        for(i = 0; i <= M; i++)
        {
            kappa[i] = 0;
            p_emin_initial[i] = 0;
        }

        int interval_u = -1, interval_l = -1;

        for(i = 0; i <= (int)(tmax/tbin) + 1; i++)
        {
            for(int j = 0; j <= M; j++)
            {
                p[i][j] = 0;
                koltsegek[i][j] = 0;
            }
        }

        for(i = 1; i <= M + 1; i++)
        {
            x[i] = y[i] = 0;
            sig[i] = 0.9999;
        }

        a[1] = 1.0;
        a[2] = 1.0;

        fitted = reached = 0;

        for(i = 1; i <= ma; i++)
        {
            ia[i] = 1;
        }

        strcpy(prefix, "logs");
        createFolders(prefix, sb);
        createFolders(folderPathParts, prefix, 0, elements, logFolder, sb);
        if(random)
        {
            sprintf(logged, "%s/k%dN%da%.2lf-%d.log", logFolder, K_max, N, alpha, init[mi]);
        }
        else
        {
            sprintf(logged, "%s/%s.log", logFolder, benchmrks[mi].c_str());
        }

        f = fopen(logged, "w");
        checkFileOpening(f, logged);

        prob_start = 1;

        lowest_e = M;
        lowest_nprobb = 0;
        noOfLowest = 0;
        min_min = M;
        believe_nprobb = 0;
        final = M;

        eta = 0;
        ndata = 0;
        p2exp = 0;
        prob_start = 1;
        stop = false;
        believe = false;
        mind_change_message = false;
        believe_message = false;

        for(i = 1; i <= M; i++)
        {
            p_emin_initial[i] = 0;
        }

        strcpy(filepath, resultFolder);
        strcat(filepath, "/");
        if(random)
        {
            fprintf(f, "Initializing random maxSAT instance no. %d...\n", init[mi]);
            fflush(f);
            k.initksat(init[mi]);
            strcat(filepath, CFGfilename);
        }
        else
        {
            fprintf(f, "Initializing %s maxSAT instance...\n", benchmrks[mi].c_str());
            fflush(f);
            strcat(filepath, benchmrks[mi].c_str());
        }
        strcat(filepath, ".dat");
        results = fopen(filepath, "w");
        checkFileOpening(results, filepath);

        t_start = clock();

        do
        {
            for(prob = prob_start; prob <= Nprobb; prob++)
            {
                if(prob == 2)
                {
                    sa.b = lowest_e * 1.0/M - pow(0.5, 2 * K_max);
                    if(sa.b < 0)
                    {
                        sa.b = 0;
                    }
                    fprintf(f, "b = %lf\n", sa.b);
                }

                if(prob % 1000 == 0)
                {
                    fprintf(f, "%d\n", prob);
                    fflush(f);
                }
                totalTime = 0;
                binindex = 0;

                t1=clock();

                init_random(prob);
                sa.init_srand(&k, prob, smin, smax);

                for (i = 1; i <= k.N; i++)  ystart[i] = sa.s[i];
                for (i = 1; i <= k.M; i++)  ystart[N + i] = sa.a[i];
                steps = sa.rungekutta_maxt(ystart, k.N + k.M, 0.0, 100000.0, maxerror, h1, hmin, &goodsteps, &badsteps, &deltat, tmax, A, B, smallWindowSize, cycleThreshold, timeCycleThreshold, tbin, minE, prob, binindex, lastIndex);
                t2 = clock();
                dt = (t2 - t1) * 1E3/CLOCKS_PER_SEC;     //running time in miliseconds
                dtalltime = (t2 - t_start) * 1E3/CLOCKS_PER_SEC;

                if(prob == 1)
                {
                    first_e = lastIndex[prob];
                }

                //jumps to fit if lowest_e changes
                //if(prob > 100 && lowest_e > lastIndex[prob])
                //{
                //    Nprobb = prob;
                //    break;
                //}
                if(prob > 100)
                {
                    if(p[binindex][lowest_e] >= 100)
                    {
                        if(lowest_e > lastIndex[prob])
                        {
                            Nprobb = prob;
                            break;
                        }
                    }
                    else
                    {
                        if(lowest_e >= lastIndex[prob])
                        {
                             Nprobb = prob;
                             break;
                        }
                    }
                }

                //if(lowest_e > lastIndex[prob])
                //{
                //    lowest_nprobb = prob;
                //    lowest_e = lastIndex[prob];
                //}

                strcpy(prefix, "solutions");
                createFolders(prefix, sb);
                createFolders(folderPathParts, prefix, 0, elements, solutionsFolder, sb);
                if(lowest_e > lastIndex[prob])
                {
                    lowest_nprobb = prob;
                    lowest_e = lastIndex[prob];

                    //exporting SAT solution for a new energy level found
                    iteration = 0;
                    sprintf(correspondingEnergyLevel, "-E%d", lastIndex[prob]);
                    strcpy(solutionFile, solutionsFolder);
                    strcat(solutionFile, "/");
                    if(random)
                    {
                        strcat(solutionFile, CFGfilename);
                    }
                    else
                    {
                        strcat(solutionFile, benchmrks[mi].c_str());
                    }
                    strcat(solutionFile, correspondingEnergyLevel);
                    strcat(solutionFile, ".dat");
                    solutionExport = true;
                }
                else
                {
                    if(lowest_e == lastIndex[prob])
                    {
                        while(access(solutionFile, F_OK) != -1)
                        {
                            iteration++;
                            sprintf(iterationSuffix, "-[%d]", iteration);
                            strcpy(solutionFile, solutionsFolder);
                            strcat(solutionFile, "/");
                            if(random)
                            {
                                strcat(solutionFile, CFGfilename);
                            }
                            else
                            {
                                strcat(solutionFile, benchmrks[mi].c_str());
                            }
                            strcat(solutionFile, correspondingEnergyLevel);
                            strcat(solutionFile, iterationSuffix);
                            strcat(solutionFile, ".dat");
                        }
                        solutionExport = true;
                    }
                }
                if(solutionExport)
                {
                    solution = fopen(solutionFile, "w");
                    checkFileOpening(solution, solutionFile);
                    sa.printtofile_sanda(solution);
                    fclose(solution);
                }
            }

            for(i = prob_start; i <= Nprobb; i++)
            {
                for(int j = M; j >= lastIndex[i]; j--)
                {
                    if(minE[j][i] != -1 && minE[j - 1][i] == -1)
                    {
                        interval_u = j;
                    }
                    if(minE[j][i] == -1 && minE[j - 1][i] != -1)
                    {
                        interval_l = j;
                    }
                    if(interval_l != -1 && interval_u != -1)
                    {
                        for(int jj = interval_u - 1; jj >= interval_l; jj--)
                        {
                            minE[jj][i] = minE[interval_l - 1][i];
                        }
                        interval_l = -1;
                        interval_u = -1;
                    }
                }
            }

            for(int j = prob_start; j <= Nprobb; j++)
            {
                if(minE[lastIndex[j]][j] > 0) //>=1 //originally only >
                {
                    for(int m = M; m >= lastIndex[j]; m--)
                    {
                        for(int binindex = (int)(minE[lastIndex[j]][j]/tbin); binindex <= (int)(tmax/tbin); binindex++)
                        {
                            p[binindex][m]++;
                        }
                    }
                }
            }

            int bin_check = (int)(tmax/tbin);
            double t_check = bin_check * tbin;

            for(i = M; i >= 0; i--)
            {
                if(p[bin_check][i] > 0 && p[bin_check][i] < Nprobb)
                {
                    kappa[i] = -1.0 * log(1 - p[bin_check][i] * 1.0/Nprobb)/t_check;
                    p_nprobb[i] = p[bin_check][i] * 1.0/Nprobb;
                }
            }

            int ii = 1; //because x, y index start from 1

            for(i = M; i >= 0; i--)
            {
                if(kappa[i] != 0)
                {
                    x[ii] = kappa[i];
                    y[ii] = i;
                    lowest_e = i;
                    ii++;
                }
            }

            if(container_i <= container_size)
            {
                lowest_e_container[container_i] = lowest_e;
                //container_i++; //There is no need to increment this variable here, because in paralell we should register the corresponding fitted values and we will increment the index there
            }
            else
            {
                moveContainerElementsForward(lowest_e_container, container_size);
                lowest_e_container[container_size] = lowest_e;
            }

            ndata = ii - 1; //counting no. of variables should start from 0

            if(ndata == 1)
            {
                //stop condition no. 1
                min_min = lowest_e;
                believe_nprobb = Nprobb;
                stop = true;
                stop_code = 1;
                break;
            }

            eta = -1.0 * log(kappa[lowest_e])/log(N);
            fprintf(f, "eta = %lf\n", eta);

            fprintf(f, "%d\t\t%lf\t\t%lf\n", Nprobb, dt, dtalltime);
            fflush(f);

            double chisq_prev = 0, chisq_prev_prev = 0;
            double a1_prev, a1_prev_prev, a2_prev, a2_prev_prev;
            bool chisq_prev_set = false, chisq_prev_prev_set = false;
            bool fit = true;
            int l = 1; //because we start one level over the minimum energy level reached

            f_constant = lowest_e + 2 * offset; //it has a factor of 2, because as a first step inside the upcoming loop we decrement it with the offset

            do
            {
                alamdba = -1;
                a[1] = 1.0;
                a[2] = 1.0;

                f_constant -= offset;
                if(f_constant < -1)
                {
                    //stop condition no. 2
                    fit = false;
                    stop_code = 2;
                    break;
                }

                for(i = 1; i <= 100 && error != -1; i++)
                {
                    error = mrqmin(x, y, sig, ndata, a, ia, ma, covar, al_pha, &chisq, &alamdba, f_constant);
                }

                if(error == -1)
                {
                    fprintf(f, "gaussj: Singular Matrix\n");
                    fit = false;
                }
                else
                {
                    fprintf(f, "---Fit no. %d---\n", l);
                    l++;
                    fprintf(f, "const = %lf a1 = %lf a2 = %lf chi2 = %lf lambda = %lf | %d [%d]\n", f_constant, a[1], a[2], chisq, alamdba, p[bin_check][lowest_e], lowest_e);
                    fflush(f);

                    if(chisq_prev_set && chisq_prev_prev_set && chisq > chisq_prev && chisq_prev > chisq_prev_prev)
                    {
                        fconst_good = f_constant + 2 * offset;
                        a1 = a1_prev_prev;
                        a2 = a2_prev_prev;
                        chisq_good = chisq_prev_prev;
                        predicted = fconst_good;
                        min_min = fconst_good;
                        stop_code = 3;
                        fit = false;
                    }
                    else
                    {
                        if(chisq_prev_set)
                        {
                            chisq_prev_prev = chisq_prev;
                            a1_prev_prev = a1_prev;
                            a2_prev_prev = a2_prev;
                            chisq_prev_set = false;
                        }
                    }

                    if(!chisq_prev_prev_set)
                    {
                        chisq_prev_prev = chisq;
                        a1_prev_prev = a[1];
                        a2_prev_prev = a[2];
                        chisq_prev_prev_set = true;
                    }
                    else
                    {
                        if(!chisq_prev_set)
                        {
                            chisq_prev = chisq;
                            a1_prev = a[1];
                            a2_prev = a[2];
                            chisq_prev_set = true;
                        }
                    }
                }
            }
            while(fit);

            if(error != -1)
            {
                double p_exp = 0;
                double power_base = (lowest_e - 1 - fconst_good) * 1.0/a1; //it contained fconst_good - offset for some misteriously 'unknown' reason(s), but it has no logical point
                if(power_base <= 0)
                {
                    p_exp = 0; //Nprobb;
                }
                else
                {
                    double power_pow = 1.0/a2;
                    double power = pow(power_base, power_pow);
                    double exponential = exp(-1.0 * t_check * power);
                    p_exp = Nprobb * (1 - exponential);
                }
                predicted = fconst_good;
                min_min = fconst_good;
                p2exp = p_exp;
                noOfLowest = p[bin_check][lowest_e];

                fprintf(f, "---Final desicion about the fit---\n");
                fprintf(f, "const = %lf a1 = %lf a2 = %lf chi2 = %lf | %d [%d]\n", fconst_good, a1, a2, chisq_good, p[bin_check][lowest_e], lowest_e);
                fprintf(f, "lowest_e = %d, believed_lowest_e = %lf, p_exp = %lf, ndata = %d\n", lowest_e, fconst_good, p_exp, ndata);
                fflush(f);

                //How we STOP?
                //1
                if(container_i <= container_size)
                {
                    fitted_container[container_i] = (int)(min_min) + 1;
                    container_i++; //this is the only place where we should increment the index so the lowest_e values recorded remain paired with the corresponding fitted values
                }
                else
                {
                    moveContainerElementsForward(fitted_container, container_size);
                    fitted_container[5] = (int)(min_min) + 1;
                }

                bool changedE = verifyChange(lowest_e_container, container_size); //record change in lowest_e
                bool changedF = verifyChange(fitted_container, container_size); //record change in fitted value

                if(Nprobb >= stat_threshold) //1.a & 1.b
                {
                    if(!changedF && !changedE) //1.c //1.d.ii: || (ndata < 5 && p[bin_check][lowest_e] >= 100)
                    {
                        //We start to believe :)
                        believe = true;
                        mind_change_message = false;
                        if(!believe_message)
                        {
                            fprintf(f, "I start to believe... :)\n");
                            fflush(f);
                            believe_message = true;
                            believe_nprobb = Nprobb;
                        }
                        else
                        {
                            fprintf(f, "I still believe... :)\n");
                            fflush(f);
                            mind_change_message = false;
                        }

                    }
                    else
                    {
                        believe = false;
                        believe_nprobb = Nprobb;
                        believe_message = false;
                        if(!mind_change_message && believe_message)
                        {
                            fprintf(f, "Actually I just changed my mind... Please stand by...\n");
                            fflush(f);
                            mind_change_message = true;
                        }
                        else
                        {
                            fprintf(f, "Still deciding... Please stand by...\n");
                            fflush(f);
                        }
                    }
                }

                int p_cutoff = 10;
                if(eta >= 3)
                {
                    p_cutoff = 5;
                }

                //2
                if(believe) //2.b
                {
                    if(lowest_e == 0 || (p[bin_check][lowest_e] >= p_cutoff && p_exp >= 100) || Nprobb >= Nprobb_max) //2.b.i && 2.b.iii && 2.b.v //2.b.iv: || p_exp > 100
                    {
                        stop_code = 4;
                        stop = true;
                    }
                    else //2.b.ii
                    {
                        stop = false;
                    }
                    //if(ndata <= 5 && p[bin_check][lowest_e] >= 1000)
                    if(p[bin_check][lowest_e] >= 1000)
                    {
                        fprintf(f, "Problem too easy... Assessing results...\n");
                        fflush(f);
                        stop_code = 5;
                        stop = true;
                    }
                    if(lowest_e - 1 < min_min && Nprobb > stat_threshold)
                    {
                        if(going_under != 0)
                        {
                            going_under--;
                        }
                        else
                        {
                            fprintf(f, "We reached the lowest possible energy level %d, having an assymptote in %.2lf...\n", lowest_e, min_min);
                            fflush(f);
                            stop_code = 6;
                            stop = true;
                        }
                    }
                }

                if(Nprobb >= Nprobb_max)
                {
                    stop_code = 7;
                    stop = true;
                }

                if(!believe || !stop) //2.a && 2.b.ii'
                {
                    //if((ndata <= 5 && p[bin_check][lowest_e] >= 1000) || (p[bin_check][lowest_e] >= 100 && p_exp >= 100) || (p_exp <= p_exp_threshold && (int)(min_min + 1) == lowest_e && Nprobb > stat_threshold))
                    if((ndata <= 5 && p[bin_check][lowest_e] >= 1000) || (p[bin_check][lowest_e] >= 100 && p_exp >= 100))
                    {
                        stop_code = 8;
                        stop = true;
                        min_min = lowest_e;
                        believe_nprobb = Nprobb;

                        if(Nprobb < Nprobb_min)
                        {
                            stop = false;
                        }
                        else
                        {
                            break;
                        }
                    }

                    prob_start = Nprobb + 1;

                    if(p_exp > 0)
                    {
                        if(p_exp == Nprobb)
                        {
                            Nprobb *= 1.20;
                        }
                        else
                        {
                            int jump = (int)(Nprobb * 1.0/p_exp);
                            if(jump - Nprobb < 100)
                            {
                                Nprobb += 100;
                            }
                            else
                            {
                                Nprobb += (int)(Nprobb * 1.0/p_exp);
                            }
                        }
                    }
                    else
                    {
                        Nprobb *= 1.20;
                    }

                    if(Nprobb > Nprobb_max)
                    {
                        Nprobb = Nprobb_max;
                    }

                    if(Nprobb < 0)
                    {
                        Nprobb = Nprobb_max;
                    }

                    //reallocating minE, lastIndex
                    //minE = re_dmatrix(minE, 0, M, 1, Nprobb);
                    //lastIndex = re_ivector(lastIndex, 0, Nprobb);

                    fprintf(f, "Nprobb_tng = %d\n", Nprobb);
                    fflush(f);
                }
            }
            else
            {
                stop = true;
            }
        }
        while(!stop);

        t_stop = clock();
        t_run = (t_stop - t_start) * 1E3/CLOCKS_PER_SEC;     //running time in miliseconds

        fprintf(f, "\nFinished after running: %lf\n", t_run);

        if(Nprobb >= stat_threshold)
        {
            if(lowest_e == first_e || (ndata < 5 && noOfLowest >= 100))
            {
                fprintf(f, "Finished after running %d instances...\nProblem too easy...\nLowest possible energy: %d\nExit code: %d\n", Nprobb, lowest_e, stop_code);
                final = lowest_e;
            }
            else
            {
                if(min_min >= lowest_e)
                {
                    predicted = (int)min_min;
                    fprintf(f, "Finished after running %d instances...\nOptimal expected: %d\nOptimal found: %d\nLowest energy found at: %d\nStarted to believe: %d\nExporting results, and exiting...\nExit code: %d\n", Nprobb, (int)(predicted), lowest_e, lowest_nprobb, believe_nprobb, stop_code);
                }
                else
                {
                    predicted = (int)(min_min + 1);
                    fprintf(f, "Finished after running %d instances...\nOptimal expected: %d\nOptimal found: %d\nLowest energy found at: %d\nStarted to believe: %d\nExporting results, and exiting...\nExit code: %d\n", Nprobb, (int)(predicted), lowest_e, lowest_nprobb, believe_nprobb, stop_code);
                }
                final = (int)predicted;
            }
        }
        else
        {
            predicted = (int)min_min;
            fprintf(f, "Finished after running %d instances...\nExited with smaller statistics than the threshold %d...\nOptimal expected: %d\nFitted value: %lf\nLowest possible energy: %d\nExit code: %d\n", Nprobb, stat_threshold, (int)predicted, fconst_good, lowest_e, stop_code);
            final = lowest_e;
        }
        fflush(f);
        fclose(f);

        char *kappa_out = (char *)malloc(256 * sizeof(char));
        char *kappaFolder = (char *)malloc(256 * sizeof(char));
        strcpy(prefix, "kappa");
        createFolders(prefix, sb);
        createFolders(folderPathParts, prefix, 0, elements, kappaFolder, sb);

        if(random)
        {
            sprintf(kappa_out, "%s/k%dN%da%.2lf-%d.dat", kappaFolder, K_max, N, alpha, init[mi]);
        }
        else
        {
            sprintf(kappa_out, "%s/%s.dat", kappaFolder, benchmrks[mi].c_str());
        }

        f = fopen(kappa_out, "w");
        checkFileOpening(f, kappa_out);

        for(i = M; i >= 0; i--)
        {
            if(kappa[i])
            {
                fprintf(f, "%d\t%.12lf\t%.12lf\n", i, kappa[i], p_nprobb[i]);
            }
        }
        fclose(f);
        free(kappa_out);
        free(kappaFolder);

        fprintf(results, "%d\t%s\t%d\t%d\t%d\t%d\t%.2lf\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%d\n", mi, benchmrks[mi].c_str(), lowest_e, lowest_nprobb, noOfLowest, (int)(fconst_good + 1), fconst_good, believe_nprobb, final, Nprobb, eta, ndata, p2exp, t_run, sa.b, stop_code);
        fflush(results);

        free_dvector(ystart, 1, k.N + k.M);

        free_dmatrix(minE, 0, M, 1, Nprobb_max);
        free_imatrix(p, 0, (int)(tmax/tbin) + 1, 0, M);
        free_dmatrix(koltsegek, 0, (int)(tmax/tbin) + 1, 0, M);
        free_ivector(p_emin_initial, 0, M);

        free_dvector(x, 1, M + 1);
        free_dvector(y, 1, M + 1);
        free_dvector(sig, 1, M + 1);

        free_dvector(kappa, 0, M);
        free_dvector(p_nprobb, 0, M);

        //sanda, ksat freed on exiting class
    }

    fclose(results);

    free_svector(benchmrks, 1, seed_max);

    for(int i = 0; i < folderLayers; i++)
    {
        free(folderPathParts[i]);
    }
    free(folderPathParts);
    free(today);
    free(resultFolder);
    free(filename);
    free(CFGfilename);
    free(filepath);
    free(auxfilepath);
    free(logged);
    free(benchmrksDummy);
    free(logFolder);
    free(solutionsFolder);
    free(solutionFile);
    free(iterationSuffix);
    free(correspondingEnergyLevel);
    free(prefix);

    free_ivector(bins, 0, tmax);
    free_ivector(lastIndex, 0, Nprobb_max);
    free_dvector(a, 1, ma);
    free_ivector(ia, 1, ma);
    free_dmatrix(covar, 1, ma, 1, ma);
    free_dmatrix(al_pha, 1, ma, 1, ma);

    free_ivector(lowest_e_container, 1, container_size);
    free_ivector(fitted_container, 1, container_size);

    return 0;
}
