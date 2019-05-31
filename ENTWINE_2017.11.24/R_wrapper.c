// These functions are wrapped by corresponding R functions using the .C interface, allowing computation via R.
// The R user supplies the path (directory), project (name of files), # of replicates, genotype (as an array of integers corresponding to gene IDs)
// The R wrapper functions supply the constants from hardcoded values and allocate memory for the return values, which are seen here as arguments.
// Each function has a common core which loads the genome from binary files.
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "G.h"
#include "develop.h"
#include "utilities.h"

const gsl_rng_type* rT;
gsl_rng* rg;

//Global variables used to id mutations, duplications, and novel haplotypes
int next;
int nextLocus;
int nextHaploid;

//Declarations for global constants
unsigned int sValue;
int num_T;
int num_G;
int num_P;
int num_SE;
int num_ST;
int num_SP;
int chromosomes;
int ploidy;
int sex;
double rec_rate;
int max_b;
double basic_signal_rate;
int k_mut;
int k_weak;
double mu_nuc;
int L_cis;
int n_cis;
double mu_ts;
double mu_cs;
double mu_small_K;
double mu_big_K;
double mu_del;
double mu_dup;
double mu_dbm;
double mu_trait;
double mu_signalK;
double mu_signalN;
double mu_whichSignal;
double trans_mutational_target;
double trait_mutational_target;
double dbm_mutational_target;
double cis_mutational_target;
double sig_mutational_target;
double sigma_ts;
double sigma_trait;
double sigma_cs;
double sigma_signalK;
double pos_cs;
double initial_mean_cs;
double initial_mean_ts;
double initial_sigma_trait;
double initial_mean_signalK;
double max_time;
double basal;
double transcription;
double translation;
double mRNA_decay;
double protein_decay;
double eps;
int BUFFER;
double omega_squared;
double* trait_optima;
double* trait_optima_env;
double* trait_optima_perf;
double environment;
double Env;
double Env_env;
double Env_perf;
int GeneToMutate;
//double high_environment;
int prefix;
double K_vals [4];
int dualPhenoGenes;
int treatment;
int EnvSignal;

// mutational constants

double mu_pheno;

// Geometric Mean Confidence Interval variables

double seLogFit;
double varLogFit;
double meanLogFit;

// Stuff
int nbNewMutationsNeeded;

// Just to fix a bug!
double pleiotropy_proba;

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double kth_smallest(double* a, int n, int k)
{
    register int i,j,l,m ;
    register double x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}


#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))


void MR (
    char** path,
    char** project,
    int* intConstants,
    double* doubleConstants,
    int* replicate,
    int* geno,
    int* nGenes,
    int* reps,
    int* nbmutants,
    int* nbmutations,
    int* nbStepsModel,
    double* meanW,
    double* medianW,
    double* sdW,
    double* environment,
    double* meanP,
    double* medianP,
    double* sdP,
    int* trace // Either 1 or 0. If 1, then the last mutant is used as the genome for the next nb mutations. 'nbmutations' can then likely be an array of 1.
    )
{
    // Mutational Constants
    n_cis = 12;
    L_cis = 500;
    sigma_ts = 0.5;
    sigma_trait = 0.5;
    sigma_cs = 0.5;
    pos_cs = 0.9;
    initial_mean_cs = 0.5;
    k_mut = 4;

    // Adjust Mutational constants
    mu_nuc = 0.001;
    mu_pheno = 100;

    mu_ts = mu_nuc * mu_pheno;
    mu_trait = mu_nuc * mu_pheno;
    mu_dbm = mu_nuc * mu_pheno;
    mu_cs = mu_nuc * 50;
    mu_small_K = mu_nuc * n_cis;
    mu_big_K = mu_nuc * k_mut * L_cis * factorial(n_cis) * pow(0.75, k_mut) * pow(0.25, n_cis-k_mut) / factorial(k_mut) / factorial(n_cis - k_mut);

    // Other constants
    ploidy = intConstants[0];
    max_b = intConstants[1];
    k_weak = intConstants[2];
    num_T = intConstants[3];
    chromosomes = intConstants[4];
    dualPhenoGenes = intConstants[5];
    sex = intConstants[6];
    EnvSignal = intConstants[7];
    
    max_time = doubleConstants[0];
    basal = doubleConstants[1];
    transcription = doubleConstants[2];
    translation = doubleConstants[3];
    mRNA_decay = doubleConstants[4];
    protein_decay = doubleConstants[5];
    eps = doubleConstants[6];
    omega_squared = doubleConstants[7];
    //environment = doubleConstants[8];

    setup();
    int offset = sizeof(int) * (8 + max_b + num_T) + sizeof(double) * (4 + max_b + num_T);
    //printf("%d\n", offset);
    char* geneFilename = calloc(500, sizeof(char));
    char geneDesig[] = "bin";
    char geneSuffix[] = ".bin";
    char underscore[] = "_";
    char slash[] = "/";
    char prefixChar [4];
    prefixChar[0] = (char)(replicate[0] / 100 + 48);
    prefixChar[1] = (char)((replicate[0] /10) % 10 + 48);
    prefixChar[2] = (char)(replicate[0] % 10 + 48);
    prefixChar[3] = '\0';
    strcat(geneFilename, path[0]);
    strcat(geneFilename, "/");
    strcat(geneFilename, project[0]);
    strcat(geneFilename, underscore);
    strcat(geneFilename, geneDesig);
    strcat(geneFilename, underscore);
    strcat(geneFilename, prefixChar);
    strcat(geneFilename, geneSuffix);
    //printf("%s\n", geneFilename);
    FILE* infile = fopen(geneFilename, "rb");
    assert(infile);
    Org* genome = makeOrg();
    genome->geneCounts[0] = nGenes[0];
    genome->genes[0] = malloc(sizeof(Gene*) * genome->geneCounts[0]);
    Org* mutant = makeOrg();
    mutant->geneCounts[0] = nGenes[0];
    mutant->genes[0] = malloc(sizeof(Gene*) * mutant->geneCounts[0]);
    for(int i = 0; i < nGenes[0]; i++)
    {
        fseek(infile, (long)(offset * geno[i]), SEEK_SET);
        genome->genes[0][i] = readGene(infile); 
    } 
    //double** tmp_meanP = malloc(sizeof(double*) * nbmutants[0]); //tmp_sdP[trait][mutant] that's how it works
    //double** tmp_sdP = malloc(sizeof(double*) * num_T);
    int position = 0;

    for(int i = 0; i < num_T; i++) trait_optima[i] = environment[0];
    
    for (int nbmutations_index=0;nbmutations_index<nbStepsModel[0];nbmutations_index++)
    {
        for(int mut_count = 0; mut_count < nbmutants[nbmutations_index]; mut_count++)
        {
            /*if ((mut_count+1) % (nbmutants[mut_count]/2) == 0)
            {

                printf("NbMutationsStep = %i/%i  |  mutant = %i/%i\n",
                    nbmutations_index+1,
                    nbStepsModel[0],
                    mut_count+1,
                    nbmutants[nbmutations_index]);    
                //printf("modulo -> %i\n",(mut_count+1) % (nbmutants[mut_count]/2));
            }*/
            
            double* ps = malloc(sizeof(double) * reps[nbmutations_index]);
            double* ws = malloc(sizeof(double) * reps[nbmutations_index]);

            //printf("nbmutations_index = %i  | nbmutants[nbmutations_index] = %i  |  mut_count = %i  |  position %i\n",nbmutations_index,nbmutants[nbmutations_index],mut_count,position);
            meanW[position] = 0; // mean fitness
            meanP[position] = 0; // mean pheno
            // get Mutant
            if (trace[0] & !nbmutations_index) 
            {
                nbNewMutationsNeeded = nbmutations[nbmutations_index] - nbmutations[nbmutations_index-1];
            } else
            {
                nbNewMutationsNeeded = nbmutations[nbmutations_index];    
            }
            
            if (nbNewMutationsNeeded == 0) 
            {
                mutant = cloneOrg(genome);
            } else 
            {
                int changes = -1;
                int tries = 0;
                //int randPloidy = gsl_rng_uniform_int(rg, ploidy); 
                while(changes != nbNewMutationsNeeded)
                {
                    mutant = cloneOrg(genome);

                    GeneToMutate = rand() % nGenes[0];
                    //printf("%i < %i\n",GeneToMutate,nGenes[0]);
                    //printf("genome->haplo = %i  | mutant->haplo = %i\n",genome->haplotype,mutant->haplotype);
                    changes = ReDevelop_substitutions(mutant->genes[0][GeneToMutate], 0, 0, rg, NULL, NULL, nbmutations[nbmutations_index]);

                    if(changes != nbNewMutationsNeeded) deleteOrg(mutant);

                    if (tries > 100)
                    {
                      printf("tries : %d --> Expected nb mutations = %i  |  Total Mut rate = %f  | got %i mutations at last try\n",tries,nbmutations[nbmutations_index],mu_ts+mu_trait+mu_dbm+mu_cs+mu_small_K+mu_big_K,changes);  
                    } 
                    if (tries > 1000000)
                    {
                        printf("Tried to mutate 1000000 without success. Abort!\n");abort();
                    }
                    tries++;
                }
            }

            // loop over reps
            for(int rep = 0; rep < reps[nbmutations_index]; rep++)
            {
                if(rep % 500 == 0) R_CheckUserInterrupt();
                develop(mutant, 0, NULL, -1);
                ws[rep] = mutant->fitness;
                meanW[position] += mutant->fitness;
                ps[rep] = mutant->traits[0];
                meanP[position] += mutant->traits[0];
            }

            meanW[position]  /= reps[nbmutations_index];
            meanP[position]  /= reps[nbmutations_index];
            medianW[position] = median(ws,reps[nbmutations_index]);
            medianP[position] = median(ps,reps[nbmutations_index]);
            for (int rep = 0; rep < reps[nbmutations_index]; rep++)
            {
                sdW[position] += pow(ws[rep] - meanW[position], 2);
                sdP[position] += pow(ps[rep] - meanP[position], 2);
            }

            sdW[position] = pow(sdW[position]/reps[nbmutations_index],0.5);
            sdP[position] = pow(sdP[position]/reps[nbmutations_index],0.5);    

            free(ps);
            free(ws);
            position ++;
        }
        if (trace[0]) genome = cloneOrg(mutant); // Last mutant becomes the new genome and further mutates.
    }
    free(geneFilename);
    //free(genome->genes[0]);
    //free(trait_optima); // that's what cleanup() does
    deleteOrg(mutant);
    deleteOrg(genome);
    cleanup();
    fclose(infile);
}

void DR(char** path, char** project, int* intConstants, double* doubleConstants, int* replicate, int* geno, int* nGenes, int* reps, double* meanW, double* sdW, double* meanP, double* sdP)
{
    //printf("replicate: %i", replicate[0]);
    ploidy = intConstants[0];
    max_b = intConstants[1];
    k_weak = intConstants[2];
    num_T = intConstants[3];
    chromosomes = intConstants[4];
    dualPhenoGenes = intConstants[5];
    sex = intConstants[6];
    EnvSignal = intConstants[7];
    
    max_time = doubleConstants[0];
    basal = doubleConstants[1];
    transcription = doubleConstants[2];
    translation = doubleConstants[3];
    mRNA_decay = doubleConstants[4];
    protein_decay = doubleConstants[5];
    eps = doubleConstants[6];
    omega_squared = doubleConstants[7];
    Env = doubleConstants[8];
    Env_env = doubleConstants[9];
    Env_perf = doubleConstants[10];

    setup();
    int offset = sizeof(int) * (8 + max_b + num_T) + sizeof(double) * (4 + max_b + num_T);
    //printf("%d\n", offset);
    char* geneFilename = calloc(500, sizeof(char));
    char geneDesig[] = "bin";
    char geneSuffix[] = ".bin";
    char underscore[] = "_";
    char slash[] = "/";
    char prefixChar [4];
    prefixChar[0] = (char)(replicate[0] / 100 + 48);
    prefixChar[1] = (char)((replicate[0] /10) % 10 + 48);
    prefixChar[2] = (char)(replicate[0] % 10 + 48);
    prefixChar[3] = '\0';
    strcat(geneFilename, path[0]);
    strcat(geneFilename, "/");
    strcat(geneFilename, project[0]);
    strcat(geneFilename, underscore);
    strcat(geneFilename, geneDesig);
    strcat(geneFilename, underscore);
    strcat(geneFilename, prefixChar);
    strcat(geneFilename, geneSuffix);
    //printf("%s\n", geneFilename);
    FILE* infile = fopen(geneFilename, "rb");
    assert(infile);
    Org* genome = makeOrg();
    genome->geneCounts[0] = nGenes[0];
    genome->genes[0] = malloc(sizeof(Gene*) * genome->geneCounts[0]);
    for(int i = 0; i < nGenes[0]; i++)
    {
        fseek(infile, (long)(offset * geno[i]), SEEK_SET);
        genome->genes[0][i] = readGene(infile);
    }
    double* ws = malloc(sizeof(double) * reps[0]);
    double* ps = malloc(sizeof(double*) * reps[0]);
    meanW[0] = 0;
    sdW[0] = 0;
    meanP[0] = 0;
    sdP[0] = 0;
    for(int i = 0; i < num_T; i++) 
    {
        trait_optima[i] = Env;
        trait_optima_env[i] = Env_env;
        trait_optima_perf[i] = Env_perf;
    }
    for(int i = 0; i < reps[0]; i++)
    {
        if(i % 100 == 0) R_CheckUserInterrupt();
        develop(genome, 0, NULL, -1);
        
        meanW[0] += genome->fitness;
        ws[i] = genome->fitness;
        
        meanP[0] += genome->traits[0];
        ps[i] = genome->traits[0]; 
    }

    meanW[0] /= reps[0];
    meanP[0] /= reps[0];
    for (int i = 0; i < reps[0]; i++) 
    {
        sdW[0] += pow(ws[i] - meanW[0], 2);   
        sdP[0] += pow(ps[i] - meanP[0], 2);   
    }
    sdW[0] = sqrt(sdW[0]/reps[0]);
    sdP[0] = sqrt(sdP[0]/reps[0]);

    free(ps);
    free(ws);
    //for(int i = 0; i < reps[0]; i++) free(fitness[i]);
    //free(trait_optima);
    free(geneFilename);
    //free(genome->genes[0]);
    cleanup();
    deleteOrg(genome);
    fclose(infile);
}



void MeasureFit (char** path, char** project, int* intConstants, double* doubleConstants, int* replicate, int* geno, int* nGenes, int* reps, double* ArMeanFitness, double* sdFitness, double* GeoMeanFitness, double* GeoUp95CIFitness, double* GeoLow95CIFitness, double* environments, int* nbenvironments, double* meanP, double* sdP)
{
     //Mutational constants
    mu_nuc = 0.001;
    n_cis = 12;
    L_cis = 500;
    mu_pheno = 100;
    sigma_ts = 0.5;
    sigma_trait = 0.5;
    sigma_cs = 0.5;
    pos_cs = 0.9;
    initial_mean_cs = 0.5;
    k_mut = 4;
    mu_ts = mu_nuc * mu_pheno;
    mu_trait = mu_nuc * mu_pheno;
    mu_dbm = mu_nuc * mu_pheno;
    mu_cs = mu_nuc * 50;
    mu_small_K = mu_nuc * n_cis;
    mu_big_K = mu_nuc * k_mut * L_cis * factorial(n_cis) * pow(0.75, k_mut) * pow(0.25, n_cis-k_mut) / factorial(k_mut) / factorial(n_cis - k_mut);

    ploidy = intConstants[0];
    max_b = intConstants[1];
    k_weak = intConstants[2];
    num_T = intConstants[3];
    chromosomes = intConstants[4];
    dualPhenoGenes = intConstants[5];
    sex = intConstants[6];
    EnvSignal = intConstants[7];
    
    max_time = doubleConstants[0];
    basal = doubleConstants[1];
    transcription = doubleConstants[2];
    translation = doubleConstants[3];
    mRNA_decay = doubleConstants[4];
    protein_decay = doubleConstants[5];
    eps = doubleConstants[6];
    omega_squared = doubleConstants[7];
    //environment = doubleConstants[8];
    // printf("line 247: %f\n",environment);

    
    setup();
    
    int offset = sizeof(int) * (8 + max_b + num_T) + sizeof(double) * (4 + max_b + num_T);
    //printf("%d\n", offset);
    char* geneFilename = calloc(500, sizeof(char));
    char geneDesig[] = "bin";
    char geneSuffix[] = ".bin";
    char underscore[] = "_";
    char slash[] = "/";
    char prefixChar [4];
    prefixChar[0] = (char)(replicate[0] / 100 + 48);
    prefixChar[1] = (char)((replicate[0] /10) % 10 + 48);
    prefixChar[2] = (char)(replicate[0] % 10 + 48);
    prefixChar[3] = '\0';
    strcat(geneFilename, path[0]);
    strcat(geneFilename, "/");
    strcat(geneFilename, project[0]);
    strcat(geneFilename, underscore);
    strcat(geneFilename, geneDesig);
    strcat(geneFilename, underscore);
    strcat(geneFilename, prefixChar);
    strcat(geneFilename, geneSuffix);
    //printf("%s\n", geneFilename);
    FILE* infile = fopen(geneFilename, "rb");
    assert(infile);
    Org* genome = makeOrg();
    genome->geneCounts[0] = nGenes[0];
    genome->genes[0] = malloc(sizeof(Gene*) * genome->geneCounts[0]);
    Org* mutant = makeOrg();
    mutant->geneCounts[0] = nGenes[0];
    mutant->genes[0] = malloc(sizeof(Gene*) * mutant->geneCounts[0]);
    for(int i = 0; i < nGenes[0]; i++)
    {
        fseek(infile, (long)(offset * geno[i]), SEEK_SET);
        genome->genes[0][i] = readGene(infile);
    }

    // double* meanFitness = malloc(sizeof(double)*nbmutants[0]); Memory already allocated from R
    // double* sdFitness = malloc(sizeof(double)*nbmutants[0]); Memory already allocated from R
    int totalNbReps = nbenvironments[0]*reps[0];
    
    //printf("mut_count: %i\n", mut_count);
    // Set some variables to zero
    ArMeanFitness[0] = 0; // Arithmetic mean
    GeoMeanFitness[0] = 1; // Geometric mean
    
    for(int i = 0; i < num_T; i++)
    {
        meanP[i] = 0;
        sdP[i] = 0;
        trait_optima[i] = environment;
    }
    double** ps = malloc(sizeof(double*) * num_T);
    for(int i = 0; i < num_T; i++) ps[i] = malloc(sizeof(double) * reps[0]);
    
    double* fitness = malloc(sizeof(double) * totalNbReps);
    // for Geometric Mean use only
    double* LogFit = malloc(sizeof(double) * totalNbReps);
    meanLogFit = 0;
    
    // loop over reps
    int ind_count = 0;
    for(int i = 0; i < reps[0]; i++)
    {
        if(i % 50 == 0) R_CheckUserInterrupt();
        for (int i_env = 0; i_env < nbenvironments[0] ; i_env++)
        {
            // Set environment
            environment = environments[i_env];
            for (int i = 0; i < num_T; i++) trait_optima[i] = environment;
            develop(genome, 0, NULL, -1);
             // Arithmetic mean
            ArMeanFitness[0] += genome->fitness;
            fitness[ind_count] = genome->fitness;
             // Geometric mean
            LogFit[ind_count] = log(genome->fitness);
            meanLogFit += LogFit[ind_count];
            for(int j = 0; j < num_T; j++)
            {
                meanP[j] += genome->traits[j]; // mean (see division afterward) trait value for each trait
                ps[j][i] = genome->traits[j]; // trait value for each trait j at each replicate i
            }
            ind_count++;
        }
    }
    // Arithmetic mean
    sdFitness[0] = 0;
    ArMeanFitness[0] /= totalNbReps;
    for (int i = 0; i < totalNbReps ; i++) sdFitness[0] += pow(fitness[i] - ArMeanFitness[0], 2);
    sdFitness[0] = sqrt(sdFitness[0]/totalNbReps);
    

    // Geometric mean
    meanLogFit /= totalNbReps;
    GeoMeanFitness[0] = exp(meanLogFit);
    
    varLogFit = 0;
    for (int i = 0; i < totalNbReps; i++) varLogFit += pow(LogFit[i] - meanLogFit, 2);
    varLogFit /= totalNbReps;
    seLogFit = sqrt(varLogFit / totalNbReps);
    GeoUp95CIFitness[0] = exp(meanLogFit + seLogFit);
    GeoLow95CIFitness[0] = exp(meanLogFit - seLogFit);

    // Phenotype
    for(int i = 0; i < num_T; i++) meanP[i] /= totalNbReps;    
    for(int i = 0; i < num_T; i++) // i is trait
    {
        for(int j = 0; j < totalNbReps; j++) // j is replicate
        {
            sdP[i] += pow(meanP[i] - ps[i][j], 2); // add all square differences between the mean trait and the trait of each replicate. It is therefore the variance times the number of replicates
        }
    }
    for(int i = 0; i < num_T; i++)
    {
        sdP[i] /= totalNbReps; // sdP is the average variance of each replicate
        sdP[i] = sqrt(sdP[i]); // sdP is the average standard deviation of each replicate.
        free(ps[i]);
    }
    free(ps);

    free(fitness);
    free(LogFit);

    free(geneFilename);
    //free(genome->genes[0]);
    //free(trait_optima); // that's what cleanup() does
    deleteOrg(genome);
    cleanup();
    fclose(infile);
}

