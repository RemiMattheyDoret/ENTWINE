/*
 Declares global constants and variables.
 Contains main(), which:
 - Sets up RNG
 - Handles command-line input
 - Creates files
 - Manages the major conditional (restarting=={0,1,2})
 - Finds an initial genotype/restarts from previous population
 - Runs evolution
 - Cleans up memory
*/


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_ieee_utils.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "constants.h"
#include "G.h"
#include "develop.h"
#include "utilities.h"

//Global variables for RNG
const gsl_rng_type* rT;
gsl_rng* rg;

//Global variables used to id mutations, duplications, and novel haplotypes
int next;
int nextLocus;
int nextHaploid;

//Declarations for global constants
int expressInKiloGenerations;
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
int signals;
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
double low_environment;
double high_environment;
int prefix;
double K_vals [4];
int dualPhenoGenes;
int treatment;
int EnvSignal;
int EnvHeterogeneity; // 0 = no heterogeneity, 1 = spatial heterogeneity, 2 = temporal heterogeneity
double pleiotropy_proba;

int main(int argc, char * argv[])
{
    //Variables set by parameters but without the need for global scope
    int N;
    int GENS;
    char* project = malloc(sizeof(char) * 100);
    char* treatment = malloc(sizeof(char) * 100);
    char* home = malloc(sizeof(char) * 1000);
    
    //Restarting == 0; run new populations from scratch
    //Restarting == 1; restart populations from a checkpoint, appending to files
    int restarting;
    int whichCheckpoint;
    
    // Structure for dictating how command-line arguments are received.
    struct option longopts[] = {
        { "N",                      required_argument, NULL, 1},
        {"GENS",                    required_argument, NULL, 2},
        {"num_T",                   required_argument, NULL, 3},
        {"num_G",                   required_argument, NULL, 4},
        {"num_P",                   required_argument, NULL, 5},
        {"num_SE",                  required_argument, NULL, 6},
        {"chromosomes",             required_argument, NULL, 7},
        {"ploidy",                  required_argument, NULL, 8},
        {"sex",                     required_argument, NULL, 9},
        {"rec_rate",                required_argument, NULL, 10},
        {"max_b",                   required_argument, NULL, 11},
        {"EnvSignal",               required_argument, NULL, 12},

        {"k_mut",                   required_argument, NULL, 14},
        {"k_weak",                  required_argument, NULL, 15},
        {"mu_nuc",                  required_argument, NULL, 16},
        {"L_cis",                   required_argument, NULL, 17},
        {"n_cis",                   required_argument, NULL, 18},
        {"omega_squared",           required_argument, NULL, 19},
        {"low_environment",         required_argument, NULL, 20},
        {"high_environment",        required_argument, NULL, 21},
        {"mu_del",                  required_argument, NULL, 22},
        {"mu_dup",                  required_argument, NULL, 23},
        {"trans_mutational_target", required_argument, NULL, 24},
        {"trait_mutational_target", required_argument, NULL, 25},
        {"dbm_mutational_target",   required_argument, NULL, 26},
        {"cis_mutational_target",   required_argument, NULL, 27},
        {"sig_mutational_target",   required_argument, NULL, 28},
        {"sigma_ts",                required_argument, NULL, 29},
        {"sigma_trait",             required_argument, NULL, 30},
        {"sigma_cs",                required_argument, NULL, 31},
        {"sigma_signalK",           required_argument, NULL, 32},
        {"pos_cs",                  required_argument, NULL, 33},
        {"initial_mean_cs",         required_argument, NULL, 34},
        {"initial_mean_ts",         required_argument, NULL, 35},
        {"initial_sigma_trait",     required_argument, NULL, 36},
        {"initial_mean_signalK",    required_argument, NULL, 37},
        {"max_time",                required_argument, NULL, 38},
        {"basal",                   required_argument, NULL, 39},
        {"transcription",           required_argument, NULL, 40},
        {"translation",             required_argument, NULL, 41},
        {"mRNA_decay",              required_argument, NULL, 42},
        {"protein_decay",           required_argument, NULL, 43},
        {"dual_pheno_genes",        required_argument, NULL, 44},
        {"eps",                     required_argument, NULL, 45},
        {"project",                 required_argument, NULL, 46},
        {"prefix",                  required_argument, NULL, 47},
        {"restarting",              required_argument, NULL, 48},
        {"whichCheckpoint",         optional_argument, NULL, 49},
        {"treatment",               optional_argument, NULL, 50},
        {"num_ST",                  required_argument, NULL, 51},
        {"num_SP",                  required_argument, NULL, 52},
        {"pleiotropy_proba",        required_argument, NULL, 53},
        {"home",                    required_argument, NULL, 54},
        {"EnvHeterogeneity",        required_argument, NULL, 55},
        {"seed",                    required_argument, NULL, 56},
        {"expressInKiloGenerations",required_argument, NULL, 57},
        { 0, 0, 0, 0 }
    };
    
    //Process command-line arguments
    int c;
    while( (c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1)
    {
        switch(c)
        {
            case 1:
                N = atoi(optarg);
                break;
             case 2:
                GENS = atoi(optarg);
                break;
            case 3:
                num_T = atoi(optarg);
                break;
            case 4:
                num_G = atoi(optarg);
                break;
            case 5:
                num_P = atoi(optarg);
                break;
            case 6:
                num_SE = atoi(optarg);
                break;
            case 7:
                chromosomes = atoi(optarg);
                break;
            case 8:
                ploidy = atoi(optarg);
                break;
            case 9:
                sex = atoi(optarg);
                break;
            case 10:
                rec_rate = atof(optarg);
                break;
            case 11:
                max_b = atoi(optarg);
                break;
            case 12:
                EnvSignal = atoi(optarg);
                break;
            case 14:
                k_mut = atoi(optarg);
                break;
            case 15:
                k_weak = atoi(optarg);
                break;
            case 16:
                mu_nuc = atof(optarg);
                break;
            case 17:
                L_cis = atoi(optarg);
                break;
            case 18:
                n_cis = atoi(optarg);
                break;
            case 19:
                omega_squared = atof(optarg);
                break;
            case 20:
                low_environment = atof(optarg);
                break;
            case 21:
                high_environment = atof(optarg);
                break;
            case 22:
                mu_del = atof(optarg);
                break;
            case 23:
                mu_dup = atof(optarg);
                break;
            case 24:
                trans_mutational_target = atof(optarg);
                break;
            case 25:
                trait_mutational_target = atof(optarg);
                break;
            case 26:
                dbm_mutational_target = atof(optarg);
                break;
            case 27:
                cis_mutational_target = atof(optarg);
                break;
            case 28:
                sig_mutational_target = atof(optarg);
                break;
            case 29:
                sigma_ts = atof(optarg);
                break;
            case 30:
                sigma_trait = atof(optarg);
                break;
            case 31:
                sigma_cs = atof(optarg);
                break;
            case 32:
                sigma_signalK = atof(optarg);
                break;
            case 33:
                pos_cs = atof(optarg);
                break;
            case 34:
                initial_mean_cs = atof(optarg);
                break;
            case 35:
                initial_mean_ts = atof(optarg);
                break;
            case 36:
                initial_sigma_trait = atof(optarg);
                break;
            case 37:
                initial_mean_signalK = atof(optarg);
                break;
            case 38:
                max_time = atof(optarg);
                break;
            case 39:
                basal = atof(optarg);
                break;
            case 40:
                transcription = atof(optarg);
                break;
            case 41:
                translation = atof(optarg);
                break;
            case 42:
                mRNA_decay = atof(optarg);
                break;
            case 43:
                protein_decay = atof(optarg);
                break;
            case 44:
                dualPhenoGenes = atoi(optarg);
                break;
            case 45:
                eps = atof(optarg);
                break;
            case 46:
                strcpy(project, optarg);
                break;
            case 47:
                prefix = atoi(optarg);
                break;
            case 48:
                restarting = atoi(optarg);
                break;
            case 49:
                whichCheckpoint = atoi(optarg);
                break;
            case 50:
                treatment = optarg;
                break;      
            case 51:
                num_ST = atoi(optarg);
                break;
            case 52:
                num_SP = atoi(optarg);
                break;
            case 53:
                pleiotropy_proba = atoi(optarg);
                break;
            case 54:
                strcpy(home, optarg);
                break;
            case 55:
                EnvHeterogeneity = atoi(optarg);
                break;
            case 56:
                sValue = atoi(optarg);
            case 57:
                expressInKiloGenerations = atoi(optarg) != 0;
        }
    }
    
    setup();
    
    //Set up output files
    //char home[] = "/Users/remi/Documents/1.Plasticity/WestGrid/Simulations/SimulationOutputs/";
    //char home[] = "/home/draghi/data/";
    char* filePath = calloc(500, sizeof(char));
    strcat(filePath, home);
    strcat(filePath, project);
    char* systemCall = calloc(500, sizeof(char));
    if(restarting == 0)
    {
        strcat(systemCall, "mkdir ");
        strcat(systemCall, filePath);
        system(systemCall);
    }
    strcat(filePath, "/");
	char* fitFilename = calloc(500, sizeof(char));
    char* alleleFilename = calloc(500, sizeof(char));
    char* geneFilename = calloc(500, sizeof(char));
    char* binFilename = calloc(500, sizeof(char));
    char* popFilename = calloc(500, sizeof(char));
    char* configFilename = calloc(500, sizeof(char));
    char* haploFilename = calloc(500, sizeof(char));
    char* envFilename = calloc(500, sizeof(char));
	char prefixChar [4];
	char textSuffix[] = ".txt";
    char binSuffix[] = ".bin";
    char underscore[] = "_";
    char fitDesig[] = "fitness";
    char geneDesig[] = "gene";
    char alleleDesig[] = "allele";
    char binDesig[] = "bin";
    char popDesig[] = "pop";
    char binPopDesig[] = "binPop";
    char mtSeedDesig[] = "mt_seed";
    char configDesig[] = "config";
    char haploDesig[] = "haplo";
    char envDesig[] = "env";
	prefixChar[0] = (char)(prefix / 100 + 48);
	prefixChar[1] = (char)((prefix /10) % 10 + 48);
	prefixChar[2] = (char)(prefix % 10 + 48);
	prefixChar[3] = '\0';
    
    char* name;
    char* desig;
    char* suffix;
    for(int i = 0; i < 8; i++)
    {
        if(i == 0)
        {
            name = fitFilename;
            desig =fitDesig;
            suffix = textSuffix;
        }
        if(i == 1)
        {
            name = alleleFilename;
            desig =alleleDesig;
            suffix = textSuffix;
        }
        if(i == 2)
        {
            name = geneFilename;
            desig =geneDesig;
            suffix = textSuffix;
        }
        if(i == 3)
        {
            name = binFilename;
            desig =binDesig;
            suffix = binSuffix;
        }
        if(i == 4)
        {
            name = popFilename;
            desig = popDesig;
            suffix = textSuffix;
        }
        if(i == 5)
        {
            name = configFilename;
            desig = configDesig;
            suffix = textSuffix;
        }
        if(i == 6)
        {
            name = haploFilename;
            desig = haploDesig;
            suffix = textSuffix;
        }
        if(i == 7)
        {
            name = envFilename;
            desig = envDesig;
            suffix = textSuffix;
        }
        strcat(name, filePath);
        strcat(name, project);
        strcat(name, underscore);
        strcat(name, desig);
        strcat(name, underscore);
        strcat(name, prefixChar);
        strcat(name, suffix);
    }
    
    FILE* fitFile;
    FILE* alleleFile;
    FILE* geneFile;
    FILE* binFile;
    FILE* popFile;
    FILE* configFile;
    FILE* haploFile;
    FILE* envFile;

    fitFile = fopen(fitFilename, "a");
    alleleFile = fopen(alleleFilename, "a");
    geneFile = fopen(geneFilename, "a");
    binFile = fopen(binFilename, "ab");
    popFile = fopen(popFilename, "a");
    haploFile = fopen(haploFilename, "a");
    envFile = fopen(envFilename, "a");
    printf("fitFilename = %s\n", fitFilename);
    assert(fitFile);
    assert(alleleFile);
    assert(geneFile);
    assert(binFile);
    assert(popFile);
    assert(haploFile);
    assert(envFile);
    
    //Variables for population
    Org** pop = malloc(sizeof(Org*) * N);
    Org** newPop = malloc(sizeof(Org*) * N);
    Org* ancestor;
    Org** temp;
    
    //Variables for selection
    double* intervals = malloc(sizeof(double) * N);
    double tFit, mFit;
    int picked;
    
    //Initialize these variables for record-keeping (though they automatically initilize to zero anyway).
    next = 0;
    nextLocus = 0;
    nextHaploid = 0;

    int t = 0;
    
    //Print all constants to file
    configFile = fopen(configFilename, "w");
    fprintf(configFile, "%d\n", sValue);
    for(int i = 1; i < argc; i++)
    {
        fprintf(configFile, "%s\n", argv[i]);
    }
    fclose(configFile);
    // Print output file headers
    fprintf(alleleFile, "time id count\n");
    fprintf(fitFile, "prefix rep time fitness\n");
    fprintf(geneFile, "id locus parentLocus parentAllele time type map dbm trans whichSignal signalK signalN ");
    for(int i = 0; i < num_T; i++) fprintf(geneFile, "target%d effect%d ", i+1, i+1);
    for(int i = 0; i < max_b; i++) fprintf(geneFile, "c%d ", i);
    for(int i = 0; i < max_b; i++) fprintf(geneFile, "k%d ", i);
    fprintf(geneFile, "\n");
    fprintf(popFile, "time individual haplotype_one\n");
    fprintf(haploFile, "id time parentID haplotype\n");
    //double tmpT1, tmpT2;
    double tmpFit;
    int counter = 0;
    
    for(int i = 0; i < num_T; i++)
    {
        trait_optima[i] = low_environment;
    }
    
    //Generate candidate starting genotypes
    printf("%d\n", sValue);
    do
    {
        ancestor = randomOrg(rg);
        tmpFit = 0;
        for(int i = 0; i < 1000; i++)
        {
            develop(ancestor, 0, 0, -1);
            tmpFit += ancestor->fitness;
        }            
        tmpFit /= 1000;
        printf("# %d %f\n", counter, tmpFit);
        fflush(stdout);
        counter++;
        if(tmpFit < 0.15) deleteOrg(ancestor);

    }while(tmpFit < 0.15);
    
    //Record ancestor and set counters
    for(int i = 0; i < ancestor->geneCounts[0]; i++)
    {
        for(int j = 0; j < ploidy; j++)
        {
            ancestor->genes[j][i]->number = next;
            ancestor->genes[j][i]->locus = nextLocus;
        }
        printGene(ancestor->genes[0][i], geneFile, binFile, 0, -1, -1);
        next++;
        nextLocus++;
    }
    ancestor->haplotype = nextHaploid;
    nextHaploid++;
    printHaplotype(ancestor, haploFile, 0, -1);
    for(int i = 0; i < N; i++)
    {
        pop[i] = cloneOrg(ancestor);
    }
    
    if (EnvHeterogeneity == 0) // no environmental heterogeneity
    {
        for(int tr = 0; tr < num_T; tr++)
        {
           trait_optima[tr] = low_environment;
        }
    }
    

    for(int i = 0; i < N; i++)
    {
        fprintf(popFile, "%d %d ", t, i);
        for(int j = 0; j < ploidy; j++)
        {
            for(int k = 0; k < pop[i]->geneCounts[j]; k++)
            {
                fprintf(popFile, "%d,", pop[i]->genes[j][k]->number);
            }
            fprintf(popFile, " ");
        }
        fprintf(popFile, "\n");
    }
    fflush(popFile);
    char* binPopFilename = calloc(500, sizeof(char));

    char popNum[10];
    if (expressInKiloGenerations)
    {
        int k = t / 1000;
        popNum[0] = (char)((k) / 100 + 48);
        popNum[1] = (char)(((k) /10) % 10 + 48);
        popNum[2] = (char)((k) % 10 + 48);
        popNum[3] = 'k';
        popNum[4] = '\0';
    } else
    {
        sprintf(popNum, "%d", t);
    }

    strcat(binPopFilename, filePath);
    strcat(binPopFilename, project);
    strcat(binPopFilename, underscore);
    strcat(binPopFilename, binPopDesig);
    strcat(binPopFilename, underscore);
    strcat(binPopFilename, prefixChar);
    if(restarting == 2)
    {
        strcat(binPopFilename, underscore);
        strcat(binPopFilename, treatment);
    }
    strcat(binPopFilename, underscore);
    strcat(binPopFilename, popNum);
    strcat(binPopFilename, binSuffix);
    FILE* binPopFile;
    binPopFile = fopen(binPopFilename, "wb");
    for(int i = 0; i < N; i++) writeOrgToFile(pop[i], binPopFile);
    fclose(binPopFile);
    free(binPopFilename);
    
    char* mtSeedFilename = calloc(500, sizeof(char));
    strcat(mtSeedFilename, filePath);
    strcat(mtSeedFilename, project);
    strcat(mtSeedFilename, underscore);
    strcat(mtSeedFilename, mtSeedDesig);
    strcat(mtSeedFilename, underscore);
    strcat(mtSeedFilename, prefixChar);
    strcat(mtSeedFilename, underscore);
    strcat(mtSeedFilename, popNum);
    strcat(mtSeedFilename, binSuffix);
    FILE* mtSeedFile;
    mtSeedFile = fopen(mtSeedFilename, "wb");
    gsl_rng_fwrite(mtSeedFile, rg);
    fclose(mtSeedFile);
    free(mtSeedFilename);


    // Clean up
    for(int i = 0; i < N; i++)
    {
        deleteOrg(pop[i]);
    }
    if(restarting == 0) deleteOrg(ancestor);
    
    cleanup();
    free(pop);
    free(newPop);
    free(intervals);
    free(project);
    fclose(fitFile);
    fclose(alleleFile);
    fclose(geneFile);
    fclose(binFile);
    fclose(popFile);
    fclose(haploFile);
    fclose(envFile);
    free(fitFilename);
    free(haploFilename);
    free(alleleFilename);
    free(geneFilename);
    free(envFilename);
    free(filePath);
    free(systemCall);
    free(treatment);
    return 0;
}

