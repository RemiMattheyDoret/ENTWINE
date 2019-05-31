/*
    Contains a number of functions for managing structs to contain gene and organism data, allocating/deallocating memory,
    calculating mutations, writing genotypes to files, and producing offspring.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "G.h"
#include "constants.h"
#include "utilities.h"
#include "develop.h"

void setup()
{
    if (sValue == -1)
    {
        FILE* m_random;
        m_random = fopen("/dev/random", "rb");
        fgets((char*)(&sValue), 4, m_random); 
        fclose(m_random);
    }	
    
    //Generate mersenne twister random number generator and seed
	rT = gsl_rng_mt19937;
	rg = gsl_rng_alloc(rT);
	gsl_rng_set(rg, sValue);
    
    //Calculate runtime constants
    mu_ts = mu_nuc * trans_mutational_target;
    mu_trait = mu_nuc * trait_mutational_target;
    mu_dbm = mu_nuc * dbm_mutational_target;
    mu_cs = mu_nuc * cis_mutational_target;
    mu_small_K = mu_nuc * n_cis;
    mu_big_K = mu_nuc * k_mut * L_cis * factorial(n_cis) * pow(0.75, k_mut) * pow(0.25, n_cis-k_mut) / factorial(k_mut) / factorial(n_cis - k_mut);
    mu_signalK = mu_nuc * sig_mutational_target;
    mu_signalN = mu_nuc * sig_mutational_target;
    //mu_whichSignal = mu_nuc * sig_mutational_target;
    mu_whichSignal = 0;
    
    //Hardcoded constants
    K_vals[0] = 5.00000;
    K_vals[1] = 36.94528;
    K_vals[2] = 272.99075;
    K_vals[3] = 2017.14397;
    
    basic_signal_rate = 80;
    trait_optima = malloc(sizeof(double) * num_T);
    trait_optima_env = malloc(sizeof(double) * num_T);
    trait_optima_perf = malloc(sizeof(double) * num_T);    

/*signals = 0;
    
    if(EnvSignal == 1)
    {
        signals = 1;
    }
    if(traitFeedback == 1) signals++;
    */
}

void cleanup()
{
    free(trait_optima);
}

// Allocates memory for a gene
Gene* makeGene(int type)
{
    Gene* ret = malloc(sizeof(Gene));
    assert(ret);
    ret->type = type;
    ret->Ks = malloc(sizeof(int) * max_b);
    ret->cs = malloc(sizeof(double) * max_b);
    assert(ret->Ks);
    assert(ret->cs);
    ret->number = -1;
    ret->locus = -2;
    ret->targets = malloc(sizeof(int) * num_T);
    ret->traitEffects = malloc(sizeof(double) * num_T);
    assert(ret->targets);
    assert(ret->traitEffects);

    return ret;
}

//Frees memory for a gene
void deleteGene(Gene* g)
{
    free(g->Ks);
    free(g->cs);
    free(g->targets);
    free(g->traitEffects);
    free(g);
}

//Allocates memory for an organism
Org* makeOrg()
{
    Org* ret = malloc(sizeof(Org));
    assert(ret);
    ret->geneCounts = malloc(sizeof(int) * ploidy);
    ret->genes = malloc(sizeof(Gene**) * ploidy);
    ret->traits = malloc(sizeof(double) * num_T);
    assert(ret->geneCounts);
    assert(ret->genes);
    assert(ret->traits);
    return ret;
}

//Makes a diploid homozygote from the haplotype indicated by side
Org* makeHomozygote(Org* src, int side)
{
    Org* dest = makeOrg();
    dest->genes[0] = malloc(sizeof(Gene*) * src->geneCounts[side]);
    dest->genes[1] = malloc(sizeof(Gene*) * src->geneCounts[side]);
    assert(dest->genes[0]);
    assert(dest->genes[1]);
    for(int i = 0; i < src->geneCounts[side]; i++)
    {
        dest->genes[0][i] = cloneGene(src->genes[side][i]);
        dest->genes[1][i] = cloneGene(src->genes[side][i]);
    }
    dest->geneCounts[0] = src->geneCounts[side];
    dest->geneCounts[1] = src->geneCounts[side];
    return dest;
}

//Free memory allocated by an organism
void deleteOrg(Org* o)
{
    for(int i = 0; i < ploidy; i++)
    {
        for(int j = 0; j < o->geneCounts[i]; j++)
        {
            deleteGene(o->genes[i][j]);
        }
        free(o->genes[i]);
    }
    free(o->genes);
    free(o->geneCounts);
    free(o->traits);
    free(o);
}

//Copy a gene without mutation
Gene* cloneGene(Gene* g)
{
    Gene* ret = makeGene(g->type);
    assert(ret);
    ret->map = g->map;
    for(int i = 0; i < max_b; i++)
    {
        ret->Ks[i] = g->Ks[i];
        ret->cs[i] = g->cs[i];
    }
    ret->dbm = g->dbm;
    ret->transEffect = g->transEffect;
    ret->whichSignal = g->whichSignal;
    ret->signalK = g->signalK;
    ret->signalN = g->signalN;
    for(int i = 0; i < num_T; i++)
    {
        ret->targets[i] = g->targets[i];
        ret->traitEffects[i] = g->traitEffects[i];
    }
    ret->number = g->number;
    ret->locus = g->locus;
    return ret;
}

// Create an organism with num_G regualtory genes, num_P phenotype genes, num_SE Environment sensor genes, num_ST Trait sensor genes and num_SP Perfomence sensor genes
// Also contains code for random formation of genes
Org* randomOrg(gsl_rng* rg)
{
    Org* o = makeOrg();
    for(int i = 0; i < ploidy; i++)
    {
        o->geneCounts[i] = num_G + num_P + num_SE + num_ST + num_SP ;
        o->genes[i] = malloc(sizeof(Gene*) * (num_G + num_P + num_SE + num_ST + num_SP));
        assert(o->genes[i]);
    }
    
    for(int i = 0; i < (num_G + num_P + num_SE + num_ST + num_SP); i++)
    {
        int type = 0;
        // int sensorSubType;
        if(i >= num_G && i < (num_G + num_P)) type = 1;
        else if(i >= (num_G + num_P)) type = 2;
        
        Gene* g = makeGene(type);
        g->map = gsl_rng_uniform(rg) * chromosomes + 1;
        for(int j = 0; j < max_b; j++)
        {
            g->cs[j] = gsl_ran_exponential(rg, initial_mean_cs);
            if(gsl_rng_uniform(rg) >= pos_cs) g->cs[j] *= -1;
            g->Ks[j] = (int)gsl_rng_uniform_int(rg, 5);
        }
        
        if(type == 1)
        {
            //int targets = 0;
            //while(targets == 0)
            //{
               // targets = 0;
            for(int j = 0; j < num_T; j++)
            {
                g->targets[j] = 0;
                if(gsl_rng_uniform(rg) < pleiotropy_proba)
                {
                    g->targets[j] = 1;
                    //targets++;
                }
            }
            //}
            for(int j = 0; j < num_T; j++)
            {
                g->traitEffects[j] = gsl_ran_exponential(rg, initial_sigma_trait);
            }
        }
        else
        {
            for(int j = 0; j < num_T; j++)
            {
                g->targets[j] = 0;
                g->traitEffects[j] = 0;
            }
        }

        if(dualPhenoGenes || type != 1)
        {
            g->transEffect = gsl_ran_exponential(rg, initial_mean_ts);
            if(gsl_rng_uniform(rg) < 0.5) g->transEffect *= -1;
            g->dbm = (int)(gsl_rng_uniform_int(rg, max_b));
        }
        else
        {
            g->transEffect = 0;
            g->dbm = -1;
        }
        
        if(type == 2)
        {
            // whichSignal = 0; -> Not a sensor gene
            // whichSignal = 1; -> SE
            // whichSignal = 2; -> ST
            // whichSignal = 3; -> SP
            // g->whichSignal = (int)(gsl_rng_uniform_int(rg, signals));
            g->whichSignal = 1;
            if (i >= (num_G + num_P + num_SE) && i < (num_G + num_P + num_SE + num_ST)) g->whichSignal = 2;
            else if (i >= (num_G + num_P + num_SE + num_ST )) g->whichSignal = 3;
            g->signalK = gsl_ran_exponential(rg, initial_mean_signalK);
            g->signalN = 1 + gsl_rng_uniform(rg) * 4;
        }
        else
        {
            g->whichSignal = 0;
            g->signalK = 0;
            g->signalN = 0;
        }
        
        o->genes[0][i] = g;
        for(int j = 1; j < ploidy; j++)
        {
            o->genes[j][i] = cloneGene(g);
        }
    }
    return o;
}

//Prints a text representation of a gene and also calls printBinGene()
void printGene(Gene* g, FILE* outfile, FILE* binfile, int time, int parent, int parentAllele)
{
    fprintf(outfile, "%d %d %d %d %d %d %f %d %f %d %f %f ", g->number, g->locus, parent, parentAllele, time, g->type, g->map, g->dbm, g->transEffect, g->whichSignal, g->signalK, g->signalN); 
    for(int i = 0; i < num_T; i++) fprintf(outfile, "%d %f ", g->targets[i], g->traitEffects[i]);
    
    for(int k = 0; k < max_b; k++)
    {
        fprintf(outfile, "%f ", g->cs[k]);
    }
    for(int k = 0; k < max_b; k++)
    {
        fprintf(outfile, "%d ", g->Ks[k]);
    }
    fprintf(outfile, "\n");
    fflush(outfile);
    printBinGene(g, binfile, time, parent, parentAllele);
}

//Prints a binary representation of a gene
void printBinGene(Gene* g, FILE* binfile, int time, int parent, int parentAllele)
{
    fwrite(&(g->number), sizeof(int), 1, binfile);
    fwrite(&(g->locus), sizeof(int), 1, binfile);
    fwrite(&(parent), sizeof(int), 1, binfile);
    fwrite(&(parentAllele), sizeof(int), 1, binfile);
    fwrite(&(time), sizeof(int), 1, binfile);
    fwrite(&(g->type), sizeof(int), 1, binfile);
    fwrite(&(g->map), sizeof(double), 1, binfile);
    
    fwrite(&(g->dbm), sizeof(int), 1, binfile);
    fwrite(&(g->transEffect), sizeof(double), 1, binfile);
    fwrite(&(g->whichSignal), sizeof(int), 1, binfile);
    fwrite(&(g->signalK), sizeof(double), 1, binfile);
    fwrite(&(g->signalN), sizeof(double), 1, binfile);

    fwrite(g->targets, sizeof(int), num_T, binfile);
    fwrite(g->traitEffects, sizeof(double), num_T, binfile);
    
    fwrite(g->cs, sizeof(double), max_b, binfile);
    fwrite(g->Ks, sizeof(int), max_b, binfile);
    fflush(binfile);
}

//Write an entire org to a binary file
void writeOrgToFile(Org* o, FILE* outfile)
{
    fwrite(&(o->haplotype), sizeof(int), 1, outfile);
    for(int i = 0; i < ploidy; i++)
    {
        fwrite(&(o->geneCounts[i]), sizeof(int), 1, outfile);
        for(int j = 0; j < o->geneCounts[i]; j++) printBinGene(o->genes[i][j], outfile, 0, 0, 0);
    }
}

//Prints a haplotype (comma separated list of gene ids) to a text file
void printHaplotype(Org* o, FILE* outfile, int time, int parentHaplo)
{
    char* haplotype = calloc(500, sizeof(char));
    char str[15];
    
    for(int i = 0; i < o->geneCounts[0]; i++)
    {
        sprintf(str, "%d,", o->genes[0][i]->number);
        strcat(haplotype, str);
    }
    if(o->geneCounts[0] > 0) fprintf(outfile, "%d %d %d %s\n", o->haplotype, time, parentHaplo, haplotype);
    else fprintf(outfile, "%d %d %d null\n", o->haplotype, time, parentHaplo);
    fflush(outfile);
    free(haplotype);
}

Gene* readGene(FILE* infile)
{
    int number, locus, type, dummy;
    fread(&number, sizeof(int), 1, infile);
    fread(&locus, sizeof(int), 1, infile);
    fread(&dummy, sizeof(int), 1, infile);
    fread(&dummy, sizeof(int), 1, infile);
    fread(&dummy, sizeof(int), 1, infile);
    fread(&type, sizeof(int), 1, infile);
    Gene* g = makeGene(type);
    g->number = number;
    g->locus = locus;
    fread(&(g->map), sizeof(double), 1, infile);
    
    fread(&(g->dbm), sizeof(int), 1, infile);
    fread(&(g->transEffect), sizeof(double), 1, infile);
    fread(&(g->whichSignal), sizeof(int), 1, infile);
    fread(&(g->signalK), sizeof(double), 1, infile);
    fread(&(g->signalN), sizeof(double), 1, infile);
    
    fread(g->targets, sizeof(int), num_T, infile);
    fread(g->traitEffects, sizeof(double), num_T, infile);
    
    fread(g->cs, sizeof(double), max_b, infile);
    fread(g->Ks, sizeof(int), max_b, infile);
    
    //for(int i = 0; i < max_b; i++) printf("%f ", g->cs[i]);
    //printf("\n");

    return(g);
}

//Allocates an organism and populates it with genes from the next entry in a binary file
Org* readNextOrgFromFile(FILE* infile)
{
    Org* ret = makeOrg();
    fread(&(ret->haplotype), sizeof(int), 1, infile);
    for(int i = 0; i < ploidy; i++)
    {
        fread(&(ret->geneCounts[i]), sizeof(int), 1, infile);
        ret->genes[i] = malloc(sizeof(Gene*) * ret->geneCounts[i]);
        for(int j = 0; j < ret->geneCounts[i]; j++)
        {
            //Note: repeated calls to readGene read in successive genes because the
            // variable pointed to by $infile is altered by each call to fread().
            ret->genes[i][j] = readGene(infile);
        }
    }
    return ret;
}

//Copies an organism without mutation
Org* cloneOrg(Org* o)
{
    Org* ret = makeOrg();
    ret->haplotype = o->haplotype;
    for(int i = 0; i < ploidy; i++)
    {
        ret->geneCounts[i] = o->geneCounts[i];
        ret->genes[i] = malloc(sizeof(Gene*) * ret->geneCounts[i]);
        assert(ret->genes[i]);
        for(int j = 0; j < ret->geneCounts[i]; j++)
        {
            ret->genes[i][j] = cloneGene(o->genes[i][j]);
        }
    }
    return ret;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Mutation functions: Each of these changes a value in place, rather than returning a value.
/////////////////////////////////////////////////////////////////////////////////////////////

//Mutates a cis-regulatory effect
void mutateCis(Gene* g, int site)
{
    int sign = 1;
    if(g->cs[site] < 0) sign = -1; 
    g->cs[site] = pow(2, log(fabs(g->cs[site]))/log(2) + gsl_ran_gaussian(rg, sigma_cs));
    g->cs[site] *= sign;
}

//Mutates the binding affinity for an existing site. May not change Ks due to inherent degeneracy of
// what Ks represents (i.e., # of mismatches to an ideal sequence).
void mutateSmallK(Gene* g, int site)
{
    if(gsl_rng_uniform(rg) < (double)(n_cis - g->Ks[site]) / n_cis)
    {
        g->Ks[site]++;
    }
    else if(gsl_rng_uniform(rg) < 0.33333)
    {
        g->Ks[site]--;
    }
}

//Creates an absent binding site by mutation
void mutateBigK(Gene* g, int site)
{
    g->Ks[site] = k_mut - 1;
    g->cs[site] = gsl_ran_exponential(rg, initial_mean_cs);
    if(gsl_rng_uniform(rg) >= pos_cs) g->cs[site] *= -1;
}

//Mutates the affinity of a sensor--these mutate on a continuous scale because they are not
// based on sequence matches of a binding site.
void mutateSignalK(Gene* g)
{
   g->signalK = pow(2, log(fabs(g->signalK))/log(2) + gsl_ran_gaussian(rg, sigma_signalK));
}

void mutateSignalN(Gene* g)
{
    g->signalN = gsl_rng_uniform(rg) * 10;
}

//Mutates the ligand type of a sensor. Not in use!
void mutateWhichSignal(Gene* g)
{
    printf("You called the function mutateWhichSignal!");
    /*if(signals > 1)
    {
        int oldValue = g->whichSignal;
        while(oldValue == g->whichSignal) g->whichSignal = (int)(gsl_rng_uniform_int(rg, signals));
    }*/
}

//Mutates the effect on a phenotypic trait
void mutateTrait(Gene* g)
{
    int trait = (int)(gsl_rng_uniform_int(rg, num_T));
    g->traitEffects[trait] = pow(2, log(g->traitEffects[trait]) / log(2.0) + gsl_ran_gaussian(rg, sigma_trait));
}

//Mutates the DNA binding motif of a protein--always changes the old value
void mutateDBM(Gene* g)
{
    int oldValue = g->dbm;
    while(oldValue == g->dbm) g->dbm = (int)(gsl_rng_uniform_int(rg, max_b));
}

//Mutates the trans effect of a protein
void mutateTrans(Gene* g)
{
    int sign = 1;
    if(g->transEffect < 0) sign = -1;
    g->transEffect = pow(2, log(fabs(g->transEffect)) / log(2.0) + gsl_ran_gaussian(rg, sigma_ts));
    g->transEffect *= sign;
}

//Determines which parts of a gene mutate upon replication and calls the mutation functions.
//If a gene has changed, the new allele is printed.
//Returns the number of changes to determine if a new halpotype needs to be output.
int substitutions(Gene* g, int time, int parent, gsl_rng* rg, FILE* outfile, FILE* binfile)
{
    int alleleID = g->number;
    int changed = 0;
    //Binding site mutations
    //First, draw for all c_ij effects
    int c_muts = gsl_ran_poisson(rg, mu_cs * max_b);
    while(c_muts)
    {
        int site = (int)(gsl_rng_uniform_int(rg, max_b));
        mutateCis(g, site);
        c_muts--;
        changed++;
    }
    //Second, check for mutation in each K individually
    for(int i = 0; i < max_b; i++)
    {
        if(g->Ks[i] == k_mut && gsl_rng_uniform(rg) < mu_big_K)
        {
            mutateBigK(g, i);
            changed++;
        }
        if(g->Ks[i] < k_mut && gsl_rng_uniform(rg) < mu_small_K)
        {
            mutateSmallK(g, i);
            changed++;
        }
    }
    if(g->type == 1)
    {
        //Mutations for phenotype gene
        if(gsl_rng_uniform(rg) < mu_trait)
        {
            mutateTrait(g);
            changed++;
        }
    }
    if(g->type == 2)
    {
        //Mutations for sensor gene
        if(gsl_rng_uniform(rg) < mu_signalK)
        {
            mutateSignalK(g);
            changed++;
        }
        /*if(gsl_rng_uniform(rg) < mu_whichSignal)
        {
            mutateWhichSignal(g);
            changed++;
        }*/
        /*
        if(gsl_rng_uniform(rg) < mu_signalN)
        {
            mutateSignalN(g);
            changed++;
        }
         */
    }
    if(dualPhenoGenes || g->type != 1)
    {
        //Mutations for regulatory and sensor genes
        if(gsl_rng_uniform(rg) < mu_dbm)
        {
            mutateDBM(g);
            changed++;
        }
        if(gsl_rng_uniform(rg) < mu_ts)
        {
            mutateTrans(g);
            changed++;
        }
    }
    if(changed)
    {
        g->number = next;
        if(outfile && binfile) printGene(g, outfile, binfile, time, parent, alleleID);
        next++;
    }
    return changed;
}

int ReDevelop_substitutions(Gene* g, int time, int parent, gsl_rng* rg, FILE* outfile, FILE* binfile, int nbmutations)
{
    int alleleID = g->number;
    int changed = 0;
    do
    {
        //Binding site mutations
        //First, draw for all c_ij effects
        int c_muts = gsl_ran_poisson(rg, mu_cs * max_b);
        while(c_muts)
        {
            int site = (int)(gsl_rng_uniform_int(rg, max_b));
            mutateCis(g, site);
            c_muts--;
            changed++;
        }
        //Second, check for mutation in each K individually
        for(int i = 0; i < max_b; i++)
        {
            if(g->Ks[i] == k_mut && gsl_rng_uniform(rg) < mu_big_K)
            {
                mutateBigK(g, i);
                changed++;
            }
            if(g->Ks[i] < k_mut && gsl_rng_uniform(rg) < mu_small_K)
            {
                mutateSmallK(g, i);
                changed++;
            }
        }
        if(g->type == 1)
        {
            //Mutations for phenotype gene
            if(gsl_rng_uniform(rg) < mu_trait)
            {
                mutateTrait(g);
                changed++;
            }
        }
        if(g->type == 2)
        {
            //Mutations for sensor gene
            if(gsl_rng_uniform(rg) < mu_signalK)
            {
                mutateSignalK(g);
                changed++;
            }
            /*if(gsl_rng_uniform(rg) < mu_whichSignal)
            {
                mutateWhichSignal(g);
                changed++;
            }*/
            /*
            if(gsl_rng_uniform(rg) < mu_signalN)
            {
                mutateSignalN(g);
                changed++;
            }
             */
        }
        if(dualPhenoGenes || g->type != 1)
        {
            //Mutations for regulatory and sensor genes
            if(gsl_rng_uniform(rg) < mu_dbm)
            {
                mutateDBM(g);
                changed++;
            }
            if(gsl_rng_uniform(rg) < mu_ts)
            {
                mutateTrans(g);
                changed++;
            }
        }
        if(changed)
        {
            g->number = next;
            if(outfile && binfile) printGene(g, outfile, binfile, time, parent, alleleID);
            next++;
        }
    } while (changed<nbmutations); // As long as too few mutations, just keep going. If too many then the r.wrapper will ensure to recall Redevelop_substitution. So it is important to keep an expected number of mutation lower than 1.
    return changed;
}


//Performs recombination or asexual reproduction
void makeGamete(Org* dest, int chr, Org* src, int time, gsl_rng* rg, FILE* outfile, FILE* binfile, FILE* hapFile)
{
    //Assign enough space for all genes to end up in a gamete; duplications that exceed this are handled by a realloc below
    int total = 0;
    for(int i = 0; i < ploidy; i++) total += src->geneCounts[i];
    dest->genes[chr] = malloc(sizeof(Gene*) * total);
    assert(dest->genes[chr]);
    int counter = 0;
    int changes = 0;
    dest->haplotype = src->haplotype;
    
    for(int i = 0; i < chromosomes; i++)
    {
        if(sex)
        {
            //Index of which homologous chromosome is contributing genes to the gamete
            int chosen;
            //A Poisson number of crossover points are chosen uniformly and sorted
            int crossovers = gsl_ran_poisson(rg, rec_rate);
            double* crosspoints;
            if(crossovers)
            {
                crosspoints = malloc(sizeof(double) * (crossovers + 1));
                assert(crosspoints);
                for(int j = 0; j < crossovers; j++) crosspoints[j] = gsl_rng_uniform(rg);
                crosspoints[crossovers] = 1;
                sortDouble(crosspoints, crossovers);
            }
            // A chromosome is chosen randomly and the distance from the origin to the first crossover is assigned to the interval $range
            double range [2];
            range[0] = i + 1;
            if(crossovers) range[1] = i + 1 + crosspoints[0];
            else range[1] = i + 2;
            int crossCounter = 0;
            chosen = (int)(gsl_rng_uniform_int(rg, ploidy));
            
            // Genes on the focal chromosome within $range are copied to the gamete, with substitutions and potential deletions.
            do
            {
                for(int j = 0; j < src->geneCounts[chosen]; j++)
                {
                    if(src->genes[chosen][j]->map >= range[0] && src->genes[chosen][j]->map < range[1] && gsl_rng_uniform(rg) >= mu_del)
                    {
                        dest->genes[chr][counter] = cloneGene(src->genes[chosen][j]);
                        substitutions(dest->genes[chr][counter], time, src->genes[chosen][j]->locus, rg, outfile, binfile);
                        counter++;
                    }
                }
                
                //$Range is updated to the next segment between crossovers, and the focal chromosome is switched.
                crossCounter++;
                if(crossCounter <= crossovers)
                {
                    int oldChosen = chosen;
                    while(chosen == oldChosen) chosen = (int)(gsl_rng_uniform_int(rg, ploidy));
                    range[0] = range[1];
                    range[1] = i + 1 + crosspoints[crossCounter];
                }
            }while(crossCounter <= crossovers);
            if(crossovers) free(crosspoints);
        }
        else
        {
            if(ploidy == 1)
            {
                for(int j = 0; j < src->geneCounts[0]; j++)
                {
                    if(gsl_rng_uniform(rg) >= mu_del)
                    {
                        dest->genes[chr][counter] = cloneGene(src->genes[0][j]);
                        changes += substitutions(dest->genes[chr][counter], time, src->genes[0][j]->locus, rg, outfile, binfile);
                        counter++;
                    }
                    else
                    {
                        changes++;
                    }
                }
            }
        }
    }
   
    //Duplications are possible if the offspring inherited some genes, and occur at a per-genome rate
    if(counter)
    {
        int dups = gsl_ran_poisson(rg, mu_dup);
        changes += dups;
        //Check to see if more space is needed
        if(counter + dups > total)
        {
            dest->genes[chr] = realloc(dest->genes[chr], sizeof(Gene*) * (counter + dups));
            assert(dest->genes[chr]);
        }
        
        while(dups)
        {
            int choice = (int)(gsl_rng_uniform_int(rg, counter));
            dest->genes[chr][counter] = cloneGene(dest->genes[chr][choice]);
            dest->genes[chr][counter]->map = gsl_rng_uniform(rg) * chromosomes + 1;
            dest->genes[chr][counter]->number = next;
            dest->genes[chr][counter]->locus = nextLocus;
            if(outfile && binfile) printGene(dest->genes[chr][counter], outfile, binfile, time, dest->genes[chr][choice]->locus, dest->genes[chr][choice]->number);
            next++;
            nextLocus++;
            changes += substitutions(dest->genes[chr][counter], time, dest->genes[chr][choice]->locus, rg, outfile, binfile);
            counter++;
            dups--;
        }
    }
    
    dest->geneCounts[chr] = counter;
    
    if(changes > 0)
    {
        dest->haplotype = nextHaploid;
        nextHaploid++;
        if(hapFile) printHaplotype(dest, hapFile, time, src->haplotype);
    }
}
