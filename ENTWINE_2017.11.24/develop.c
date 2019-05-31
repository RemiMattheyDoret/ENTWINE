/*
    Contains develop() and variants, which handle stochastic production of phenotypes.
    Also contains a Poisson RNG optimized for speed.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "G.h"
#include "constants.h"

// Slightly modified copy of rpois() function from R.
// Comments from original.

#define a0	-0.5
#define a1	 0.3333333
#define a2	-0.2500068
#define a3	 0.2000118
#define a4	-0.1661269
#define a5	 0.1421878
#define a6	-0.1384794
#define a7	 0.1250060

#define one_7	0.1428571428571428571
#define one_12	0.0833333333333333333
#define one_24	0.0416666666666666667
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */

#define repeat for(;;)
#define ISNAN(x) (isnan(x)!=0)
#define FALSE (0)
#define TRUE (1)


int imax2(int a, int b)
{
    if(a > b) return a;
    else return b;
}

int imin2(int a, int b)
{
    if(a < b) return a;
    else return b;
}

double fsign(double x, double y)
{
    if (ISNAN(x) || ISNAN(y)) return x + y;
    return ((y >= 0) ? fabs(x) : -fabs(x));
}


double rpois(double mu)
{
    /* Factorial Table (0:9)! */
    const static double fact[10] =
    {
        1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };
    
    /* These are static --- persistent between calls for same mu : */
    static int l, m;
    
    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold	 = 0.*/
    
    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = FALSE;
    
    //if (!R_FINITE(mu) || mu < 0)
    //    ML_ERR_return_NAN;
    
    if (mu <= 0.)
        return 0.;
    
    big_mu = mu >= 10.;
    if(big_mu)
        new_big_mu = FALSE;
    
    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */
        
        if (big_mu) {
            new_big_mu = TRUE;
            /* Case A. (recalculation of s,d,l	because mu has changed):
             * The poisson probabilities pk exceed the discrete normal
             * probabilities fk whenever k >= m(mu).
             */
            muprev = mu;
            s = sqrt(mu);
            d = 6. * mu * mu;
            big_l = floor(mu - 1.1484);
            /* = an upper bound to m(mu) for all mu >= 10.*/
        }
        else { /* Small mu ( < 10) -- not using normal approx. */
            
            /* Case B. (start new table and calculate p0 if necessary) */
            
            /*muprev = 0.;-* such that next time, mu != muprev ..*/
            if (mu != muprev) {
                muprev = mu;
                m = imax2(1, (int) mu);
                l = 0; /* pp[] is already ok up to pp[l] */
                q = p0 = p = exp(-mu);
            }
            
            repeat {
                /* Step U. uniform sample for inversion method */
                u = gsl_rng_uniform(rg);
                if (u <= p0)
                    return 0.;
                
                /* Step T. table comparison until the end pp[l] of the
                 pp-table of cumulative poisson probabilities
                 (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
                if (l != 0) {
                    for (k = (u <= 0.458) ? 1 : imin2(l, m);  k <= l; k++)
                        if (u <= pp[k])
                            return (double)k;
                    if (l == 35) /* u > pp[35] */
                        continue;
                }
                /* Step C. creation of new poisson
                 probabilities p[l..] and their cumulatives q =: pp[k] */
                l++;
                for (k = l; k <= 35; k++) {
                    p *= mu / k;
                    q += p;
                    pp[k] = q;
                    if (u <= q) {
                        l = k;
                        return (double)k;
                    }
                }
                l = 35;
            } /* end(repeat) */
        }/* mu < 10 */
        
    } /* end {initialize persistent vars} */
    
    /* Only if mu >= 10 : ----------------------- */
    
    /* Step N. normal sample */
    g = mu + s * gsl_ran_gaussian(rg, 1);/* norm_rand() ~ N(0,1), standard normal */
    
    if (g >= 0.) {
        pois = floor(g);
        /* Step I. immediate acceptance if pois is large enough */
        if (pois >= big_l)
            return pois;
        /* Step S. squeeze acceptance */
        fk = pois;
        difmuk = mu - fk;
        u = gsl_rng_uniform(rg); /* ~ U(0,1) - sample */
        if (d * u >= difmuk * difmuk * difmuk)
            return pois;
    }
    
    /* Step P. preparations for steps Q and H.
     (recalculations of parameters if necessary) */
    
    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
         because one might have exited in step I or S
         */
        muprev2 = mu;
        omega = M_1_SQRT_2PI / s;
        /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
         * approximations to the discrete normal probabilities fk. */
        
        b1 = one_24 / mu;
        b2 = 0.3 * b1 * b1;
        c3 = one_7 * b1 * b2;
        c2 = b2 - 15. * c3;
        c1 = b1 - 6. * b2 + 45. * c3;
        c0 = 1. - b1 + 3. * b2 - 15. * c3;
        c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }
    
    if (g >= 0.) {
        /* 'Subroutine' F is called (kflag=0 for correct return) */
        kflag = 0;
        goto Step_F;
    }
    
    
    repeat {
        /* Step E. Exponential Sample */
        
        E = gsl_ran_exponential(rg, 1.0);	/* ~ Exp(1) (standard exponential) */
        
        /*  sample t from the laplace 'hat'
         (if t <= -0.6744 then pk < fk for all mu >= 10.) */
        u = 2 * gsl_rng_uniform(rg) - 1.;
        
        t = 1.8 + fsign(E, u);
        if (t > -0.6744) {
            pois = floor(mu + s * t);
            fk = pois;
            difmuk = mu - fk;
            
            /* 'subroutine' F is called (kflag=1 for correct return) */
            kflag = 1;
            
        Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */
            
            if (pois < 10) { /* use factorials from table fact[] */
                px = -mu;
                py = pow(mu, pois) / fact[(int)pois];
            }
            else {
                /* Case pois >= 10 uses polynomial approximation
                 a0-a7 for accuracy when advisable */
                del = one_12 / fk;
                del = del * (1. - 4.8 * del * del);
                v = difmuk / fk;
                if (fabs(v) <= 0.25)
                    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                                          v + a3) * v + a2) * v + a1) * v + a0)
                    - del;
                else /* |v| > 1/4 */
                    px = fk * log(1. + v) - difmuk - del;
                py = M_1_SQRT_2PI / sqrt(fk);
            }
            x = (0.5 - difmuk) / s;
            x *= x;/* x^2 */
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
            if (kflag > 0) {
                /* Step H. Hat acceptance (E is repeated on rejection) */
                if (c * fabs(u) <= py * exp(px + E) - fy * exp(fx + E))
                    break;
            } else
            /* Step Q. Quotient acceptance (rare case) */
                if (fy - u * fy <= py * exp(px - fx))
                    break;
        }/* t > -.67.. */
    }
    return pois;
}


double getFitness(Org* o)
{
    return o->fitness;
}

void develop(Org* o, int plotFlag, char* outfile, int protToCount)
{
    //printf("haplotype = %d\n",o->haplotype);
    // If plotFlag > 1, write protein abundances to outfile
    // If protToCount != -1, write total expression of protein $protToCount to Org->protCount
    o->protCount = 0;
    o->rnaCount = 0;
    FILE* oFile;
    if(plotFlag) oFile = fopen(outfile, "a");
    o->fitness = 0;
    double sigRate = basic_signal_rate;
    // int signalsToTrack = signals;
    // if(traitFeedback) signalsToTrack--;
    //double* signalProducts = malloc(sizeof(double) * signalsToTrack);
    //double* signalRates = malloc(sizeof(double) * signalsToTrack);
    double signalRates = 0;
    double signalProducts = 0;
    if(EnvSignal == 1)
    {
        signalRates = trait_optima_env[0] * protein_decay;
        signalProducts = 0;
    }
    double diff;
    int counter;
    for(int i = 0; i < num_T; i++) o->traits[i] = 0;
    int basic_signal = 0;
    //Determine total number of genes
    int total = 0;
    for(int i = 0; i < ploidy; i++)
    {
        total += o->geneCounts[i];
    }
    //Vectors for sensor protein information
    int* whichSignals = malloc(sizeof(int) * total);
    double* signalKS = malloc(sizeof(double) * total);
    double* signalNS = malloc(sizeof(double) * total);
    //Vector of pointers to phenotype genes from protein index; pointer is null for regulatory and sensor genes
    Gene** pGenes = malloc(sizeof(Gene*) * total);
    assert(pGenes);
    counter = 0;
    for(int j = 0; j < ploidy; j++)
    {
        for(int k = 0; k < o->geneCounts[j]; k++)
        {
            if(o->genes[j][k]->type == 1)
            {
                pGenes[counter] = o->genes[j][k];
            }
            else
            {
                pGenes[counter] = 0;
            }
            if(o->genes[j][k]->type == 2)
            {
                whichSignals[counter] = o->genes[j][k]->whichSignal;
                signalKS[counter] = o->genes[j][k]->signalK;
                signalNS[counter] = o->genes[j][k]->signalN;
            }
            else
            {
                whichSignals[counter] = -1;
                signalKS[counter] = 0;
                signalNS[counter] = 0;
            }
            counter++;
        }
    }
    //Set up a lookup table for reading out which protein products go with which binding sites, as well as their trans effects
    int* counts = malloc(sizeof(int) * max_b);
    int** links = malloc(sizeof(int*) * max_b);
    double** trans = malloc(sizeof(double*) * max_b);
    int* tmpProt = malloc(sizeof(int) * total);
    double* tmpTrans = malloc(sizeof(double) * total);
    assert(counts && links && trans && tmpProt && tmpTrans);
    for(int i = 0; i < max_b; i++)
    {
        //For each type of binding site, count the number of interacting proteins ($counts), their indices ($links) and trans effects ($trans)
        //Use a temporary set of vectors since we don't know how a priori how many proteins interact with each site
        counts[i] = 0;
        int proteinCounter = 0;
        for(int j = 0; j < ploidy; j++)
        {
            for(int k = 0; k < o->geneCounts[j]; k++)
            {
                if((dualPhenoGenes || o->genes[j][k]->type != 1) && o->genes[j][k]->dbm == i)
                {
                    tmpProt[counts[i]] = proteinCounter;
                    tmpTrans[counts[i]] = o->genes[j][k]->transEffect;
                    counts[i]++;
                }
                proteinCounter++;
            }
        }
        if(counts[i])
        {
            //Now make vectors of the appropriate length and copy over
            links[i] = malloc(sizeof(int) * counts[i]);
            trans[i] = malloc(sizeof(double) * counts[i]);
            assert(links[i] && trans[i]);
            for(int j = 0; j < counts[i]; j++)
            {
                links[i][j] = tmpProt[j];
                trans[i][j] = tmpTrans[j];
            }
        }
        else
        {
            links[i] = NULL;
        }
    }
    free(tmpProt);
    free(tmpTrans);
    // 0 = mRNA production
    // 1 = protein production
    // 2 = mRNA decay
    // 3 = protein decay
    //Keep track of old expression rates to set time-steps.
    double* rates = malloc(sizeof(double) * total);
    double* oldRates = malloc(sizeof(double) * total);
    //Store the calcuated acceptable tau with respect to each reaction
    double* taus = malloc(sizeof(double) * total);
    double* deltaTrait = malloc(sizeof(double) * num_T);
    assert(rates && oldRates && taus);
    //Vectors for all gene products
    int* products = malloc(sizeof(int) * total);
    int* transcripts = malloc(sizeof(int) * total);
    double* eProducts = malloc(sizeof(double) * total);
    // All products start at zero
    for(int i = 0; i < total; i++)
    {
        products[i] = 0;
        transcripts[i] = 0;
    }
    //Keep track of old deltaT to set a proportional max rate of change
    double deltaT = 0.1;
    double oldDeltaT = deltaT;
    
    double t = 0;
    //Look up table for significant binding sites
    int*** lut_ID = malloc(sizeof(int**) * ploidy);
    int** lut_size = malloc(sizeof(int*) * ploidy);
    assert(lut_ID && lut_size);
    for(int i = 0; i < ploidy; i++)
    {
        lut_ID[i] = malloc(sizeof(int*) * o->geneCounts[i]);
        lut_size[i] = malloc(sizeof(int) * o->geneCounts[i]);
        assert(lut_ID[i] && lut_size[i]);
        for(int j = 0; j < o->geneCounts[i]; j++)
        {
            lut_ID[i][j] = malloc(sizeof(int) * max_b);
            assert(lut_ID[i][j]);
            lut_size[i][j] = 0;
            for(int k = 0; k < max_b; k++)
            {
                if((counts[k] || k == 0) && o->genes[i][j]->Ks[k] < k_weak)
                {
                    lut_ID[i][j][lut_size[i][j]] = k;
                    lut_size[i][j]++;
                }
            }
        }
    }
    //Look up tables to access protein effects in convenient form:
    // i == chromosome, j == gene, b == binding site, p = protein
    double**** cts;
    int**** signs;
    cts = malloc(sizeof(double***) * ploidy);
    signs = malloc(sizeof(int***) * ploidy);
    assert(cts && signs);
    for(int i = 0; i < ploidy; i++)
    {
        cts[i] = malloc(sizeof(double**) * o->geneCounts[i]);
        signs[i] = malloc(sizeof(int**) * o->geneCounts[i]);
        assert(cts[i] && signs[i]);
        for(int j = 0; j < o->geneCounts[i]; j++)
        {
            cts[i][j] = malloc(sizeof(double*) * lut_size[i][j]);
            signs[i][j] = malloc(sizeof(int*) * lut_size[i][j]);
            assert(cts[i][j] && signs[i][j]);
            for(int b = 0; b < lut_size[i][j]; b++)
            {
                int k = lut_ID[i][j][b];
                cts[i][j][b] = malloc(sizeof(double) * counts[k]);
                signs[i][j][b] = malloc(sizeof(int) * counts[k]);
                assert(cts[i][j][b] && signs[i][j][b]);
                for(int p = 0; p < counts[k]; p++)
                {
                    signs[i][j][b][p] = 1;
                    if(o->genes[i][j]->cs[k] * trans[k][p] < 0) signs[i][j][b][p] = -1;
                    cts[i][j][b][p] = 1- exp(- fabs(o->genes[i][j]->cs[k] * trans[k][p]));
                }
            }
        }
    }
    while(t < max_time)
    {
        counter = 0;
        
        // Make vector of effective protein abundances
        for(int i = 0; i < total; i++)
        {
            eProducts[i] = products[i];
            if(whichSignals[i] >= 1)
            {
                if(whichSignals[i] == 1) // SE
                {
                    eProducts[i] = products[i] * pow(signalProducts,signalNS[i]) / (pow(signalKS[i],signalNS[i]) + pow(signalProducts,signalNS[i]));
                }
                else if(whichSignals[i] == 2) // ST Works only with one (the first one) trait for the moment.
                {
                    eProducts[i] = products[i] * pow(o->traits[0],signalNS[i]) / (pow(signalKS[i],signalNS[i]) + pow(o->traits[0],signalNS[i]));
                }
                else if(whichSignals[i] == 3) // SP
                {
                    eProducts[i] = products[i] * pow(fabs(o->traits[0] - trait_optima_perf[0]),signalNS[i]) / (pow(signalKS[i],signalNS[i]) + pow(fabs(o->traits[0] - trait_optima_perf[0]),signalNS[i]));
                }
            }
        }
        for(int i = 0; i < ploidy; i++)
        {
            for(int j = 0; j < o->geneCounts[i]; j++)
            {
                double enhancing = 1 - basal;
                double repressing = 1;
                //lut_size[i][j] == how many significant binding sites gene ij has
                for(int b = 0; b < lut_size[i][j]; b++)
                {
                    //k == which protein product affects binding site b
                    int k = lut_ID[i][j][b];
                    //if more than one protein affects binding site b of gene ij
                    // OR if k == 0 (signal-binding site)

                    double totalProtein = 0;
                    for(int p = 0; p < counts[k]; p++) totalProtein += eProducts[links[k][p]];
                    if(k == 0 && o->genes[i][j]->Ks[k] < k_weak)
                    {
                        totalProtein += basic_signal;
                        if(o->genes[i][j]->cs[k] > 0) enhancing *= 1 - (1 - exp(-o->genes[i][j]->cs[k])) * basic_signal / (totalProtein + K_vals[o->genes[i][j]->Ks[k]]);
                        else repressing *= 1 - (1 - exp(o->genes[i][j]->cs[k])) * basic_signal / (totalProtein + K_vals[o->genes[i][j]->Ks[k]]);
                    }
                    if(counts[k] && totalProtein)
                    {
                        for(int p = 0; p < counts[k]; p++)
                        {
                            if(signs[i][j][b][p] == 1)
                            {
                                enhancing *= 1 - cts[i][j][b][p] * (double)(eProducts[links[k][p]]) / (totalProtein + K_vals[o->genes[i][j]->Ks[k]]);
                            }
                            else
                            {
                                repressing *= 1 - cts[i][j][b][p] * (double)(eProducts[links[k][p]]) / (totalProtein + K_vals[o->genes[i][j]->Ks[k]]);
                            }
                        }
                    }
                }
                assert(repressing <= 1 && repressing >= 0);
                assert(enhancing <= 1 && enhancing >= 0);
                rates[counter] = transcription * (1 - enhancing) * repressing;
                counter++;
            }
        }
        if(t == 0)
        {
            deltaT = 1;
        }
        else
        {
            //Find largest tolerable deltaT
            for(int i = 0; i < total; i++)
            {
                if((diff = fabs(oldRates[i] - rates[i])))
                {
                    taus[i] = eps * oldDeltaT * oldRates[i] / diff;
                }
                else taus[i] = 5;
                //Note: this safeguard obviates the binomial approximation for mRNAs, below.
                if(taus[i] >= 5) taus[i] = 5;
                if((i == 0) || taus[i] < deltaT) deltaT = taus[i];
            }
            
            //Check for slowdown condition (protein levels far from equilibrium)
            for(int i = 0; i < total; i++)
            {
                if( fabs(transcripts[i] * translation - products[i] * protein_decay) / (1 + products[i]) > 0.25)
                {
                    deltaT /= 2;
                    break;
                }
            }
        }
        //////////////////////////////////////////////////////////////////
        //Change abundances of mRNAs and proteins
        //Use binomial processes if the expected decay is greater than 10%
        //Correct binomials for decay during deltaT
        //////////////////////////////////////////////////////////////////
        if(deltaT > 1.25)
        {
            basic_signal -= gsl_ran_binomial(rg, 1 - exp(-protein_decay * deltaT), basic_signal);
            basic_signal += gsl_ran_binomial(rg, (1 - exp(-protein_decay * deltaT)) / (protein_decay * deltaT), rpois(sigRate * deltaT));
        }
        else
        {
            basic_signal -= rpois(protein_decay * deltaT * basic_signal);
            basic_signal += rpois(sigRate * deltaT);
        }
        
        if(basic_signal < 0) basic_signal = 0;
        
    
        signalProducts += (1 - exp(-protein_decay * deltaT)) * (signalRates / protein_decay - signalProducts);
        if(signalProducts < 0) signalProducts = 0;
        
        
        for(int i = 0; i < num_T; i++) deltaTrait[i] = 0;
        for(int i = 0; i < total; i++)
        {
            int deltaR = 0;
            if(deltaT > 6)
            {
                deltaR -= gsl_ran_binomial(rg, 1 - exp(-mRNA_decay * deltaT), transcripts[i]);
                int transcription = rpois(deltaT * rates[i]);
                deltaR += gsl_ran_binomial(rg, (1 - exp(-mRNA_decay * deltaT)) / (mRNA_decay * deltaT), transcription);
                if(i == protToCount) o->rnaCount += transcription;
            }
            else
            {
                deltaR -= rpois(mRNA_decay * deltaT * transcripts[i]);
                int transcription = rpois(deltaT * rates[i]);
                deltaR += transcription;
                if(i == protToCount) o->rnaCount += transcription;
            }
            
            int deltaP = 0;
            if(deltaT > 1.25)
            {
                deltaP -= gsl_ran_binomial(rg, 1 - exp(-protein_decay * deltaT), products[i]);
                int production = rpois(deltaT * translation * transcripts[i]);
                deltaP += gsl_ran_binomial(rg, (1 - exp(-protein_decay * deltaT)) / (protein_decay * deltaT), production);
                if(i == protToCount) o->protCount += production;
            }
            else
            {
                deltaP -= rpois(protein_decay * deltaT * products[i]);
                int production = rpois(deltaT * translation * transcripts[i]);
                deltaP += production;
                if(i == protToCount) o->protCount += production;
            }
            if(deltaP + products[i] < 0) deltaP = -products[i];
            if(deltaR + transcripts[i] < 0) deltaR = -transcripts[i];
            if(pGenes[i])
            {
                for(int j = 0; j < num_T; j++)
                {
                    if(pGenes[i]->targets[j]) deltaTrait[j] += pGenes[i]->traitEffects[j] * (products[i] + 0.5f * deltaP) * deltaT;
                }
            }
            products[i] += deltaP;
            transcripts[i] += deltaR;
            assert(products[i] >= 0 && transcripts[i] >= 0);
        }
        for(int i = 0; i < num_T; i++) o->traits[i] += deltaTrait[i];
        
        t += deltaT;
        double* temp = rates;
        rates = oldRates;
        oldRates = temp;
        oldDeltaT = deltaT;
    }
    double SSE = 0;
    for(int i = 0; i < num_T; i++)
    {
        SSE += pow(trait_optima[i] - o->traits[i], 2);
    }
    
    o->fitness = 0.01 + 0.99 * exp( -SSE / omega_squared);
    
    for(int i = 0; i < ploidy; i++)
    {
        for(int j = 0; j < o->geneCounts[i]; j++)
        {
            for(int b = 0; b < lut_size[i][j]; b++)
            {
                free(cts[i][j][b]);
                free(signs[i][j][b]);
            }
            free(cts[i][j]);
            free(signs[i][j]);
        }
        free(cts[i]);
        free(signs[i]);
    }
    free(cts);
    free(signs);
    for(int i = 0; i < ploidy; i++)
    {
        for(int j = 0; j < o->geneCounts[i]; j++)
        {
            free(lut_ID[i][j]);
        }
        free(lut_ID[i]);
        free(lut_size[i]);
    }
    free(lut_ID);
    free(lut_size);
    for(int i = 0; i < max_b; i++)
    {
        if(counts[i])
        {
            free(links[i]);
            free(trans[i]);
        }
    }
    free(whichSignals);
    free(signalKS);
    free(counts);
    free(links);
    free(trans);
    free(rates);
    free(oldRates);
    free(taus);
    free(pGenes);
    free(products);
    free(transcripts);
    free(deltaTrait);
    free(eProducts);
    //free(signalProducts);
    //free(signalRates);
    free(signalNS);
    if(plotFlag) fclose(oFile);
}


