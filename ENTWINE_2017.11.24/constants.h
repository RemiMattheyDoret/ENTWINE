#ifndef gene_networks_constants_h
#define gene_networks_constants_h

//Number of traits
extern int num_T;

extern const gsl_rng_type* rT;
extern gsl_rng* rg;
extern unsigned int sValue;

//*******************
//Genotype parameters
//*******************

//Initial number of regulatory genes
extern int num_G;
//Initial number of phenotype genes
extern int num_P;
//Initial number of environmental sensor genes
extern int num_SE;
//Initial number of Trait sensor genes
extern int num_ST;
//Initial number of Performance sensor genes
extern int num_SP;
//Number of chromosomes
extern int chromosomes;
//Ploidy
extern int ploidy;

extern int sex;

//Probability of a crossover per unit chromosome length
extern double rec_rate;
//Maximum number of cis-regulatory binding site possibilities
extern int max_b;
//Do we want an environmental signal?
extern int EnvSignal;

// what is the probability of any given phenotypic gene to affect any given trait? If one, every phenotypic gene will affect every trait.
extern double pleiotropy_proba;

extern int k_mut;
extern int k_weak;

extern double K_vals [4];

//***********************
//Evolutionary parameters
//***********************
extern double mu_nuc;
extern int L_cis;
extern int n_cis;
extern int n_trans;

extern double omega_squared;
extern double* trait_optima;
extern double* trait_optima_env;
extern double* trait_optima_perf;
extern double low_environment;
extern double high_environment;



//Per-active-gene mutation rate for trans effects
extern double mu_ts;
//Per-active_gene mutation rate for cis effects
extern double mu_cs;
//Per-active_gene mutation rate for affinity constants
extern double mu_small_K;
extern double mu_big_K;
//Per-active-gene deletion rate
extern double mu_del;
//Per-genome duplication rate
extern double mu_dup;
//Per-active-gene mutation rate for DNA binding type
extern double mu_dbm;
//Mutation rate for trait effects
extern double mu_trait;

extern double mu_signalK;
extern double mu_signalN;
extern double mu_whichSignal;

//Target sizes (effective number of base pairs) for mutation rates
extern double trans_mutational_target;
extern double trait_mutational_target;
extern double dbm_mutational_target;
//Note: target size for the cis effect, which is distinct from the target for its binding (n_cis).
extern double cis_mutational_target;
extern double sig_mutational_target;



//Mutation standard deviation for trans effects -- same for negative and positive
extern double sigma_ts;
//Mutation sd for trait effects
extern double sigma_trait;
//Mutation standard deviation for cis effects -- same for negative and positive, but see below
extern double sigma_cs;

extern double sigma_signalK;

//Probability of switching cis effects signs in both directions
extern double pos_cs;

//Initial distribution parameters for the mean of exponentials for G_cs and G_ts
extern double initial_mean_cs;
extern double initial_mean_ts;
extern double initial_sigma_trait;
extern double initial_mean_signalK;

//************************************
//Development parameters -- biological
//************************************
//Maximum development time
extern double max_time;
//Expression propensity of unstimulated, unrepressed gene
extern double basal;
//Expression rate coefficient
extern double transcription;
extern double translation;
//Per-capita mRNA decay rate
extern double mRNA_decay;
//Per-capita protein decay rate
extern double protein_decay;
//Can phenotype genes also regulate?
extern int dualPhenoGenes;

//Production rate for "on" signal
extern double basic_signal_rate;

//***************************************
//Development parameters -- error control
//***************************************
extern double eps;

//***************************************
//Derived constants -- set at runtime
//***************************************

extern int prefix;

extern int next;

extern int nextLocus;

extern int nextHaploid;

extern int treatment;

#endif
