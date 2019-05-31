
#ifndef gene_networks_G_h
#define gene_networks_G_h

double* trait_optima_env;
double* trait_optima_perf;

typedef struct Gene
{
    //0 indicates regulatory gene, 1 indicates phenotype gene, 2 indicates sensor gene
    int type;
    // Index for DNA binding motif
    int dbm;
    // Regulatory effect
    double transEffect;
    
    //Variables used by sensor genes only
    //Flags which signal is processed by this gene 
    int whichSignal;
    //Effect on each trait
    double signalK;
    //Exponent for signal
    double signalN;
    
    //Variables used by phenotype genes only
    //Flags whether the gene affects each trait
    int* targets;
    //Effect on each trait
    double* traitEffects;
    
    //Recombination map location
    double map;
    //Locus id
    int locus;
    //Array of cis-regulatory effects
    double* cs;
    //Array of half-saturation constant indices
    int* Ks;
    
    //Gene id
    int number;
} Gene;

typedef struct Org
{
    int* geneCounts;
    Gene*** genes;
    
    double* traits;
    double fitness;
    int protCount;
    int rnaCount;
    int* allRNAs;
    int haplotype;
    
} Org;

void setup();

void cleanup();

double mutate_cs(double x, gsl_rng* rg);

int mutate_Ks(int x, gsl_rng* rg);

Gene* makeGene(int type);

void deleteGene(Gene* g);

Org* makeOrg();

Org* makeHomozygote(Org* src, int side);

void deleteOrg(Org* o);

Gene* cloneGene(Gene* g);

Org* randomOrg(gsl_rng* rg);

void printGene(Gene* g, FILE* outfile, FILE* binfile, int time, int parent, int parentAllele);

void printBinGene(Gene* g, FILE* binfile, int time, int parent, int parentAllele);

void writeOrgToFile(Org* o, FILE* outfile);

void printHaplotype(Org* o, FILE* outfile, int time, int parentHaplo);

Gene* readGene(FILE* infile);

Org* readNextOrgFromFile(FILE* infile);

Org* cloneOrg(Org* o);

void recordChange(Gene* g, int time, int mT, int sL);

void mutateCis(Gene* g, int site);

void mutateSmallK(Gene* g, int site);

void mutateBigK(Gene* g, int site);

void mutateSignalK(Gene* g);

void mutateSignalN(Gene* g);

void mutateWhichSignal(Gene* g);

void mutateTrait(Gene* g);

void mutateDBM(Gene* g);

void mutateTrans(Gene* g);

int substitutions(Gene* g, int time, int parent, gsl_rng* rg, FILE* outfile, FILE* binfile);

int ReDevelop_substitutions(Gene* g, int time, int parent, gsl_rng* rg, FILE* outfile, FILE* binfile, int nbmutations);

void makeGamete(Org* dest, int chr, Org* src, int time, gsl_rng* rg, FILE* outfile, FILE* binfile, FILE* hapFile);

#endif
