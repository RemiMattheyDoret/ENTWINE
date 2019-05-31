#ifndef gene_networks_develop_h
#define gene_networks_develop_h

double* trait_optima_env;
double* trait_optima_perf;

double getFitness(Org* o);

double* develop(Org* g, int plotFlag, char* outfile, int protToCount);
double* develop2(Org* g, int plotFlag, char* outfile, int protToCount);

double* developWithTraitFeedback(Org* g, int plotFlag, char* outfile, double startingTrait);

double* developWithFeedback(Org* g, int plotFlag, char* outfile, int protToCount, double feedIn);

void developAllRNAs(Org* o);

double* developMap(Org* g, int plotFlag, char* outfile, double startingTrait);

#endif
