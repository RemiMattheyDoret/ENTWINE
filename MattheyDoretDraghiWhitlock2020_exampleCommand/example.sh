# This example is for the simulation constant environment with the environmental signal treatment. 

/home/matthey/Plasticity/ENTWINE_2017.11.24/makePopulation -seed=4527123 -N=10000 -GENS=5005 -num_T=1 -num_G=3 -num_P=3 -num_SE=3 -chromosomes=1 -ploidy=1 -sex=0 -rec_rate=1.000000000000 -max_b=20 -EnvSignal=1 -k_mut=4 -k_weak=4 -mu_nuc=0.00000001f -L_cis=500 -n_cis=12 -low_environment=1000. -high_environment=3000. -mu_del=0.000010 -mu_dup=0.000001 -trans_mutational_target=100 -trait_mutational_target=100 -dbm_mutational_target=100 -cis_mutational_target=100 -sig_mutational_target=100 -sigma_ts=0.500000 -sigma_trait=0.500000 -sigma_cs=0.500000 -sigma_signalK=0.500000 -pos_cs=0.900000 -initial_mean_cs=0.500000 -initial_mean_ts=0.500000 -initial_sigma_trait=0.100000 -initial_mean_signalK=2000 -max_time=300.000000 -basal=0.000000 -transcription=0.680000 -translation=3.900000 -mRNA_decay=0.016700 -protein_decay=0.080000 -eps=0.100000 -omega_squared=500000.000000 -prefix=190 -dual_pheno_genes=0 -project=constEnv_envSignal -restarting=0 -num_SP=0 -num_ST=0 -pleiotropy_proba=1 -home=/global/scratch/matthey/Plasticity/outputs/ -EnvHeterogeneity=1 -expressInKiloGenerations=0