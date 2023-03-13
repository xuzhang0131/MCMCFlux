# MCMCFlux
### This repository is served as the source code for Bayesian kinetic modeling for tracer-based metabolomic data.
### myMCMC.cpp is the main function.
### genDRAM_GibbsWithMH.cpp is for the MCMC sampling.
### PurineSynthesis.cpp provides all the ODE equations.
### To compile the cpp functions, run "g++ -o myMCMC ./myMCMC.cpp data_process.cpp genDRAM_GibbsWithMH.cpp PurineSynthesis.cpp"
### final.R is the R code for the summarization of the MCMC samples.
