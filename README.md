# BNP-FRET: Analyze single photon smFRET data in a Bayesian nonparametrics (BNP) paradigm.

BNP-FRET is a suite of software tools to analyze single photon smFRET data under continuous and pulsed illumination. It helps learn distributions over parameters of interest: the number of states a biomolecule transitions through and associated transition rates. The tools can be used in a simple plug and play manner. Check the following set of papers to see details of all the mathematics involved in the development of BNP-FRET:

https://biorxiv.org/cgi/content/short/2022.07.20.500887v1

https://biorxiv.org/cgi/content/short/2022.07.20.500888v1

https://biorxiv.org/cgi/content/short/2022.07.20.500892v1

All the codes are written in Julia language for high performance and its open-source/free availability. Julia also allows easy parallelization of all the codes.

Each specialization of BNP-FRET (continuous or pulsed illumination) is organized in such a way that all the user input is accomplished via the "input_parameters.jl" file. It can be used to provide file names for experimental FRET data and sampler output, background rates for each detection channel, crosstalk probabilities, and plotting options. See the respective files for more details (they are well-commented).

The functions used to perform all the computations are organized in a hierarchical structure. The file "sampler_continuous_illumination.jl" contains the main sampler. All the functions called in that file are written in the file "functions_layer_1.jl". Similarly, all the functions called in the file "functions_layer_1.jl" are in file "functions_layer_2.jl". Any new set of functions can be introduced by following the same heirarchical structure. A brief description of all the functions is given below in the sequence they are called:

1. get_FRET_data(): Used to obtain photon arrival data and corresponding detection channels from input files in HDF5 format. It can be easily modified if other file formats are desired.

2. sampler():

3. check_for_existing_mcmc_data():

4. initialize_params():

5. log_likelihood():

   get_generator():
   
   get_rho():
   
   non_radiative_propagator:
   
   radiative_propagator:

6. log_prior_rates():

7. get_log_full_posterior():

8. save_mcmc_data():

9. print_and_plotting():

10. propose_params():
