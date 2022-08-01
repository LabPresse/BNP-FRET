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



The code output below shows the MCMC (Markov Chain Monte Carlo) iteration number (number of samples generated), number of active system states, labels for all the active loads, the excitation rate (related to laser power), rate matrix for the biomolecule of interest with FRET efficiencies on the diagonal instead of zeros, logarithm of the full joint posterior, and acceptance rates. All the samples until the full joint posterior reaches convergence (maximum) should be ignored. This initial period is also known as burn-in. For efficient exploration of the parameter space, tune the covariance values for the proposal distributions in such a way that the acceptance rates stay between 30 and 60% approximately.

```
=================================================================
 
Iteration: 1217
 
n_system_states = 4
 
loads_active = [1, 3, 5, 6]
 
excitation rate (s^(-1)) = 4111.3703223987695
 
rate_matrix (s^(-1)) (diagonal elements show FRET efficiencies) =
 
4×4 Matrix{Float64}:
 0.0944715  0.650713    0.0648545  0.10982
 0.411844   0.709043    0.729846   1.80119
 2.97991    2.91667     0.487081   0.432123
 0.869362   0.00541285  0.700872   0.054417 
 
log_full_posterior = -13304.721733757791
 
minimum, median, maximum acceptance rate = 60.72308956450287%  61.750205423171735%  88.82497945768283%  
 
=================================================================
```

Furthemore, the sampler output can be visualized using the plotting options in the "input_parameters.jl" file. Visualization of the sampler output makes it easy to identify when the posterior converges and the most probable model (number of system states). Furthermore, it shows a bivariate distribution for FRET efficiencies and escape rates (sum of all the rates out of a state or in a rate matrix row).

![Screenshot from 2022-08-01 02-30-28](https://user-images.githubusercontent.com/87823118/182118887-f2f7426d-0508-4e8f-8bf3-dd0846466f22.png)


