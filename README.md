# BNP-FRET: Analyze single photon smFRET data in a Bayesian nonparametrics (BNP) paradigm.

BNP-FRET is a suite of software tools to analyze single photon smFRET data under continuous and pulsed illumination. It implements Markov chain Monte Carlo (MCMC) algorithms to learn distributions over parameters of interest: the number of states a biomolecule transitions through and associated transition rates. The tools can be used in a simple plug and play manner. Check the following set of papers to see details of all the mathematics involved in the development of BNP-FRET:

https://biorxiv.org/cgi/content/short/2022.07.20.500887v1

https://biorxiv.org/cgi/content/short/2022.07.20.500888v1

https://biorxiv.org/cgi/content/short/2022.07.20.500892v1

All the codes are written in Julia language for high performance/speed and its open-source/free availability. Julia also allows easy parallelization of all the codes. To install julia, please download and install julia language from their official website (see below) for your operating system or use your package manager. The current version of the code has been successfully tested on Ubuntu 20.04, macOS 12, and Windows ..

https://julialang.org/

Once the julia language has been installed, some essential julia packages are required to be added that help simplify linear algebra and statistical calculations, and plotting. To add these package via julia REPL, first enter the julia package manager by executing "]" command in the REPL. Then simply execute the following command to install all these packages at the same time. 

```add Distributions LinearAlgebra Statistics Plots StatsPlots KernelDensity HDF5```

Also, see the image below for an example of the package installation process in julia REPL.

![Screenshot from 2022-08-01 13-00-40](https://user-images.githubusercontent.com/87823118/182234995-db174ea5-3157-4b8c-98b9-dd0aeabc4399.png)


Each specialization of BNP-FRET (continuous or pulsed illumination) is organized in such a way that all the user input is accomplished via the "input_parameters.jl" file. It can be used to provide file names for experimental FRET data and sampler output, background rates for each detection channel, crosstalk probabilities, and plotting options. See the respective files for more details (they are well-commented).

The functions used to perform all the computations are organized in a hierarchical structure. The file "sampler_continuous_illumination.jl" contains the main sampler. All the functions called in that file are written in the file "functions_layer_1.jl". Similarly, all the functions called in the file "functions_layer_1.jl" are in file "functions_layer_2.jl". Any new set of functions can be introduced by following the same heirarchical structure. A brief description of all the functions is given below in the sequence they are called:

1. get_FRET_data(): Used to obtain photon arrival data and corresponding detection channels from input files in HDF5 format. It can be easily modified if other file formats are desired.

2. sampler(): Used to generate samples for parameters of interest using Gibbs algorithm.

3. check_for_existing_mcmc_data(): Called by sampler() to searche for previously generate MCMC samples stored in HDF5 format files in the working directory.

4. initialize_params(): Called by sampler() to initialize all the parameters of interest to constitute the first set of MCMC samples.

5. get_log_likelihood(): Called by sampler() to computes the logarithm of the likelihood function for the FRET data.

   get_generator(): Called by get_log_likelihood() to obtain the full generator matrix containing photophysical and biomolecular transition rates.
   
   get_rho(): Called by get_log_likelihood() to obtain the initial probability vector associated with the generator matrix.
   
   non_radiative_propagator: Called by get_log_likelihood() to compute propagators during periods when no photons are detected.
   
   radiative_propagator: Called by get_log_likelihood() to compute propagators at photon arrival times.

6. get_log_prior_rates(): Called by sampler() to obtain logarithm of prior density at a parameter's value.

7. get_log_full_posterior(): Called by sampler() to get full joint posterior. Sums logarithms of likelihood and all the priors.

8. save_mcmc_data(): Called by sampler() to save output to files.

9. print_and_plotting(): Called by sampler() to print values on standard out (terminal/screen). It also generates plots (shown in the example below).

10. propose_params(): Called by sampler() in the for loop to propose new samples for each parameter that may or may not be accepted in the Metropolis-Hastings (MH) step.



An example of code output below shows the MCMC (Markov Chain Monte Carlo) iteration number (number of samples generated), number of active system states, labels for all the active loads, the excitation rate (related to laser power), rate matrix for the biomolecule of interest with FRET efficiencies on the diagonal instead of zeros, logarithm of the full joint posterior, and acceptance rates. We used smFRET data from experiments studying binding and unbinding of intrinsically disordered proteins (immobilized ACTR and freely-diffusing NCBD) in presence of 36 % ethylene glycol to generate this output.

https://www.pnas.org/doi/full/10.1073/pnas.1921617117


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






Furthemore, the sampler output can be visualized using the plotting options in the "input_parameters.jl" file. Visualization of the sampler output makes it easy to identify when the posterior converges and the most probable model (number of system states). Furthermore, it shows a bivariate distribution for FRET efficiencies and escape rates (sum of all the rates out of a state or in a rate matrix row). All the samples until the full joint posterior reaches convergence (maximum value) should be ignored. This initial period is also known as burn-in. For efficient exploration of the parameter space, tune the covariance values for the proposal distributions in such a way that the acceptance rates stay between 30 and 60% approximately.


![Screenshot from 2022-08-01 02-30-28](https://user-images.githubusercontent.com/87823118/182118887-f2f7426d-0508-4e8f-8bf3-dd0846466f22.png)


