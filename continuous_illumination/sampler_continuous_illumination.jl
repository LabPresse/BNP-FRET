# BNP-FRET (Continuous Illumination)

###########################################################
# Copyright (C) 2022 Presse Lab - All Rights Reserved
# 
# Author: Ayush Saurabh
#
# You may use, distribute and modify this code under the
# terms of the MIT license, which unfortunately won't be
# written for another century.
###########################################################

using Random, Distributions
using LinearAlgebra
using Statistics
using Plots, Plots.PlotMeasures
using StatsPlots, KernelDensity
using HDF5


include("input_parameters.jl")
include("functions_layer_1.jl")
include("functions_layer_2.jl")

# Globally defined FRET data.
photon_arrival_times, detection_channels = get_FRET_data()

# Sampler utilizing Gibbs algorithm with Metropolis-Hastings steps
function sampler(total_draws)

 	# Create arrays storing all the generated MCMC samples, log posterior
	# values, and acceptance rates.
	mcmc_samples = zeros(n_parameters+1, total_draws)
	# +1 above is for storing the "n_system_states".
	mcmc_log_posteriors = zeros(n_parameters+1, total_draws)
	# +1 above for storing joint posterior.
	mcmc_acceptance_rates = zeros(n_parameters, total_draws)

	# Since MCMC chains can be restarted at any time, we define the variable
	# last_draw to tell the sampler how many samples have been collected before.
	# We initialized this variable here.
	last_draw = check_for_existing_mcmc_data(mcmc_samples, mcmc_log_posteriors,
													mcmc_acceptance_rates)

	# Initialize old_parameters array for the sampler. This array contains the
	# values for all the parameters of interest at the current draw in the MCMC
	# chain. The excitation rate, donor relaxation rate, and acceptor relaxation
	# rate constitute the first three elements of this array. Next, to store the
	# FRET rates and system transitions rates, we first replace the diagonal
	# elements of the system transition rate matrix with FRET rates associated
	# with the specific rows' system states. Then we store these parameters in
	# the old_parameters array row-by-row.
	#
	# In summary, with max_n_system_states be the maximum number of system
	# states used for nonparametrics, we have:
	#
	# old_parameters[1] = excitation rate,
	# old_parameters[2] = donor relaxation rate,
	# old_parameteres[3] = acceptor relaxation rates, and
	#
	# with ij = (i-1)*max_n_system + j, where i represents row and j represents
	# column of the rate matrix,
	#
	# old_parameters[3 + ij] = FRET rates, if i = j,
	# old_parameters[3 + ij] = system transition rates for system state "i",
	# if i != j.
	#
	# All the parameters are stored this way to keep the structure of the
	# sampler as general as possible.
	#
	if last_draw == 1

		old_parameters = initialize_params()

		# Compute active loads and number of system states from parameters array.
		loads_active = convert.(Integer, old_parameters[n_parameters -
								max_n_system_states+1:n_parameters])
		n_system_states = size(filter(x-> x != 0, loads_active))[1]

		# Initialize.
  		mcmc_samples[1:n_parameters, last_draw] = old_parameters[1:n_parameters]
  		mcmc_samples[n_parameters+1, last_draw] = n_system_states

		# Define an array to count the accepted proposals by the Metropolis-Hastings
		# (MH) step in the Gibbs sampler.
		accepted_samples_MH_counter = ones(n_parameters)

	else

		old_parameters = mcmc_samples[1:n_parameters, last_draw]

		#Compute active loads and number of system states from parameters array.
		loads_active = convert.(Integer, old_parameters[n_parameters -
											max_n_system_states+1:n_parameters])
		n_system_states = size(filter(x-> x != 0, loads_active))[1]

		# Initiate the accepted samples counter for Metropolis-Hastings step.
		accepted_samples_MH_counter = round.(last_draw .*
											mcmc_acceptance_rates[:, last_draw])

	end

	# Compute the log likelihood for the initializing sample (last draw)
	old_log_likelihood = get_log_likelihood(old_parameters, loads_active,
														n_system_states)

	# Below, the variable offset is needed because the first three elements
	# of the mcmc_* arrays defined above correspond to excitation rate and
	# relaxation rates of dyes.
	#
	# Also, we use the Bayes' rule to compute the posterior, that is,
	# posterior = likelihood * prior or log(posterior) = log(likelihood)
	# + log(prior)
	#
	offset = 3
	for i in 1:max_n_system_states
		for j in 1:max_n_system_states
		 	param = offset + (i-1)*max_n_system_states + j
			old_log_prior =	get_log_prior_rates(param, old_parameters[param],
														max_n_system_states)
			old_log_conditional_posterior = old_log_likelihood + old_log_prior
			mcmc_log_posteriors[param, last_draw] =
												old_log_conditional_posterior
			mcmc_acceptance_rates[param, last_draw] =
				accepted_samples_MH_counter[param]/convert(Float64, last_draw)
		end
	end

 	# log-Posterior associated with non-parametrics parameters
	offset = 3 + max_n_system_states^2 #new offset for loads
    	prior_success_probability = 1.0/(1.0 + ((max_n_system_states - 1)/
													expected_n_system_states))
	for i in 1:max_n_system_states
		param = offset + i
		if loads_active[i] != 0
			old_log_prior = log(prior_success_probability)
		elseif loads_active[i] == 0
			old_log_prior = log(1.0 - prior_success_probability)
		end

		old_log_conditional_posterior = old_log_prior + old_log_likelihood
		mcmc_log_posteriors[param, last_draw] = old_log_conditional_posterior
		mcmc_acceptance_rates[param, last_draw] =
				accepted_samples_MH_counter[param]/convert(Float64, last_draw)
	end

	# Add contribution from all the priors to get full joint posterior
	log_full_posterior = get_log_full_posterior(last_draw, old_log_likelihood,
													loads_active, mcmc_samples)
	mcmc_log_posteriors[n_parameters+1, last_draw] = log_full_posterior

	# This function checks for existing MCMC samples (files in the working
	# directory) and initiates the sampler from the last samples in those files.
	# When not doing that, it saves the new samples and other parameters in
	# newly generated files.
	#
	save_mcmc_data(last_draw, mcmc_samples,
									mcmc_log_posteriors, mcmc_acceptance_rates)
	print_and_plotting(last_draw, mcmc_samples, mcmc_log_posteriors,
														mcmc_acceptance_rates)

	# MCMC Part Starts Here.
	for draw in last_draw+1:total_draws

		# Initialize the parameters at current draw
 		mcmc_samples[:, draw] = mcmc_samples[:, draw-1]
 		mcmc_log_posteriors[:, draw] = mcmc_log_posteriors[:, draw-1]
 		mcmc_acceptance_rates[:, draw] = mcmc_acceptance_rates[:, draw-1]

		# Gibbs Sampling Step with Metropolis-Hastings starts here. Sequentially
		# generate samples from conditional distributions associated with each
		# parameter of interest.

		# Proposed new values (from multivariate normal distribution)
		proposed_parameters = propose_params(old_parameters, draw)


		# For FRET rates and system transition rates
		offset = 3
		for i in 1:max_n_system_states
			for j in 1:max_n_system_states

				param = offset + (i-1)*max_n_system_states + j

				if n_system_states > 0

					# Likelihood computation is only needed for active pair of
					# loads. This saves computational time. For an inactive
					# pair, the values are drawn simply from the priors and not
					# the posterior.

					if i in loads_active && j in loads_active

	 					# log-Posterior associated with old parameters
	 					old_log_prior = get_log_prior_rates(param,
									old_parameters[param], max_n_system_states)
						old_log_conditional_posterior = old_log_prior +
															old_log_likelihood

						# Define a new array with only the parameter
						# corresponding to current Gibbs step updated to
						# proposed value. This is needed because the proposal
						# may be rejected in the Metropolis-Hastings step.
						new_parameters = copy(old_parameters)
						new_parameters[param] = proposed_parameters[param]

	 					# log-Posterior associated with proposed parameters
						new_log_prior = get_log_prior_rates(param,
									new_parameters[param], max_n_system_states)
						new_log_likelihood = get_log_likelihood(new_parameters,
												loads_active, n_system_states)
						new_log_conditional_posterior = new_log_prior +
															new_log_likelihood

	 					# logarithm of Hastings Ratio
						log_hastings = new_log_conditional_posterior -
												old_log_conditional_posterior +
												log(new_parameters[param]) -
												log(old_parameters[param])

					else
 						# No likelihood computation needed this time
						old_log_prior = get_log_prior_rates(param,
									old_parameters[param], max_n_system_states)
						old_log_conditional_posterior = old_log_likelihood +
																old_log_prior

						new_parameters = copy(old_parameters)
						new_parameters[param] = proposed_parameters[param]

						new_log_likelihood = old_log_likelihood
						new_log_prior = get_log_prior_rates(param,
									new_parameters[param], max_n_system_states)
						new_log_conditional_posterior = new_log_likelihood +
																new_log_prior

						log_hastings = new_log_conditional_posterior -
												old_log_conditional_posterior +
													log(new_parameters[param]) -
													log(old_parameters[param])
					end

	 				# Move the chain forward
	 				if log_hastings >= log(rand())
						#if accepted, update all the arrays with new values
						accepted_samples_MH_counter[param] =
										accepted_samples_MH_counter[param] + 1.0
						old_parameters = copy(new_parameters)
						old_log_likelihood = new_log_likelihood
	 					old_log_conditional_posterior =
												new_log_conditional_posterior
	 				end

					# Store the samples and associated log_posterior, acceptance
					# ratio for the current draw
	 				mcmc_samples[param, draw] = old_parameters[param]
					mcmc_log_posteriors[param, draw] =
												old_log_conditional_posterior
					mcmc_acceptance_rates[param, draw] =
							accepted_samples_MH_counter[param]/
													convert(Float64, draw)


				elseif n_system_states == 0

					# Likelihood is 0 for a zero-state model. Automatic
					# rejection as Hastings ratio is 0.Furthermore, since we
					# already initialized the storage arrays earlier with values
					# from previous draw, we do not need to store them again.
					# Only acceptance ratio needs to be updated.
					#
					mcmc_acceptance_rates[param, draw] =
							accepted_samples_MH_counter[param]/
													convert(Float64, draw)

				end
			end
		end

		# Gibbs part for the active loads (model). We follow a similar approach
		# as for the parameters. No Metropolis-Hastings required since direct
		# sampling is available.
		#
		offset = 3 + max_n_system_states^2
		for i in 1:max_n_system_states

			param = offset + i

 			# Initialize proposed loads array
			proposed_loads_active = copy(loads_active)

			if proposed_loads_active[i] == 0

				log_likelihood_inactive = old_log_likelihood
				log_p_inactive = log_likelihood_inactive + log(1.0 -
													prior_success_probability)

				proposed_loads_active[i] = i
				n_system_states = size(filter(x-> x != 0,
													proposed_loads_active))[1]

				log_likelihood_active = get_log_likelihood(vec(old_parameters),
								      proposed_loads_active, n_system_states)

				log_p_active = log_likelihood_active +
												log(prior_success_probability)

			elseif proposed_loads_active[i] == i

				log_likelihood_active = old_log_likelihood
				log_p_active = log_likelihood_active +
												log(prior_success_probability)

 				# log-Posterior associated with old parameters
				proposed_loads_active[i] = 0
				n_system_states = size(filter(x-> x != 0,
													proposed_loads_active))[1]

				if n_system_states > 0

					log_likelihood_inactive = get_log_likelihood(
									vec(old_parameters),proposed_loads_active,
																n_system_states)

					log_p_inactive = log_likelihood_inactive +
											log(1.0 - prior_success_probability)

				elseif n_system_states == 0

					# Posterior is 0 for a zero-state model
					log_p_inactive = -Inf

				end

			end

			# The following procedure helps avoid overflow issues
			max_val = max(log_p_active, log_p_inactive)
			p_active = exp(log_p_active-max_val)
			p_inactive = exp(log_p_inactive-max_val)

			# Probability vector to activate or deactivate a load
			p_load = [p_active, p_inactive]
			# Normalize this probability vector
			p_load = p_load ./ sum(p_load)

			loads_active[i] =  rand(Categorical(p_load), 1)[1]
			if loads_active[i] == 1 #Active load

				loads_active[i] = i
				old_log_likelihood = log_likelihood_active
				old_log_conditional_posterior = log_p_active

			elseif loads_active[i] == 2 #Inactive load

				loads_active[i] = 0
				old_log_likelihood = log_likelihood_inactive
				old_log_conditional_posterior = log_p_inactive

			end


			n_system_states = size(filter(x-> x != 0, loads_active))[1]
			accepted_samples_MH_counter[param] =
										accepted_samples_MH_counter[param] + 1.0
			old_parameters[param] = convert(Float64, loads_active[i])

 			mcmc_samples[param, draw] = old_parameters[param]
			mcmc_log_posteriors[param, draw] = old_log_conditional_posterior
			# All samples accepted due to direct sampling
			mcmc_acceptance_rates[param, draw] = 1.0
		end

 		# Full joint posterior. Sequentially add log-priors to log-likelihood
		log_full_posterior = get_log_full_posterior(draw, old_log_likelihood,
													loads_active, mcmc_samples)

		# Update MCMC data arrays
		mcmc_log_posteriors[n_parameters+1, draw] = log_full_posterior
  		mcmc_samples[n_parameters+1, draw] = n_system_states

		# Save MCMC Data
		save_mcmc_data(draw, mcmc_samples, mcmc_log_posteriors,
														mcmc_acceptance_rates)
		# Print and plot MCMC data
		print_and_plotting(draw, mcmc_samples, mcmc_log_posteriors,
														mcmc_acceptance_rates)

 	end
	return nothing
end


# Start the sampler
total_draws = 10000
sampler(total_draws)
