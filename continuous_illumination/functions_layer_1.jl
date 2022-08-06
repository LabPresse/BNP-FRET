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


# Function to obtain the FRET data in the format containing an array full of
# photon arrival times and associated detection channels
function get_FRET_data()
		file_name = string(working_directory, filename_prefix, 316,
															filename_suffix)
		fid = h5open(file_name,"r")
		photon_arrival_times = read(fid, "photon_arrival_time_(ns)")
		detection_channels = read(fid, "detection_channel")
		close(fid)
		return photon_arrival_times, detection_channels
end

# Function to initialize parameters of interest to start the MCMC chain
function initialize_params()

	# Initialize initial_parameters array for the sampler. This array contains
	# the values for all the parameters of interest at the current draw in the
	# MCMC chain. The excitation rate, donor relaxation rate, and acceptor
	# relaxation rate constitute the first three elements of this array. Next,
	# to store the FRET rates and system transitions rates, we first replace
	# the diagonal elements of the system transition rate matrix with FRET rates
	# associated with the specific rows' system states. Then we store these
	# parameters in the old_parameters array row-by-row.
	#
	# In summary, with max_n_system_states be the maximum number of system
	# states used for nonparametrics we have:
	#
	# initial_parameters[1] = excitation rate,
	# initial_parameters[2] = donor relaxation rate,
	# initial_parameteres[3] = acceptor relaxation rates, and
	#
	# with ij = (i-1)*max_n_system + j, where i represents row and j represents
	# column of the rate matrix,
	#
	# initial_parameters[3 + ij] = FRET rates, if i = j,
	# initial_parameters[3 + ij] = system transition rates for system state "i",
	# if i != j.
	#
	# All the parameters are stored this way to keep the structure of the
	# sampler as general as possible.
	#
	initial_parameters = zeros(n_parameters)

	# All rates here are in ns^-1 units.

	# To obtain excitation rate, we take reciprocal of average interphoton
	# arrival time while correcting for crosstalk and background.
	#
	data_donor = []
	data_acceptor = []
	for i in 1:total_photons

     	if detection_channels[i] == 1 # Acceptor
			data_acceptor = vcat(data_acceptor, photon_arrival_times[i])
       	end

       	if detection_channels[i] == 2 # Donor
			data_donor = vcat(data_donor, photon_arrival_times[i])
       	end

   	end

	# Duration of FRET trace.
	duration_FRET_trace = photon_arrival_times[total_photons] -
											photon_arrival_times[1]

	# First correct for background.
	corrected_acceptor_photons = size(data_acceptor)[1] - lambda_bg[1] *
											duration_FRET_trace
	corrected_donor_photons = size(data_donor)[1] - lambda_bg[2] *
											duration_FRET_trace

	# Then correct for detection efficiency/crosstalk.
	corrected_total_photons = (RCM[1, 1]*corrected_acceptor_photons +
								(RCM[1, 2] + RCM[2,2])*corrected_donor_photons)

	# Effective excitation rate.
	lambda_ex =  corrected_total_photons/duration_FRET_trace
	initial_parameters[1] = lambda_ex

	# Initalize donor and acceptor relaxation rates.
	initial_parameters[2] = lambda_d
	initial_parameters[3] = lambda_a


	# Initialize all FRET rates and system transition rates.
	#
	# Below, the variable offset is needed because the first three elements
	# of the mcmc_* arrays defined above correspond to excitation rate and
	# relaxation rates of dyes.
	#
	offset = 3 #new offset for loads
	for i in 1:max_n_system_states
		for j in 1:max_n_system_states
			ij = offset + (i - 1) * max_n_system_states + j
 			if i == j # FRET rates
				initial_parameters[ij] = 1.0
			elseif i != j # System transition rates
				initial_parameters[ij] = typical_system_transition_rate
			end

		end
	end

	# Initialize loads for nonparametrics
	offset = 3 + max_n_system_states^2
	if modeling_choice == "nonparametric"
		prior_success_probability = 1.0/(1.0 + ((max_n_system_states - 1)/
													expected_n_system_states))
		p_load = [prior_success_probability, 1.0 - prior_success_probability]
		n_system_states = 0
		for i in 1:max_n_system_states
			initial_parameters[offset+i] = rand(Categorical(p_load), 1)[1]
			if  initial_parameters[offset+i]== 1 #Active
				initial_parameters[offset+i] = i
				n_system_states = n_system_states + 1
			elseif initial_parameters[offset+i] == 2 #Inactive
				initial_parameters[offset+i] = 0
				if i == max_n_system_states && n_system_states == 0
					initial_parameters[offset+i] = i
				end
			end
		end
	else
		for i in 1:expected_n_system_states
			initial_parameters[offset+i] = i
		end
		n_system_states = expected_n_system_states
	end

	return initial_parameters
end

# Function to compute next values in the chain from normal proposal
# distributions.
function propose_params(old_parameters, draw)

	#Initialize
	proposed_parameters = copy(old_parameters)

	# Set the covariance for excitation_rate, donor, and acceptor rates
	# covariance_ex = 1.0e-5
	# proposed_parameters[1] = rand(Normal(log(proposed_parameters[1]),
	# 													covariance_ex),1)[1]
	# proposed_parameters[1] = exp(proposed_parameters[1])

	# Set covariances for FRET rates and system transtion rates. Alternate
	# between two sets of values to explore the parameter space more efficiently
	# for multiscale analysis
	if draw % 2 == 0

		covariance_FRET = 0.5
		covariance_system = 5.0

	else

		covariance_FRET = 1.0e-2
		covariance_system = 1.0e-1
	end


	# Sample new proposals from normal distributions in the log space
	offset = 3
	for i in 1:max_n_system_states
		for j in 1:max_n_system_states
			ij = offset + (i - 1) * max_n_system_states + j
			if i == j # FRET rates
				proposed_parameters[ij] = rand(Normal(log(
								proposed_parameters[ij]), covariance_FRET),1)[1]
				proposed_parameters[ij] = exp(proposed_parameters[ij])
			elseif i != j # System transition rates
				proposed_parameters[ij] = rand(Normal(log(
							proposed_parameters[ij]), covariance_system),1)[1]
				proposed_parameters[ij] = exp(proposed_parameters[ij])
			end

		end
	end

	return proposed_parameters
end

# Logarithm of likelihood
function get_log_likelihood(parameters, loads_active, n_system_states)

	# Get full generator matrix and initial probability vector.
 	generator_FRET, generator_bg = get_generator(parameters, loads_active)
	rho = get_rho(generator_FRET, generator_bg)
	p = 1.0
	logL = 0.0

	for i in 1:total_photons

		if i == 1
			inter_arrival_time = photon_arrival_times[i]
		else
			inter_arrival_time = photon_arrival_times[i] -
													photon_arrival_times[i-1]
		end

		# Non-radiative part
		tau = inter_arrival_time
		Q = non_radiative_propagator(tau*generator_FRET, tau*generator_bg,
															n_system_states)
		rho = rho * Q
		p = sum(rho)
		logL = logL + log(p)
		rho = rho/p

		# Radiative part
		Q = radiative_propagator(generator_FRET, generator_bg,
										detection_channels[i], n_system_states)
		rho = rho * Q
		p = sum(rho)
 		logL = logL + log(p)
		rho = rho/p

	end

	return logL
end

# log-prior for rates.
function get_log_prior_rates(param, value, max_n_system_states)

	# For excitation rate
	if param == 1
		log_prior = logpdf(Gamma(1, 1.0e-6), value)
	# For donor and acceptor dyes' relaxation rates
	elseif param == 2 || param == 3
		log_prior = logpdf(Gamma(1, 1.0), value)
	end

	offset = 3
	j = (param - offset) % max_n_system_states
	if j == 0
		j = max_n_system_states
	end
	i = (param - offset - j)/max_n_system_states + 1


	if i == j
		# For FRET rates
		log_prior = logpdf(Gamma(1, 1.0), value)
	elseif i != j
		# For system transition rates
		log_prior = logpdf(Gamma(1, typical_system_transition_rate), value)
	end

	return log_prior
end

# Compute full joint posterior which is equal to likelihood plus all
# the prior densities corresponding to current values of transition
# rates. This includes rates for active and inactive loads.
function get_log_full_posterior(draw, old_log_likelihood, loads_active,
																mcmc_samples)
 	# Full joint posterior. Sequentially add log-priors to log-likelihood
	log_full_posterior = old_log_likelihood
	offset = 3
	for i in 1:max_n_system_states
		for j in 1:max_n_system_states
			param = offset + (i-1)*max_n_system_states + j
			log_prior = get_log_prior_rates(param, mcmc_samples[param, draw],
															max_n_system_states)
			log_full_posterior = log_full_posterior + log_prior
		end
	end

	offset = 3 + max_n_system_states^2 #new offset for loads
    	prior_success_probability = 1.0/(1.0 + ((max_n_system_states - 1)/
													expected_n_system_states))
	for i in 1:max_n_system_states
		param = offset + i
		if loads_active[i] != 0
			log_prior = log(prior_success_probability)
		elseif loads_active[i] == 0
			log_prior = log(1.0 - prior_success_probability)
		end
		log_full_posterior = log_full_posterior + log_prior
	end
	return log_full_posterior
end




# Print Stuff and Plot MCMC output
function print_and_plotting(current_draw, mcmc_samples, mcmc_log_posteriors,
														mcmc_acceptance_rates)

 	if current_draw == 1

		println( "    ")
		println( "Modeling choice = ", modeling_choice)
		println( "    ")
		println( "Maximum number of system states = ", max_n_system_states)
		println( "    ")
		println( "Expected number of system states = ",
													expected_n_system_states)
		println( " ")
		println( "Success probability = ", 1.0/(1.0 + ((max_n_system_states -
																	1)/2.0)))
		println( " ")

	end

	n_system_states = convert(Int, mcmc_samples[n_parameters+1, current_draw])
	loads_active = convert.(Int, mcmc_samples[n_rates+1:n_rates +
											max_n_system_states, current_draw])
	filtered_loads_active = filter(x -> x != 0, loads_active)

	rate_matrix = zeros(n_system_states, n_system_states)
	acceptance_rates = ones(n_system_states, n_system_states)
	offset = 3
	for i in 1:n_system_states
		for j in 1:n_system_states
			param = offset + (filtered_loads_active[i]-1)*max_n_system_states +
														filtered_loads_active[j]
			#for conversion to s^(-1) units
			rate_matrix[i, j] = mcmc_samples[param, current_draw] * 1.0e9
			if i == j
				rate_matrix[i, j] = rate_matrix[i, j]/(rate_matrix[i, j] +
															lambda_d * 1.0e9)
			end
			acceptance_rates[i, j] = mcmc_acceptance_rates[param, current_draw]
		end
	end

	log_full_posterior = mcmc_log_posteriors[n_parameters+1, current_draw]


	# Print current values to standard out (terminal)
	println("=================================================================")
	println(" ")
	println("Iteration: ", current_draw)
	println(" ")
	println("n_system_states = ", n_system_states)
	println(" ")
	println("loads_active = ", filtered_loads_active)
	println(" ")
	println("excitation rate (s^(-1)) = ", mcmc_samples[1, current_draw] * 1.0e9)
	println(" ")
	println("rate_matrix (s^(-1)) (diagonal elements show FRET efficiencies) =")
	println(" ")
	display(rate_matrix)
	println(" ")
	println(" ")
	println("log_full_posterior = ", log_full_posterior)
	println(" ")
	println("minimum, median, maximum acceptance rate = ",
										minimum(acceptance_rates)*100, "%  ",
										median(acceptance_rates)*100, "%  ",
										maximum(acceptance_rates)*100, "%  ")
	println(" ")

	# Plotting
	if (turn_on_plotting == true) && (current_draw % plotting_frequency == 0)

		plot_log_full_posterior = plot(mcmc_log_posteriors[n_parameters+1,
							1:current_draw],
			     			xlabel="Iterations", ylabel="log(posterior)",
			     			legend = false, linewidth= 1.5,
						xlims = (1, current_draw),
						xticks = 1:Int(ceil((current_draw-1)/5)):current_draw,
						xtickfontsize=18, ytickfontsize=18,
						xguidefontsize = 20, yguidefontsize = 20,
						fontfamily = "Computer Modern",
						yformatter = :scientific,
						right_margin=5mm, bottom_margin = 5mm,
						left_margin = 5mm,
						dpi = 300, size = (1000, 600))

		plot_n_system_states = histogram(mcmc_samples[n_parameters+1,
						1:current_draw],
						xlabel="Number of System States", ylabel="Probability",
						legend = false, linewidth= 1.5,
						normalize = :probability,
						bins = 0.5:max_n_system_states+0.5, bar_width = 0.2,
						xticks = 0:Int(ceil((max_n_system_states-1)/5)):
															max_n_system_states,
						xtickfontsize=18, ytickfontsize=18,
						xguidefontsize = 20, yguidefontsize = 20,
						fontfamily = "Computer Modern",
						right_margin=5mm, bottom_margin = 5mm,
						left_margin = 5mm,
						dpi = 300, size = (1000, 600))


		# Gather all samples for escape rates and FRET rates for active loads,
		# and collect the values of n_system_states.
		#
		data_FRET = []
		data_Escape = []

		offset = 3
		for draw in 1:current_draw
			loads_active = convert.(Int, mcmc_samples[n_rates+1:n_rates +
													max_n_system_states, draw])
			for i in 1:max_n_system_states
				escape_rate = 0.0
				for j in 1: max_n_system_states
					param = offset + (loads_active[i]-1)*max_n_system_states +
																loads_active[j]
					if loads_active[i] * loads_active[j] != 0
						if i == j
							data_FRET = vcat(data_FRET,
													mcmc_samples[param, draw])
						elseif i != j
							escape_rate = escape_rate +
													mcmc_samples[param, draw]
						end
					end
				end

 				if loads_active[i] != 0
					data_Escape = vcat(data_Escape, escape_rate)
				end
			end
		end

		# Gather all FRET efficiencies.
		data_FRET = convert.(Float64, data_FRET)
		data_FRET_Efficiency = data_FRET ./ (data_FRET .+ lambda_d)
		# Gather all escape rates.
		data_Escape = convert.(Float64, data_Escape)
		data_Escape = data_Escape .* 1.0e9

		dens = kde((data_FRET_Efficiency, data_Escape))
		plot_bivariate_posterior = heatmap(dens,
					xlabel="\$\\epsilon_{FRET}\$",
					ylabel="\$\\lambda_{esc} \\, (s^{-1})\$",
					linewidth= 1.5,
					xtickfontsize=18, ytickfontsize=18,
					xguidefontsize = 20, yguidefontsize = 20,
					fontfamily = "Computer Modern",
					right_margin=5mm, bottom_margin = 3mm, left_margin = 0mm,
					dpi = 300, size = (600, 600))


		l = @layout [[a; b] c]
		display(plot(plot_log_full_posterior, plot_n_system_states,
					plot_bivariate_posterior,
 			     	plot_title = "BNP-FRET (Continuous)",
					plot_titlefontsize = 34,
					right_margin=10mm, bottom_margin = 10mm, left_margin = 10mm,
					size=(2000,1000), layout = l))

	end
	return nothing
end

function check_for_existing_mcmc_data(mcmc_samples,
									mcmc_log_posteriors, mcmc_acceptance_rates)

	file_name = string(working_directory, "mcmc_output_", filename_prefix,
											total_photons, filename_suffix)
	if isfile(file_name) == true

		# Check for existing MCMC data when starting the sampler. If it exists,
		# read the data  and compute the necessary parameters to initiate
		# the sampler. The files are assumed to be in HDF5 format.

		fid = h5open(file_name,"r")
		old_mcmc_samples = read(fid, "mcmc_samples")
		old_mcmc_log_posteriors = read(fid, "mcmc_log_posteriors")
		old_mcmc_acceptance_rates = read(fid, "mcmc_acceptance_rates")
		close(fid)

		last_draw = size(old_mcmc_samples)[2]
		mcmc_samples[:, 1:last_draw] = old_mcmc_samples[:, 1:last_draw]
		mcmc_log_posteriors[:, 1:last_draw] =
										old_mcmc_log_posteriors[:, 1:last_draw]
		mcmc_acceptance_rates[:, 1:last_draw] =
									old_mcmc_acceptance_rates[:, 1:last_draw]

	else

		last_draw = 1

	end

	return  last_draw
end

function save_mcmc_data(current_draw, mcmc_samples,
									mcmc_log_posteriors, mcmc_acceptance_rates)

	# Save the data in HDF5 format.
	file_name = string(working_directory, "mcmc_output_", filename_prefix,
											total_photons, filename_suffix)

	fid = h5open(file_name,"w")
	write_dataset(fid, "mcmc_samples",
									mcmc_samples[:, 1:current_draw])
	write_dataset(fid, "mcmc_log_posteriors",
								mcmc_log_posteriors[:, 1:current_draw])
	write_dataset(fid, "mcmc_acceptance_rates",
								mcmc_acceptance_rates[:, 1:current_draw])
	close(fid)

	return nothing
end
