# Compute generator matrix from the MCMC parameters. Generator matrix
# consists of all the photophysical rates (for donor and acceptor relaxation,
# and FRET) and transition rates for the biomolecule of interest.
function get_generator(parameters, loads_active)

	filtered_loads_active = filter(x -> x != 0, loads_active)
	n_system_states = size(filtered_loads_active)[1]

	# Excitation rate
	lambda_ex = parameters[1]


	# Rate matrix is created by adding one-by-one block matrices containing
	# photophysical rates and system transition rates. This is done to reduce
	# computational cost by keeping the rate matrix limited to active loads.
	# The blocks on the diagonal correspond to photophysical rates, and all
	# the off-diagonal blocks correspond to biomolecular transition rates.
	#
	# Furthermore, a "local" rate matrix needs to be defined because each "for"
	# loop iteration is independent from other iterations with regards to
	# scope of variables. So rate_matrix initialized in iteration 1 is undefined
	# in other iterations. You cannot really initialize variables inside a
	# for loop.
	#
 	offset = 3
	local rate_matrix
	for i in filtered_loads_active
		local rate_matrix_partial
		for j in filtered_loads_active
			# Here the rate matrix is built block by block. Blocks are added
			# horizontally and then vertically to the  matrix.
			#
			ij = offset + (i-1)*max_n_system_states+j

			if j == filtered_loads_active[1] # First active load
				# Photophysical rates
				if i == j
					lambda_FRET = parameters[ij]
 					rate_matrix_partial = [0.0 lambda_ex 0.0;
							      lambda_d 0.0 lambda_FRET;
							      lambda_a 0.0 0.0]
				end

				# System transition rates
				if i != j
					lambda_system = parameters[ij]
					rate_matrix_partial = lambda_system*Matrix(1.0I,3,3)
				end
			else
				# Photophysical rates
				if i == j
					lambda_FRET = parameters[ij]
 					rate_matrix_partial = hcat(rate_matrix_partial,
											[0.0 lambda_ex 0.0;
											lambda_d 0.0 lambda_FRET;
											lambda_a 0.0 0.0])
				end

				# System transition rates
				if i != j
					lambda_system = parameters[ij]
					rate_matrix_partial = hcat(rate_matrix_partial,
												lambda_system*Matrix(1.0I,3,3))
				end
			end
		end

		# Add the partial matrices to form the full rate matrix.
		if i == filtered_loads_active[1]
			rate_matrix = copy(rate_matrix_partial)
		end
		if i != filtered_loads_active[1]
			rate_matrix = vcat(rate_matrix, rate_matrix_partial)
		end
	end

	# Replace the zero diagonal elements with negative row-sums to get the
	# generator matrix that is used to compute the propagator
	#
	generator_FRET = copy(rate_matrix)
	length_generator_FRET = size(generator_FRET)[1]
	for i in 1:length_generator_FRET
 		generator_FRET[i,i] = -sum(generator_FRET[i,:])
	end

	# Background rate matrices for each channel.
	rate_matrix_bg_acceptor = lambda_bg[1]*[0.0 1.0;
						1.0 0.0]
	rate_matrix_bg_donor = lambda_bg[2]*[0.0 1.0;
						1.0 0.0]
	# Combine the two background matrices to get full background matrix
	id2 = Matrix(1.0*I, 2, 2) #Identity matrix
 	rate_matrix_bg = kron(rate_matrix_bg_acceptor, id2) +
											kron(id2, rate_matrix_bg_donor)

	# As before, replace zeroes in the diagonal with negative row-sums to get
	# the background generator matrix
	generator_bg = copy(rate_matrix_bg)
	length_generator_bg = size(generator_bg)[1]
	for i in 1:length_generator_bg
 		generator_bg[i,i] = -sum(generator_bg[i,:])
	end

	return generator_FRET, generator_bg
end

# Compute probability vector from the MCMC parameters. As shown in the
# research articles, this is done by assuming that the biomolecule and
# FRET pair composite are in equilibrium just before the experiment
# begins. This allows computing the initial probability vector needed
# to compute the likelihood by simply computing the kernel or the 
# eigenvector with zero eigenvalue of the generator matrix.
#
function get_rho(generator_FRET, generator_bg)

	length_generator_bg = size(generator_bg)[1]
	#Identity Matrix
	id_bg = Matrix(1.0*I, length_generator_bg, length_generator_bg)

	length_generator_FRET = size(generator_FRET)[1]
	#Identity Matrix
	id_FRET = Matrix(1.0*I, length_generator_FRET, length_generator_FRET)

	generator = kron(generator_bg, id_FRET) + kron(id_bg, generator_FRET)
	rho = permutedims(nullspace(Transpose(generator)))
	rho = rho/sum(rho)

	return rho
end

# Non-radiative Propagator.
function non_radiative_propagator(generator_FRET, generator_bg, n_system_states)

	length_generator_bg = size(generator_bg)[1]
	#Identity Matrix.
	id_bg = Matrix(1.0*I, length_generator_bg, length_generator_bg)

	length_generator_FRET = size(generator_FRET)[1]
	#Identity Matrix.
	id_FRET = Matrix(1.0*I, length_generator_FRET, length_generator_FRET)

	# Number of photo states (donor, acceptor): (ground, ground),
	# (ground, excited), (excited, ground).
	n_photo_states = 3

	# Get non-radiative detection matrix. Account for crosstalk/detection
	# efficiency. Set the probabilities for radiative transitions to zero
	# or (1.0-detection efficiency).
	detection_matrix_non_FRET = ones(length_generator_FRET,
														length_generator_FRET)
	for  i in 0:n_system_states-1
		detection_matrix_non_FRET[i*n_photo_states+2, i*n_photo_states+1] = 0.0
		detection_matrix_non_FRET[i*n_photo_states+3, i*n_photo_states+1] =
																1.0 - phi[1, 1]
	end

	# .* is the element-by-element product to get reduced propagators.
	generator_non_FRET = generator_FRET .* detection_matrix_non_FRET
	reduced_propagator_FRET = exp(generator_non_FRET)

 	detection_matrix_non_bg = id_bg
	generator_non_bg = generator_bg .* detection_matrix_non_bg
	reduced_propagator_bg = exp(generator_non_bg)

	reduced_propagator = kron(reduced_propagator_bg, reduced_propagator_FRET)
	return reduced_propagator
end

# Radiative propagator
function radiative_propagator(generator_FRET, generator_bg, detection_channel,
																n_system_states)

	length_generator_bg = size(generator_bg)[1]
	#Identity Matrix
	id_bg = Matrix(1.0*I, length_generator_bg, length_generator_bg)

	length_generator_FRET = size(generator_FRET)[1]
	#Identity Matrix
	id_FRET = Matrix(1.0*I, length_generator_FRET, length_generator_FRET)

	id2 = Matrix(1.0*I, 2, 2)

	# Number of photo states (donor, acceptor): (ground, ground),
	# (ground, excited), (excited, ground)
	n_photo_states = 3

	# Get reduced generator matrices
 	if detection_channel == 1

		# Get radiative detection matrix if detection in acceptor channel.
		# Apply crosstalk/detection efficiency. Set probabilities
		# for nonradiative transitions to zero.
		#
		detection_matrix_acceptor_FRET = zeros(length_generator_FRET,
														length_generator_FRET)
		for  i in 0:n_system_states-1
			# Probability that donor emitted the photon detected in
			# acceptor channel.
			#
			detection_matrix_acceptor_FRET[i*n_photo_states+2,
												i*n_photo_states+1] = phi[2, 1]
			# Probability that acceptor emitted the photon detected in
			# acceptor channel.
			#
			detection_matrix_acceptor_FRET[i*n_photo_states+3,
												i*n_photo_states+1] = phi[1, 1]
		end

 		generator_rad_FRET = generator_FRET .* detection_matrix_acceptor_FRET
		detection_matrix_rad_bg = kron([0.0 1.0;
				 		1.0 0.0], id2)

 		generator_rad_bg = generator_bg .* detection_matrix_rad_bg

	elseif detection_channel == 2

		# Get radiative detection matrix if detection in donor channel.
		# Apply crosstalk/detection efficiency.Set probabilities
		# for nonradiative transitions to zero.
		#
		detection_matrix_donor_FRET = zeros(length_generator_FRET,
														length_generator_FRET)
		for  i in 0:n_system_states-1
			# Probability that donor emitted the photon detected in donor
			# channel.
			#
			detection_matrix_donor_FRET[i*n_photo_states+2,
												i*n_photo_states+1] = phi[2, 2]
			# Probability that acceptor emitted the photon detected in donor
			# channel.
			#
			detection_matrix_donor_FRET[i*n_photo_states+3,
												i*n_photo_states+1] = phi[1, 2]
		end


 		generator_rad_FRET = generator_FRET .* detection_matrix_donor_FRET
		detection_matrix_rad_bg = kron(id2, [0.0 1.0;
				      			1.0 0.0])
 		generator_rad_bg = generator_bg .* detection_matrix_rad_bg

 	end

 	reduced_propagator = kron(generator_rad_bg, id_FRET) + kron(id_bg,
															generator_rad_FRET)

	return reduced_propagator
end
