# Global input parameters for continuous illumination edition of BNP-FRET.
# All the time scales are set in nanosecond units.

# Directory containing FRET traces
#working_directory =
#        "/home/singularity/dropbox/bayesianCodes/julia/multiFRET/BNP-FRET/"
working_directory =
        "/home/mbgnjasb/Dropbox (ASU)/bayesianCodes/julia/multiFRET/BNP-FRET/"


# Input data files are named using the convention "biomolecule_label.h5".
# When running multiple parallelized MCMC chains, it is convenient to use
# the process ID itself as the label. Therefore labels are typically set in
# the function "get_data()"
filename_prefix = "EG_"
filename_suffix = ".h5"


# Total photons from the FRET trace to be used for analysis
total_photons = 1000

# Experimental Parameters
n_channels = 2

# Acceptor channels are labeled with odd numbers (1, 3, ..)
# and Donor channels are labeled with even numbers (2, 4, ..)

# Background rates (in ns^(-1) units) for each channel
lambda_bg = zeros(n_channels)
lambda_bg[1] = 0.4672e-6 #Acceptor channel
lambda_bg[2] = 0.2828e-6 #Donor channel


# ============== Crosstalk factors ===================
# phi[i, j] represents probability for photons intended
# for channel i to enter channel j. It incorporates
# effects of crosstalk, detection efficiency, and
# quantum yield.
# ====================================================
#
phi = zeros(n_channels, n_channels)
phi[1, 1] = 0.835897
phi[1, 2] = 0.0
phi[2, 2] = 0.820513
phi[2, 1] = 1.0 - phi[2, 2]
# Route correction matrix.
RCM = inv(permutedims(phi))


n_dyes = 2

# Acceptor relaxaton rate
lifetime_a = 2.8 #nanoseconds
lambda_a = 1.0/lifetime_a

# Donor relaxaton rate
lifetime_d = 1.22 #nanoseconds
lambda_d = 1.0/lifetime_d

# Prior information on system transition rates. Typical time scale for
# system transtions. Use default value if not known.
typical_time_scale = 1.0e9 #nanoseconds
typical_system_transition_rate = 1.0/typical_time_scale



# Nonparametrics parameters.
max_n_system_states = 6
expected_n_system_states = 2


# Markov Chain Monte Carlo (MCMC) Parmeters. These parameters set
# the size of the array containg all random variables (parameters of interest),
# which are excitation rate, relaxation rates of the dyes, FRET rates, system
# transtion rates, and loads.
n_rates = 1 + n_dyes + max_n_system_states^2
n_parameters = n_rates + max_n_system_states

# Plotting options.
turn_on_plotting = true # Choose true or false
plotting_frequency = 2
