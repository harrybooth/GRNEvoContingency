########## GRN Model Setup ######### 

# Experiment description: In this experiment, a mutual inhibition network is chosen. We run repeated evolution, selecting for a single central stripe.

using Distances

########## GRN Model Setup ######### 

const Nc = 100
const Ng = 3
const L = 1.

const m0 = 1.
const λm = 0.2*L

const θ = 5.
const deg_rate_g = 0.05
const init_conc_g = 0.1 

const tissue = range(0,L,length = Nc)

morph(x) = m0*exp(-x/λm)

σ(I) = 1/(1+exp(θ-θ*I))  # σ(0.) > 0 ?

include(srcdirx("TissueModel_ND.jl"))

########## data load ######### 

start_network = [0.0 0.0 0.0 1.2490335893436255; 0.0 0.0 0.0 0.0; -0.21577059555519695 0.0 0.0 0.0]

start_top = [0 0 0 1; 0 0 0 0; -1 0 0 0]

########## Topologies ###########

# These are taken from: Cotterell, J., & Sharpe, J. (2010). An atlas of gene regulatory networks reveals multiple three‐gene mechanisms for interpreting morphogen gradients. Molecular systems biology, 6(1), 425.

w_feed_forward = [0 0 0 1 ; 1 0 0 0 ; 1 -1 1 0];
w_mutual_inh = [0 0 0 1 ; 1 0 -1 0 ; 1 -1 0 0];
w_frozen_osc = [1 0 0 0; -1 0 1 0; -1 -1 1 1];
w_overlap_dom = [0 0 -1 1 ; 1 0 0 0 ; -1 1 0 0];
w_bistable = [0 0 0 1; 0 1 -1 0; -1 1 0 0];
w_classical = [0 0 0 1 ; -1 1 0 0 ; -1 -1 1 0];

network_topology_dict = Dict("feed_forward"=>w_feed_forward,"mutual_inh"=>w_mutual_inh,"frozen_osc"=>w_frozen_osc,"overlap_dom"=>w_overlap_dom,"bistable"=>w_bistable,"classical"=>w_classical)

########## Evolutionary Setup ######### 

β = (1.,10000)

grn_parameters = DefaultGRNParameters();

max_w = 10.

output_gene = 3

n_stripe = 1

min_width = 5

lower_bound = 5.

upper_bound = 10.

max_conc = 20.

fitness_function = s -> fitness_evaluation(s,x->(stripe_indicator(x,min_width,lower_bound,upper_bound),malt_fitness_absolute(x,n_stripe,max_conc)),output_gene);

# fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_stripe),output_gene);

tolerance = 0.9

mut_prob = 0.025

min_affinity = 1e-3

sign_flip_probability = 0.5

pm_prob = 0.5

viable_mutations = ones(Int,Ng,Ng+1)

viable_mutations[2,4] = 0
viable_mutations[3,4] = 0

mutation_weights = findall(viable_mutations .> 0)

mult_σ = 0.2

mult_noise_distribution = LogNormal(-0.1,mult_σ)

add_σ = 1.

additive_noise_distribution = Normal(0.,add_σ)

n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

mutation_op = MutationOperatorDual(mult_noise_distribution,additive_noise_distribution,n_sample_func,pm_prob,min_affinity,max_w,mutation_weights,sign_flip_probability);

mutate_function = i -> noise_mtype_mult_add(i,mutation_op)

########## Dyn Setup ######### 

save_id = [CartesianIndex(1,25),CartesianIndex(2,25),CartesianIndex(3,25),CartesianIndex(1,50),CartesianIndex(2,50),CartesianIndex(3,50),CartesianIndex(1,100),CartesianIndex(2,100),CartesianIndex(3,100)]
n_segments = 4
n_steps = 10

d_metric = Euclidean()
relative_dyn = true

fundamental_networks_dict = load(datadirx("networks/FindNetworks_CentreStripe_Full_RawData.jld2"));

fundamental_topologies =  ["feed_forward","mutual_inh","frozen_osc","bistable","classical"]

fundamental_networks = reduce(vcat,[fundamental_networks_dict[top_choice * "_networks"] for top_choice in fundamental_topologies])
fundamental_networks_t2s = reduce(vcat,[fundamental_networks_dict[top_choice * "_t2s"] for top_choice in fundamental_topologies])
fundamental_labels = reduce(vcat,[[top_choice for _ in 1:length(fundamental_networks_dict[top_choice * "_networks"])] for top_choice in fundamental_topologies])

n_fundamental_networks = length(fundamental_networks)

######### Simulation setup ######### 

n_trials = 10000
max_gen = 500000