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

# start_networks_dict =  load(datadirx("networks/FindNetworks_HalfStripeRight_Full_RawData.jld2"));

# topology_choice = "bistable"

# choice = 1 # full

start_network = [0.0 0.0 0.0 0.07073458757269772; 0.6889253020855619 0.0 -1.042136950503097 0.0; 4.172293675117961 -2.684028190193488 0.0 0.0] # mutual_inhibition gradient

########## Topologies ###########

# These are taken from: Cotterell, J., & Sharpe, J. (2010). An atlas of gene regulatory networks reveals multiple three‐gene mechanisms for interpreting morphogen gradients. Molecular systems biology, 6(1), 425.

w_feed_forward = [0 0 0 1 ; 1 0 0 0 ; 1 -1 1 0];
w_mutual_inh = [0 0 0 1 ; 1 0 -1 0 ; 1 -1 0 0];
w_frozen_osc = [1 0 0 0; -1 0 1 0; -1 -1 1 1];
w_overlap_dom = [0 0 -1 1 ; 1 0 0 0 ; -1 1 0 0];
w_bistable = [0 0 0 1; 0 1 -1 0; -1 1 0 0];
w_classical = [0 0 0 1 ; -1 1 0 0 ; -1 -1 1 0];

########## Evolutionary Setup ######### 
 
β = 1.

noise_cv = 1.

mut_prob = 0.1

deletion_prob = 0.1

grn_parameters = DefaultGRNParameters();

max_w = 10.

output_gene = 3

n_stripe = 1

min_width = 2.

lower_bound = 5.

upper_bound = 10.

fitness_function_left = s -> fitness_evaluation(s,x->(nstripe_fitness(x,0,min_width,lower_bound,upper_bound),malt_fitness_left(x)),output_gene);

fitness_function_stripe = s -> fitness_evaluation(s,x->(nstripe_fitness(x,1,min_width,lower_bound,upper_bound),malt_fitness(x,1)),output_gene);

fitness_function_list = [fitness_function_left,fitness_function_stripe]

tolerance_list = [0.9,0.8]

viable_mutations = ones(Int,Ng,Ng+1)

viable_mutations[2,4] = 0
viable_mutations[3,4] = 0

mutation_weights = findall(viable_mutations .> 0)

n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),n_sample_func,deletion_prob,max_w,mutation_weights)

mutate_function = i -> noise(i,mutation_op);

########## Dyn Setup ######### 

save_id = [CartesianIndex(1,25),CartesianIndex(2,25),CartesianIndex(3,25),CartesianIndex(1,50),CartesianIndex(2,50),CartesianIndex(3,50),CartesianIndex(1,100),CartesianIndex(2,100),CartesianIndex(3,100)]
n_segments = 4
n_steps = 10

d_metric = Euclidean()
relative_dyn = true

fundamental_topologies =  ["feed_forward","mutual_inh","frozen_osc","bistable","classical"]

fundamental_networks_dict_1 = load(datadirx("networks/FindNetworks_HalfStripeLeft_Full_RawData.jld2"));

fundamental_networks_1 = reduce(vcat,[fundamental_networks_dict_1[top_choice * "_networks"] for top_choice in fundamental_topologies])
fundamental_networks_t2s_1 = reduce(vcat,[fundamental_networks_dict_1[top_choice * "_t2s"] for top_choice in fundamental_topologies])
fundamental_labels_1 = reduce(vcat,[[top_choice for _ in 1:length(fundamental_networks_dict_1[top_choice * "_networks"])] for top_choice in fundamental_topologies])

fundamental_networks_dict_2 = load(datadirx("networks/FindNetworks_CentreStripe_Full_RawData.jld2"));

fundamental_networks_2 = reduce(vcat,[fundamental_networks_dict_2[top_choice * "_networks"] for top_choice in fundamental_topologies])
fundamental_networks_t2s_2 = reduce(vcat,[fundamental_networks_dict_2[top_choice * "_t2s"] for top_choice in fundamental_topologies])
fundamental_labels_2 = reduce(vcat,[[top_choice for _ in 1:length(fundamental_networks_dict_2[top_choice * "_networks"])] for top_choice in fundamental_topologies])

fundamental_networks_list = [fundamental_networks_1,fundamental_networks_2]
fundamental_networks_t2s_list = [fundamental_networks_t2s_1,fundamental_networks_t2s_2]
fundamental_labels_list = [fundamental_labels_1,fundamental_labels_2]

n_fundamental_network_list = [length(fundamental_networks_1),length(fundamental_networks_2)]

########## LMC Setup ######### 

N_interp_points = 10

development = DefaultGRNSolver()

######### Simulation setup ######### 

n_trials = 2500
max_gen = 75000