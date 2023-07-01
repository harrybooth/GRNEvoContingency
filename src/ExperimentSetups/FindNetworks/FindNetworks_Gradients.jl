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

######### GRN Topology setup ######### 

# These are taken from: Cotterell, J., & Sharpe, J. (2010). An atlas of gene regulatory networks reveals multiple three‐gene mechanisms for interpreting morphogen gradients. Molecular systems biology, 6(1), 425.

w_feed_forward = [0 0 0 1 ; 1 0 0 0 ; 1 -1 1 0];
w_mutual_inh = [0 0 0 1 ; 1 0 -1 0 ; 1 -1 0 0];
w_frozen_osc = [1 0 0 0; -1 0 1 0; -1 -1 1 1];
w_overlap_dom = [0 0 -1 1 ; 1 0 0 0 ; -1 1 0 0];
w_bistable = [0 0 0 1; 0 1 -1 0; -1 1 0 0];
w_classical = [0 0 0 1 ; -1 1 0 0 ; -1 -1 1 0];

network_topology_dict = Dict("feed_forward"=>w_feed_forward,"mutual_inh"=>w_mutual_inh,"frozen_osc"=>w_frozen_osc,"overlap_dom"=>w_overlap_dom,"bistable"=>w_bistable,"classical"=>w_classical)

# networks_to_search = ["feed_forward","mutual_inh","frozen_osc","overlap_dom","bistable","classical"]

networks_to_search = ["feed_forward","mutual_inh"]

########## Evolutionary Setup ######### 

β = 1.

noise_cv = 1.

mut_prob = 0.1

deletion_prob = 0.

grn_parameters = DefaultGRNParameters();

max_w = 10.

output_gene = 3

n_stripe = 1

fitness_function = s -> fitness_evaluation(s,x->gradient_fitness_r(x),output_gene);

tolerance = -1.

########## To define in script ######### 

# viable_mutations = Int64.(start_network .!= 0)

# mutation_weights = findall(viable_mutations .> 0)

# n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

# mutation_op = MutationOperator(Normal,(μ = 0.0,σ = 1.),n_sample_func,deletion_prob,max_w,mutation_weights)

# mutate_function = i -> noise_no_additions(i,mutation_op);

########## Simulation Setup ######### 

n_networks_required = 50

max_gen = 10000

full_networks_req = false

param_N = 1
