########## GRN Model Setup ######### 

# Experiment description: In this experiment, a mutual inhibition network is chosen. We run repeated evolution, selecting for a single central stripe, whilst restricting mutation to only take place within the topology. 

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

include(srcdir("TissueModel_ND.jl"))

########## data load ######### 

networks = load(datadir("networks/FindNetworks_HalfStripeLeft_RawData.jld2"));

topology_choice = "mutual_inh"

start_network = networks[topology_choice * "_networks"][argmin(networks[topology_choice * "_t2s"])]

########## Evolutionary Setup ######### 

β = 1.

noise_cv = 1.

mut_prob = 0.1

deletion_prob = 0.

grn_parameters = DefaultGRNParameters();

max_w = 10.

output_gene = 3

n_stripe = 1

fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_stripe),output_gene);

tolerance = 0.95

viable_mutations = Int64.(start_network .!= 0)

mutation_weights = findall(viable_mutations .> 0)

n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),n_sample_func,deletion_prob,max_w,mutation_weights)

mutate_function = i -> noise_no_additions(i,mutation_op);

######### Simulation setup ######### 

n_traj = 10000
max_gen = 40000