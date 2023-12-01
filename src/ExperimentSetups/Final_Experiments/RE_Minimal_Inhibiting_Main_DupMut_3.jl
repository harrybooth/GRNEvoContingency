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

##########

grn_parameters = DefaultGRNParameters();

max_w = 10.

output_gene = 3

n_stripe = 1

min_width = 5

lower_bound = 5.

upper_bound = 10.

max_conc = 20.

fitness_function = s -> fitness_evaluation(s,x->(nstripe_fitness(x,n_stripe,min_width,lower_bound,upper_bound),malt_fitness_absolute(x,n_stripe,max_conc)),output_gene);

# fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_stripe),output_gene);

tolerance = 0.9

###############

vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

viable_mutations = ones(Int,Ng,Ng+1)

viable_mutations[2,4] = 0
viable_mutations[3,4] = 0

mutation_weights = findall(viable_mutations.> 0)

weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in Tuple.(mutation_weights)]

mut_prob = 0.2
pm_prob = 0.5

reg_flip_probability = 0.1

##############

N = 10^5

μ = 1e-7 * N

min_affinity = 1e-5
max_affinity = 10. 

v = log(min_affinity /max_affinity)

Ls = 10

ϵ = -v/Ls

k = 4

#############

mean_m = mean([ϵ*((Ls-d0)*μ - d0*μ/(k-1)) for d0 in 1:Ls])
mean_v = mean([ϵ^2*((Ls-d0)*μ + d0*μ/(k-1)) for d0 in 1:Ls])

noise_distribution = Normal(mean_m,mean_v)

############

TF_A = TF("A",[CartesianIndex(1, 1),CartesianIndex(2, 1),CartesianIndex(3, 1)],mut_prob,pm_prob,noise_distribution,noise_distribution,min_affinity,reg_flip_probability)
TF_B = TF("B",[CartesianIndex(1, 2),CartesianIndex(2, 2),CartesianIndex(3, 2)],mut_prob,pm_prob,noise_distribution,noise_distribution,min_affinity,reg_flip_probability)
TF_C = TF("C",[CartesianIndex(1, 3),CartesianIndex(2, 3),CartesianIndex(3, 3)],mut_prob,pm_prob,noise_distribution,noise_distribution,min_affinity,reg_flip_probability)

TFBS_list = [TFBS(name,ci,mut_prob,pm_prob,noise_distribution,noise_distribution,min_affinity,reg_flip_probability) for (name,ci) in zip(weight_names,mutation_weights)];

all_sites = vcat([TF_A,TF_B,TF_C],TFBS_list);

n_sample_func() = rand(Binomial(length(all_sites),mut_prob))

mutate_function = i -> noise_mtype_dup(i,all_sites,n_sample_func)

#############

β = (1.,N)

n_trials = 10000
max_gen = 250000

