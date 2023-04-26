########## GRN Model Setup ######### 

const Nc = 100
const Ng = 3
const L = 1.

const m0 = 10.
const λm = 0.4

const h_a = -1.
const h_b = 0.1
const deg_rate_g = 0.05
const init_conc_g = 0.01

const tissue = range(0,L,length = Nc)

morph(x) = m0*exp(-x/λm)

σ(I) = 0.5*(((I + h_a)/sqrt((I + h_a)^2+h_b)) + 1) # σ(0.) > 0 ?

include(srcdir("TissueModel_ND.jl"))

########## Evolutionary Setup ######### 

β = 1.

noise_cv = 0.25

mut_prob = 0.1

deletion_prob = 0.05

grn_parameters = DefaultGRNParameters();

viable_mutations = ones(Int,Ng,Ng+1)

mutation_weights = findall(viable_mutations .> 0)

n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

max_w = 1.

mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),n_sample_func,deletion_prob,max_w,mutation_weights)

mutate_function = i -> noise(i,mutation_op);

output_gene = 3

fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_target_stripe),output_gene);

tolerance = 0.9

######### Simulation setup ######### 

n_traj = 20
max_gen = 10000

n_target_stripe = 1

topology = "bistable"

if topology == "feed_forward"
    start_network = [0.0 0.0 0.0 0.28368795845354794; 0.09693796878733349 0.0 0.0 0.0; 0.02660150950444218 -0.26272166357617865 0.6146272196396064 0.0] # right handed feed forward
elseif topology == "bistable"
    start_network = [0.0 0.0 0.0 0.12728709721871537; -0.014075930837614938 0.0 0.49938778625866675 0.0; 0.0 0.1997636901215515 0.05755668522756788 0.0] # right handed bistable
elseif topology == "mutual_inh"
    start_network = [0.0 0.0 0.0 -0.06486368943640441; -1.0 0.0 0.7990661288117087 0.0; 0.8329310528122276 0.05802424255402265 0.0 0.0]  # right handed mutual_inh
elseif topology == "classical"
    start_network = [0.0 0.0 0.0 0.31336321352475677; 0.034733909122486334 0.011312260670141676 0.0 0.0; -0.002116043070043973 -0.1751097063592794 0.5246935289524377 0.0] # right handed classical
else
    start_network = [0.0 0.0 0.0 0.28368795845354794; 0.09693796878733349 0.0 0.0 0.0; 0.02660150950444218 -0.26272166357617865 0.6146272196396064 0.0] # right handed feed forward
end
