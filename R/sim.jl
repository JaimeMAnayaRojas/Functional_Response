using Distributions
using DataFrames
using Random
using StatsBase
using Plots
# using Turing
using CSV
# using ReverseDiff
# Turing.setadbackend(:reversediff)

using MCMCChains
using LinearAlgebra

using Stan


# set the path to the cmdstan folder



# using CmdStan
# set_cmdstan_home!("/home/mao/.cmdstan/cmdstan-2.31.0/")

# Functional response model
function FR_Holling(R, a, h, b=2)
    return (a .* (R.^b)) ./ (1 .+ a .* h .* (R.^b))
end

## Funcitons for analyses


function HDI(samples; credible_mass=0.95)
	# Computes highest density interval from a sample of representative values,
	# estimated as the shortest credible interval
	# Takes Arguments posterior_samples (samples from posterior) and credible mass (normally .95)
	# Originally from https://stackoverflow.com/questions/22284502/highest-posterior-density-region-and-central-credible-region
	# Adapted to Julialang
	sorted_points = sort(samples)
	ciIdxInc = Int(ceil(credible_mass * length(sorted_points)))
	nCIs = length(sorted_points) - ciIdxInc
	ciWidth = repeat([0.0],nCIs)
	for i in range(1, stop=nCIs)
		ciWidth[i] = sorted_points[i + ciIdxInc] - sorted_points[i]
	end
	HDImin = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)]
	HDImax = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)+ciIdxInc]
	return([HDImin, HDImax])
end


function LOS(v, b = 0)
	return 100*length(findall(v .> b)) ./length(v)
end


x = 0:0.1:5
y = FR_Holling(x, 10, 0.5, 2)
# plot the functional response
plot(x, y, xlabel="Resource density", ylabel="Eaten", legend=:topleft, title="Functional response")
y = FR_Holling(x, 10, 0.5, 1)
plot!(x, y)
y = FR_Holling(x, 10, 0.1, 1)
plot!(x, y)



# Define the number of individuals per group
N = 1000

# Define the means of a and h
mean_aC = 10
mean_aI = mean_aC * 0.9
mean_ai = mean_aC

mean_hC = log(0.5)
mean_hI = mean_hC
mean_hi = mean_hC


# Define the covariance matrix for as and hs

# create a covariance matrix with 0.01 on the diagonal and 0.005 on the off-diagonal
sigma = fill(0.005, 2, 2)
sigma[diagind(sigma)] .= 0.01
sigma

#sigma = [0.01 0.005; 0.005 0.01]

# Generate the parameters for each group using a multivariate normal distribution for the mean as and hs

parametersC = rand(MvNormal([mean_aC, mean_hC], sigma), N) 
parametersI = rand(MvNormal([mean_aI, mean_hI], sigma), N)
parametersi = rand(MvNormal([mean_ai, mean_hi], sigma), N)


# Extract the a and h values from the generated parameters
a = vcat(parametersC[1, :], parametersI[1, :], parametersi[1, :])
h = exp.(vcat(parametersC[2, :], parametersI[2, :], parametersi[2, :]))

# Put parameters in a data frame
df = DataFrame(ID = 1:length(a), a = a, h = h)

# Simulate values for R
R = 0:25

# Apply the function to the data frame
df2 = DataFrame()

for j in 1:length(R)
    Eaten = FR_Holling.(fill(R[j], length(df.a)), df.a, df.h, 1.0)
    ID = df.ID
    Exp = repeat(["C", "I+", "I-"], inner=N)
    R2 = fill(R[j], length(ID))
    df_temp = DataFrame(ID=ID, Exp=Exp, R2=R2, Eaten=Eaten)
    append!(df2, df_temp)
end

# add random variation to the data

df2.Eaten = exp.(log.(df2.Eaten) .+ rand(Normal(0, 0.05), length(df2.Eaten)))

minimum(df2.Eaten)


# Calculate mean and sd for each Exp and R2 combination
gdf = groupby(df2, [:Exp, :R2])
df3 = combine(gdf, :Eaten => mean, :Eaten => std)

# Plot the the functional response for each group and R2

plot(df3.R2, df3.Eaten_mean, group=df3.Exp, yerr=df3.Eaten_std, xlabel="Resource density", ylabel="Eaten", legend=:topleft, title="Functional response")
ylims!(0, 4)



### Select data for testing

# put in function

select_sim_data = function(data, n)

    # which individuals are in each Exp group?

    df3 = combine(groupby(data, [:Exp, :ID]), nrow => :n)


    #  for each group select 100 random individuals


    cIDs = sample(filter(row -> row.Exp == "C", df3)[:,:ID], n, replace=false)
    iIDs = sample(filter(row -> row.Exp == "I-", df3)[:,:ID], n, replace=false)
    IIDs = sample(filter(row -> row.Exp == "I+", df3)[:,:ID], n, replace=false)

    # select 60 individuals

    sID = vcat(cIDs, iIDs, IIDs)

    # find individuals with the selected IDs
    df4 = filter(row -> row.ID in sID, data)

    # sort data by ID

    df4 = sort(df4, :ID)

    # Define the data
    # remove all values withere R2 == 0 (because the model is not defined for R2 == 0)
    #df4 = filter(row -> row.R2 != 0, df4)

    # make sure that the ID is an integer and they range from 1 to length(unique(df4.ID))

    for i in 1:length(unique(df4.ID))
        println(i)
        df4.ID[df4.ID .== unique(df4.ID)[i]] .= i
    end
    df4.ID


    # Keep only six densities
    df4 = filter(row -> row.R2 in [0, 1, 3, 6, 12, 25], df4)
    df4 = sort(df4, [:ID, :R2])


    # save df4 as a csv file
    

    return df4


end


df4 = select_sim_data(df2, 100)

CSV.write("data/100sim_data.csv", df4)
CSV.write("data/10sim_data.csv",select_sim_data(df2, 10) )
CSV.write("data/15sim_data.csv",select_sim_data(df2, 15) )
CSV.write("data/20sim_data.csv",select_sim_data(df2, 20) )
CSV.write("data/30sim_data.csv",select_sim_data(df2, 30) )

# plot the functional response for each selected individual, colord by group

scatter(df4.R2, df4.Eaten, group=df4.Exp, xlabel="Resource density", ylabel="Eaten", legend=:topleft, title="Functional response")
xlims!(0, 26)

using DynamicHMC


# Hamiltonian Monte Carlo (HMC) sampler parameters
ϵ = 0.05
τ = 10

# filter data set by removing all values with R2 == 0
df4 = filter(row -> row.R2 != 0, df4)


# Define the model usint Turing to estimate the parameters a and h for each individual



# @model function FR_Holling_model(Eaten, R2, ID, Exp)

#     I = ifelse.(Exp .== "I+", 1., 0)
#     im = ifelse.(Exp .== "I-", 1., 0)
#     # Define the priors for the parameters
#     a ~ Normal(0, 0.5)
#     h ~ Normal(0, 1)
#     σ ~ Exponential(1)
#     σ_ID ~ Exponential(1)
    
#     # Add parameters for individual variation
    
#     a_i ~ filldist(Normal(0, σ_ID), length(unique(ID)))
#     h_i ~ filldist(Normal(0, σ_ID), length(unique(ID)))

#     # parameters for the parasite effect
#     α_I ~ Normal(0, 0.5)
#     α_i ~ Normal(0, 0.5)

#     h_I ~ Normal(0, 0.5 ) 
#     h_im ~ Normal(0, 0.5 )


#     # Define the likelihood
#     for i in 1:length(Eaten)

#         A = exp(a + a_i[ID[i]] + α_I * I[1] + α_i * im[i])
#         H = exp(h + h_i[ID[i]] + h_I * I[1] + h_im * im[i])
#         FR = FR_Holling(R2[i], A, H, 1.0)
#         Eaten[i] ~ Normal(FR, σ)
#     end


# end

# model = FR_Holling_model(df4.Eaten, df4.R2, df4.ID, df4.Exp)
# Run the model
#chain = sample(model, NUTS(), MCMCThreads(), 1000, 4)

# summarize the results
# describe(chain)

# summarize
# summarystats(chain)

# get posteriors

# post = DataFrame(chain)


# run stan model

# define the model


# define the data

data = Dict("N" => length(df4.Eaten), 
    "R" => df4.R2, 
    "Eaten" => df4.Eaten, 
    "ID" => df4.ID, 
    "I" => ifelse.(df4.Exp .== "I+", 1, 0),
    "nI" => ifelse.(df4.Exp .== "I-", 1, 0),
    "Nid" => length(unique(df4.ID))
    
    )

# run the model

using Distributions
using DataFrames
using MonteCarloMeasurements, AxisKeys
using StanSample

# read file with the model

model = "// Functional response model
// We want to estiate the effects of parasites on it. 
// Expose and Infected = E_If
// Exposed and not infected = E_nIf
data{
  int N; 
  int Nid;
  vector[N] Eaten;
  int ID[N];
  int I[N];
  int nI[N];
  int R[N];
}
parameters{
  real a;
  // real b;
  real h;
  
  real a_I;
  real h_I;
  real a_i;
  real h_i;
 
  
  
  vector[Nid] a_Intercept;
// vector[Nid] b_Intercept;
  vector[Nid] h_Intercept;
  real<lower=0> sigma_ID;
  real sigma;
}

model{
  vector[N] FR;
  vector[N] A;
// vector[N] B;
  vector[N] H;
  
  sigma ~ cauchy( 0 , 1 );
  sigma_ID ~ cauchy( 0 , 1 );
  a_Intercept ~ normal( 0 , sigma_ID );
// b_Intercept ~ normal( 0 , sigma_ID );
  h_Intercept ~ normal( 0 , sigma_ID );

  a_I ~ normal( 0 , 0.5 );
  h_I ~ normal( 0 , 0.5 );
  
  a_i ~ normal( 0 , 0.5 );
  h_i ~ normal( 0 , 0.5 );

 
  for ( i in 1:N ) {
    A[i] = exp(a + a_I * I[i] + a_i * nI[i] + a_Intercept[ID[i]]);
    H[i] = exp(h + h_I * I[i] + h_i * nI[i] + h_Intercept[ID[i]]);

    FR[i] = A[i]*(R[i]) / (1 + A[i]*H[i]*(R[i]) );
  }
  
  Eaten ~ normal( FR , sigma );

 

}
";


ms = SampleModel("fr_mod", model);
using StanSample

# run stan model

fit = stan_sample(ms; data);


# assmble the estimated and true parameters in a data frame
using Statistics

est_aC = median(exp.(post.a))
est_aC95 = HDI(exp.(post.a))

est_aI = median(exp.(post.a + post.α_I))
est_aI95 = HDI(exp.(post.a + post.α_I))

est_ai = median(exp.(post.a + post.α_i))
est_ai95 = HDI(exp.(post.a + post.α_i))

est_hC = median(exp.(post.h))
est_hC95 = HDI(exp.(post.h))

est_hI = median(exp.(post.h + post.h_I))
est_hI95 = HDI(exp.(post.h + post.h_I))

est_hi = median(exp.(post.h + post.h_im))
est_hi95 = HDI(exp.(post.h + post.h_im))


true_aC = mean_aC

true_aI = mean_aI

true_ai = mean_ai

true_hC = exp(mean_hC)

true_hI = exp(mean_hI)

true_hi = exp(mean_hi)

10^(-310)

# put in a data frame

df6 = DataFrame(Parameter = ["aC", "aI", "ai", "hC", "hI", "hi"], 
    True = [true_aC, true_aI, true_ai, true_hC, true_hI, true_hi],
    Est = [est_aC, est_aI, est_ai, est_hC, est_hI, est_hi],
    Est95 = [est_aC95, est_aI95, est_ai95, est_hC95, est_hI95, est_hi95])

