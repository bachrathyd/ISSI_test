
using DifferentialEquations
using LinearAlgebra
using Random

using Plots

using StaticArrays

using BenchmarkTools

using Statistics
using Interpolations


using Revise
#push!(LOAD_PATH,"C:/Users/Bachrathy/.julia/dev/ISSI_test/sandbox")
push!(LOAD_PATH,"C:/Users/Mechanics/Documents/Bachrathy/Git/ISSI_test/sandbox")
using SimulationBasedStability


function delay_mathieu_model(du,u,h,p,t)
    δ, ϵ, b, κ, τ = p
    du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ)[2] 
    du[2] = u[1]
end
    
δ=0.50
ϵ=0.60
b=0.07
κ=0.2
τ=2pi 
T=2pi

p=( δ, ϵ, b, κ, τ )
h(p, t) = [1.0,1.0im] 
lags = [τ]
taumax=maximum(lags)
u0=[1.0+1.0im,1.0+1.0im]

alg = MethodOfSteps(BS3());

tspan = (0.0, T)
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags);

#-----------------------------------------------------------------------------
numerical_integraion=solve(prob,alg);
plot(numerical_integraion)
#-----------------------------------------------------------------------------
DEsolution=spectralRadiusOfMapping(prob,alg,taumax);
DEsolution.eigs #eigen values
plot(DEsolution.SolutionSet[1]) #"continuouse" soltuion of eigen vectors 


