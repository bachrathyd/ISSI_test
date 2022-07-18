5+5


#import Pkg

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


aaa=0
function delay_mathieu_model(du,u,h,p,t)
    global aaa
    aaa += 1;
    δ, ϵ, b, κ, τ = p
    du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ)[2] + 1.0*cos(4*t)
    du[2] = u[1]
end


δ=0.50
ϵ=0.003
b=0.05
κ=0.1
τ=2pi 
T=4pi

p=( δ, ϵ, b, κ, τ )
h(p, t) = [0.0,0.0] .* (1.0+1.0im)
lags = [τ]
taumax=maximum(lags)
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

u0=[1.0+1.0im,1.0+1.0im]
#u0=[0.0,0.0];
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags,reltol=1e-2,abstol=1e-2,dtmax=T/20);#

alg = MethodOfSteps(Tsit5());

@time just_a_test_solution=solve(prob,alg);
plot(just_a_test_solution)
@show aaa
just_a_test_solution.t

#@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=50,eigN=11);
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=50,eigN=11,zerofixpont=false);
@time ei=spectralRadiusOfMapping(DelayMathieu_egi_prblem)
@show ei
ei-0.849519055

DelayMathieu_egi_prblem.Ai
DelayMathieu_egi_prblem.StateCombinations
DelayMathieu_egi_prblem.Si'*DelayMathieu_egi_prblem.Si

abs.(DelayMathieu_egi_prblem.eigs) .- 0.849519




@show ei;
plot(DelayMathieu_egi_prblem.SolutionSet[1])
for k=2:DelayMathieu_egi_prblem.eigN-1
plot!(DelayMathieu_egi_prblem.SolutionSet[k])
end
plot!(DelayMathieu_egi_prblem.SolutionSet[DelayMathieu_egi_prblem.eigN])
