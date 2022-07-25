
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


Neval=0    
function delay_mathieu_model(du,u,h,p,t)
    global Neval
    Neval += 1;
    δ, ϵ, b, κ, τ = p
    du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ)[2] + 1.0* 0.1*cos(2.0 *t)
    du[2] = u[1]
end

    
δ=0.50
ϵ=0.60#0.003
b=0.07#0.05
κ=0.2
τ=2pi 
T=2pi

p=( δ, ϵ, b, κ, τ )
h(p, t) = [0.0,0.0] .* (1.0+1.0im)
lags = [τ]
taumax=maximum(lags)
u0=[1.0+1.0im,1.0+1.0im]

#for Nmax in 60:100.0
Nmax=60
tspan = (0.0, Nmax*T)
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags,reltol=1e-5,abstol=1e-5);#

#alg = MethodOfSteps(Tsit5());
alg = MethodOfSteps(BS3());

@time just_a_test_solution=solve(prob,alg);
plot(just_a_test_solution)

@show Neval    

#periodic only
finalind= (just_a_test_solution.t .-(Nmax-1)*T) .> 0;

tlast=(just_a_test_solution.t[finalind] .-(Nmax-1)*T);
ulast=(just_a_test_solution.u[finalind]);

u1=  [real.(u[1]) for u in ulast]
plot(tlast,u1)
#end


#=

tspan = (0.0, T)
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags);#
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=50,eigN=8,zerofixpont=true,maxiteration=6);
spectralRadiusOfMapping!(DelayMathieu_egi_prblem);
abs.(DelayMathieu_egi_prblem.eigs)



δ=0.50
ϵ=0.60#0.003
b=0.07#0.05
κ=0.2
τ=2pi 
T=2pi
0.8918815190927252
 0.8918814676421285
 0.10121518355160673
 0.041138727803797015
 0.0411386468377296
 0.012116777236042138
 0.01211635335709965
 0.0047864544574736665
 =#
#----------------------------
tspan = (0.0, T)
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags);#
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=50,eigN=8,zerofixpont=false,maxiteration=0);

iterate!(DelayMathieu_egi_prblem);
iterate!(DelayMathieu_egi_prblem);
iterate!(DelayMathieu_egi_prblem);
iterate!(DelayMathieu_egi_prblem);
compute_eig!(DelayMathieu_egi_prblem);
@show abs.(DelayMathieu_egi_prblem.eigs)

DelayMathieu_egi_prblem.Ai
norm(DelayMathieu_egi_prblem.Ai[:,1])
sum(DelayMathieu_egi_prblem.Ai[:,1])


dp=DelayMathieu_egi_prblem;
abs.(dp.eigs)
SS,VV=SVi1real(dp::dynamic_problem,1)

k=1
#plot!(DelayMathieu_egi_prblem.SolutionSet[k])
plot(dp.StateSmaplingTime,real.(SS[:,k]) .+ 0.0)
plot!(dp.StateSmaplingTime .+ 2*pi,real.(VV[:,k])-real.(SS[:,1]) .* 0.0)
plot!(tlast .- 2pi,(u1))
plot!(tlast,(u1))


k=0
k+=1
plot!(tlast,(u1))
plot(DelayMathieu_egi_prblem.SolutionSet[k])


#DelayMathieu_egi_prblem.Ai
#DelayMathieu_egi_prblem.StateCombinations
#DelayMathieu_egi_prblem.Si'*DelayMathieu_egi_prblem.Si









#----------------------------Repell test-------------------

prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags,reltol=1e-3,abstol=1e-3);#
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=50,eigN=12,zerofixpont=false,maxiteration=1);



iterate!(dp);
compute_eig!(dp);


dp=DelayMathieu_egi_prblem;
Tmap=dp.DDEdynProblem.tspan[end]::Float64

dp.Si .= [getvalues(dp.SolutionSet[solind],t) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];
dp.Vi .= [getvalues(dp.SolutionSet[solind],(t+Tmap)) for t in dp.StateSmaplingTime, solind in 1:dp.eigN];#initialization of the Starting Vector
    
println("fixpoint test 1")
            
ϵP=0.01
AP=zeros(Float64,dp.eigN,dp.eigN)
AP[1,:] .= 1.0
AP[2:end,2:end].=ϵP*diagm(0 => ones(dp.eigN-1));
APinv=inv(AP)
dp.Si .= dp.Si * APinv;
dp.Vi .=  dp.Vi * APinv;
Hi=(dp.Si'*dp.Si)\(dp.Si'*dp.Vi)

EigSol=eigen((dp.Si'*dp.Si)\(dp.Si'*dp.Vi), sortby = x -> -abs(x));

dp.eigs .= EigSol.values;
dp.Ai .= EigSol.vectors;
dp.StateCombinations[:] .= (APinv*(dp.Ai/ diagm(dp.eigs)) *AP )[:];#length(dp.StateSmaplingTime)*  ./ Snorm   ./ μs

dp.StateCombinations[:] .= ((dp.Ai) )[:];


iterate!(dp);
compute_eig!(dp);



dp=DelayMathieu_egi_prblem;

SS,VV=SVi1real(dp::dynamic_problem,1)

k=1
#plot!(DelayMathieu_egi_prblem.SolutionSet[k])
plot(dp.StateSmaplingTime,real.(SS[:,k]) .+ 0.0)
plot!(dp.StateSmaplingTime .+ 2*pi,real.(VV[:,k]) .+ 0.0)
plot!(tlast,(u1))

#for k=2:DelayMathieu_egi_prblem.eigN-1
#plot!(DelayMathieu_egi_prblem.SolutionSet[k])
#end
#plot!(DelayMathieu_egi_prblem.SolutionSet[DelayMathieu_egi_prblem.eigN])




k=0
k+=1
plot!(tlast,(u1))
plot!(DelayMathieu_egi_prblem.SolutionSet[k])
