
import Pkg


using DifferentialEquations
using LinearAlgebra
using Random
using Plots

using StaticArrays

using BenchmarkTools

using Statistics
using Interpolations


using Revise
include("SimulationBasedStability.jl")
include("structures.jl")
#include("functinos.jl")

function delay_mathieu_model(du,u,h,p,t)
    δ, ϵ, b, κ, τ = p
    du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ)[2]
    du[2] = u[1]
end



δ=3.0
ϵ=0.8
b=-0.05
κ=0.01
τ=2pi 
T=2pi

p=( δ, ϵ, b, κ, τ )
#h(p, t; idxs) = (rand()-0.5) .+ 1im*(rand()-0.5)*10.0
h(p, t) = [0.0,0.0] .* (1.0+1.0im)
lags = [τ]
taumax=maximum(lags)
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

u0=[1.0+1.0im,1.0+1.0im]
#u0=[0.0,0.0];
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags,reltol=1e-1,abstol=1e-1,dtmax=T/5);#
#TODO: kell ez egyáltalán? , dtmax=T/10.0
#TODO:  saveat - testing

#TODO: , dense=false - ettől kétszer gyorsabb, de mintha más lenne a gyök utána!!! - ettől nem lesz meg a continuouse interpolation!!!
#Mert lineáris interpolálást csinál!!!!

#TODO: solve(..., alg_hints=[:stiff])
alg = MethodOfSteps(Tsit5());
#alg = MethodOfSteps(BS3());
#alg = MethodOfSteps(Vern6()); XXX
#alg = MethodOfSteps(Rosenbrock23()); #stiff
#alg = MethodOfSteps(Rodas4()); #stiff
#@btime
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=10,eigN=4);
@time ei,eis=spectralRadiusOfMapping(DelayMathieu_egi_prblem);

plot(DelayMathieu_egi_prblem.SolutionSet[1])

for k=1:20
    k=0
    k+=1
iterate!(DelayMathieu_egi_prblem);
ei,si,vi,aii=compute_eig!(DelayMathieu_egi_prblem);
Aμs,Ai=eigen(aii, sortby = x -> -abs(x))
@show Ai_norm=norm(abs.(Aμs) .-1 )
@show abs.(ei)
#scatter(ones(DelayMathieu_egi_prblem.eigN) .*k ,log.(abs.(ei)))
scatter!(ones(DelayMathieu_egi_prblem.eigN) .*k ,log.(abs.(ei)))
end

function profile_test(n)
    for i = 1:n
    DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=20,eigN=4);
    ei,eis=spectralRadiusOfMapping(DelayMathieu_egi_prblem);
end
end

# compilation
@profview profile_test(1)
# pure runtime
@profview profile_test(10)

#  ===================================================

#=
dfhj
spectr
gzu
=#
function foo(x)
    ploc=(x[1],x[2],p[3:end]...);

    pr = DDEProblem(delay_mathieu_model, u0, h, tspan, ploc; constant_lags = lags, reltol=1e-3,abstol=1e-3);
    
    DM=dynamic_problem(pr,MethodOfSteps(Tsit5()),taumax,Historyresolution=10,eigN=2);
    
    return spectralRadiusOfMapping(DM)[1]-1
end

function foo(x,y)
    foo([x,y])
end

#δi=-1.0:0.1:3.0
#ϵi=-1.0:0.1:0.8

#SpecR=[foo([d,e]) for d in δi, e in ϵi]

#plot(δi,ϵi,log.(SpecR)')
#contour(δi,ϵi,log.(SpecR)',levels=-2.0:0.1:0.0,fill=true)
#contour!(δi,ϵi,log.(SpecR)',levels=0.0:0.2:2.0,fill=false)








#-------------------------------------------

using MDBM

using Plots
gr();

axis=[Axis(-1:0.4:5.,:δ),
    Axis(-2:0.4:1.5,:ϵ)]

iteration=4;#0- csak a kezdeti háló,1,2
mdbmsol=MDBM.solve!(MDBM_Problem(foo,axis),0);
for k=1:iteration
   
end
