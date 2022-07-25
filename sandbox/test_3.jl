
#import Pkg


#using LinearAlgebra
using DifferentialEquations
using StaticArrays
#using Random
using Plots

#using Statistics
#using Interpolations

using BenchmarkTools

using Statistics
using Interpolations


using Revise
#include("SimulationBasedStability.jl")
push!(LOAD_PATH,"C:/Users/Bacharthy/.julia/dev/ISSI_test/sandbox")
using SimulationBasedStability
#include("structures.jl")
#include("functinos.jl")

function delay_mathieu_model!(du,u,h,p,t)
    δ, ϵ, b, κ, τ = p
    
    #u2_past = h(p, t-τ; idxs=2)
    #du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * u2_past+0.0*cos(5.0*t)
    du[1] = - κ * u[1] - (δ + ϵ * sign(cos(t))) * u[2] + b * h(p, t-τ)[2]+0.0*cos(5.0*t)
    du[2] = u[1]
    nothing
end


δ=3.0
#ϵ=2.95
ϵ=0.8
b=-0.05
κ=0.01
τ=2pi 
T=1.0*2pi

p=( δ, ϵ, b, κ, τ )

h(p, t) = [0.0+0.0im,0.0+0.0im] 
h(p,t) = CompRand([1.0,1.0])
#h(p, t; idxs=nothing) = typeof(idxs) <: Number ? CompRand(1.0) : CompRand([1.0,1.0])

lags = [τ]
taumax=maximum(lags)
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

u0=[1.0+1.0im,1.0+1.0im];
#u0=[0.0,0.0];
prob = DDEProblem(delay_mathieu_model!, u0, h, tspan, p; constant_lags = lags,reltol=1e-2,abstol=1e-2,dtmax=T/25);#

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
aaa=solve(prob,alg);
plot(aaa)
@benchmark solve(prob,alg)

# --------------------------------------------

#@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=10,eigN=4);
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,eigN=4);
@time spectralRadiusOfMapping!(DelayMathieu_egi_prblem);


@benchmark  DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,eigN=4)
@benchmark  spectralRadiusOfMapping(DelayMathieu_egi_prblem) setup=(DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,eigN=4))

plot(DelayMathieu_egi_prblem.SolutionSet[1])

plot(DelayMathieu_egi_prblem.SolutionSet[1].t)


# -----------------root iteration -----------------------------
EigN=10
@time DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,eigN=EigN,Historyresolution=500);
NN=15;
Spect=zeros(EigN,NN)
timeing=zeros(NN)
plot()
for k=1:NN
   # k=0
  #  k+=1
  timeing[k]=@elapsed begin
        iterate!(DelayMathieu_egi_prblem);
        ei,si,vi,aii=compute_eig!(DelayMathieu_egi_prblem);
  end 
  @show timeing[k]
Aμs,Ai=eigen(aii, sortby = x -> -abs(x))
#@show Ai_norm=norm(abs.(Aμs) .-1 )
@show Spect[:,k]=abs.(ei)
#scatter(ones(DelayMathieu_egi_prblem.eigN) .*k ,log.(abs.(ei)))
scatter!(ones(DelayMathieu_egi_prblem.eigN) .*k ,log.(abs.(ei)))
end
plot!()
plot(timeing)
plot(log.(abs.(Spect[:,1:NN ] .- Spect[:,end]))')
#------------------------- Profiling ---------------

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

#  ----------------------MDBM test--------------------------

function foo(x)
    ploc=(x[1],x[2],p[3:end]...);

    pr = DDEProblem(delay_mathieu_model!, u0, h, tspan, ploc; constant_lags = lags, reltol=1e-3,abstol=1e-3);
    
    DM=dynamic_problem(pr,MethodOfSteps(Tsit5()),taumax,Historyresolution=10,eigN=2);
    spectralRadiusOfMapping!(DM)
    return abs.(DM.eigs[1])-1
end

function foo(x,y)
    foo([x,y])
end



using MDBM

using Plots
gr();

axis=[Axis(-1:0.4:5.,:δ),
    Axis(-2:0.4:1.5,:ϵ)]

iteration=4;#0- csak a kezdeti háló,1,2
mdbmsol=MDBM.solve!(MDBM_Problem(foo,axis),0);
for k=1:iteration
   
@time mdbmsol=MDBM.solve!(mdbmsol,1);
stab_border_points=getinterpolatedsolution(mdbmsol);

scatter(stab_border_points...,xlim=(-1.,5),ylim=(-2.,1.5),
    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)

end
