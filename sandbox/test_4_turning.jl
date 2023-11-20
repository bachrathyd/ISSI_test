
using DifferentialEquations

using Plots

using BenchmarkTools

using Statistics
using Interpolations
using LinearAlgebra

using Revise
push!(LOAD_PATH,"C:/Users/Bacharthy/.julia/dev/ISSI_test/sandbox")
using SimulationBasedStability


function turning_modell!(du,u,h,p,t)
    Ω,kw,κ = p
    du[1] = - 2.0 * κ * u[1] - (1+kw) * u[2] -kw * h(p, t-2*pi/Ω)[2]
    du[2] = u[1]
    nothing
end



Ω=0.34;
kw=0.1;
κ =0.2;

p=(Ω,kw,κ)
h(p,t) = [1.0,1.0]
τ=2π/Ω
lags = [τ]
taumax=maximum(lags)
T=τ;#0.1*τ;0.1
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

u0=[1.0,0.0];

prob = DDEProblem(turning_modell!, u0, h, tspan, p; constant_lags = lags,reltol=1e-2,dtmax=T/5);#

alg = MethodOfSteps(Tsit5());

just_a_test_solution=solve(prob,alg);
plot(just_a_test_solution)





@time DM=dynamic_problem(prob,MethodOfSteps(Tsit5()),taumax,Historyresolution=10,eigN=4);
compute_eig!(DM)
iterate!(DM)
@time spectralRadiusOfMapping(DM)

# ---------------------

function foo(x)
    println(x)
    
    pr=remake(prob;p=ploc=(x[1],x[2],0.1...),tspan=[0,T],constant_lags=[2π/x[1]])
    DM=dynamic_problem(pr,MethodOfSteps(Tsit5()),taumax,Historyresolution=10,eigN=4);
    
    return spectralRadiusOfMapping(DM)[1]-1
end

function foo(x,y)
    foo([x,y])
end



using MDBM

using Plots
gr();

axis=[Axis(0.2:0.2:2.0,:Ω),
    Axis(-0.5:0.2:0.5,:kw)]

iteration=2;#0- csak a kezdeti háló,1,2
mdbmsol=MDBM.solve!(MDBM_Problem(foo,axis),0);

for k=1:iteration
   
@time mdbmsol=MDBM.solve!(mdbmsol,1);
stab_border_points=getinterpolatedsolution(mdbmsol);

scatter(stab_border_points...,xlim=(-1.,5),ylim=(-2.,1.5),
    guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)

end
