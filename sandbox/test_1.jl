
import Pkg


using DifferentialEquations
using LinearAlgebra
using Random
using Plots


using Statistics
#using Interpolations



using Revise
include("SimulationBasedStability.jl")
include("structures.jl")
#include("functinos.jl")


function delay_mathieu_model(du,u,h,p,t)
    δ, ϵ, b, κ, τ = p
  
    du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ)[2]
    du[2] = u[1]
end



#h(p, t; idxs) = (rand()-0.5) .+ 1im*(rand()-0.5)*10.0
h(p, t) = [1.0,1.0]# .* (1.0+1.0im)
δ=3.0
ϵ=2.0
b=-0.15
κ=0.1
τ=2pi 
T=2pi
p=( δ, ϵ, b, κ, τ )
lags = [τ]
taumax=maximum(lags)
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

#u0=[1.0+1.0im,1.0+1.0im]
u0=[1.0,1.0];
prob = DDEProblem(delay_mathieu_model, u0, h, tspan, p; constant_lags = lags);
alg = MethodOfSteps(Tsit5());

#abc=solve(prob,alg,reltol=1e-3,abstol=1e-3)


DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=100,eigN=2);


compute_eig!(DelayMathieu_egi_prblem);

for k=1:10
iterate!(DelayMathieu_egi_prblem);
ei=compute_eig!(DelayMathieu_egi_prblem)
@show maximum(abs.(ei.values))
end

#alg=Rosenbrock23(autodiff=false)
#alg=Rodas4(autodiff=false)

#big_u0=big.(u0)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

DelayMathieu_egi_prblem=dynamic_problem(prob,alg,taumax,Historyresolution=100,eigN=8);


Nshift=0;
Si,Vi =SVi1real(DelayMathieu_egi_prblem,1)
plot(DelayMathieu_egi_prblem.StateSmaplingTime.+ T*(Nshift-1),real.(Si),linewidth=2)
plot!(DelayMathieu_egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(Vi),linewidth=2)

ei=compute_eig!(DelayMathieu_egi_prblem)
@show (abs.(ei))
#DelayMathieu_egi_prblem.StateCombinations[:] .= [1.0,-10.0,0.0,1.0]

A =DelayMathieu_egi_prblem.StateCombinations
plot!(DelayMathieu_egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(Vi*A),linewidth=5)


Nshift += 1;
iterate!(DelayMathieu_egi_prblem);
Si,Vi =SVi1real(DelayMathieu_egi_prblem,1) #Si2,Vi2
plot!(DelayMathieu_egi_prblem.StateSmaplingTime.+ T*(Nshift-1),real.(Si),linewidth=2)
plot!(DelayMathieu_egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(Vi),linewidth=2)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Vi*DelayMathieu_egi_prblem.StateCombinations .- Si2




A' * diagm(μs) * A
-Hi

Snorm=[norm(Si[:,k]) for k in 1:DelayMathieu_egi_prblem.eigN ];



solset=DelayMathieu_egi_prblem.SolutionSet;
plot(solset[1])
plot!.(solset[2:end-1])
plot!(solset[end])


#<<<<<<<<<<<<<<<<<<<<<<<<<<<
plot()
for kk=1:6
s=DelayMathieu_egi_prblem.SolutionSet[kk]
ti=collect(-taumax:0.01:T)
SS=getvalues(s,ti)
Nshift = 0;
plot!(ti .+ Nshift*T,real.(getvalues(s,ti)))
end

iterate!(DelayMathieu_egi_prblem);
ei=compute_eig!(DelayMathieu_egi_prblem)
#@show ei.values
@show maximum(abs.(ei.values))
Nshift += 1;
for kk=1:6
s=DelayMathieu_egi_prblem.SolutionSet[kk]
SS=getvalues(s,ti)
plot!(ti .+ Nshift*T,real.(getvalues(s,ti)))
end
plot!()

#<<<<<<<<<<<<<<<<<<<<<<<<<<<

SV=Si'*Vi;
SS=Si'*Si;

Snorm=[norm(Si[:,k]) for k in 1:dp.eigN ];


solset=DelayMathieu_egi_prblem.SolutionSet;
plot(solset[1])
plot!.(solset[2:end-1])
plot!(solset[end])

#using StaticArrays
#A = @SMatrix rand(3,2)
# f(u,p,t)=A*u

Ai=diagm(0 => ones(Float64,eigvecnum))

u0=1.0+.5im

rand(Complex{typeof(u0[1])}, 10)

Ntime=10
taumax=maximum(prob.constant_lags)
timelocations=-taumax  :taumax/(Ntime-1) :0

#Vi=[solset[solind](t+T) for t in timelocations, solind in 1:eigvecnum];#initialization of the Starting Vector
#------------------------------

for kiter in 1:25
#Si=Vi*Ai;
solsetnew=deepcopy(solset);
for k in 1:eigvecnum 
    #h(p, t; idxs) = histopryremapp(t+T,Ai,solset,k)[idxs]#+(rand(size(idxs)[1]) .+ 1im*rand(size(idxs)[1]) .- 0.5 .- 0.5im)*10.0 # overwrites the original function even inside the prob.variable
    h(p, t;) = histopryremapp(t+T,Ai,solset,k)
    prob.u0 .= h(p, 0)
   solsetnew[k] = solve(prob,alg,reltol=1e-3,abstol=1e-3);
end
solset=solsetnew;
Si=[solset[solind](t) for t in timelocations, solind in 1:eigvecnum];

Vi=[solset[solind](t+T) for t in timelocations, solind in 1:eigvecnum];#initialization of the Starting Vector







SV=Si'*Vi;
SS=Si'*Si;

#VV=Vi'*Vi;
#SV=Si'*Vi;

Snorm=[norm(Si[:,k]) for k in 1:eigvecnum ];

Hi=SS\SV
#Hi=SV\VV
#Si\Vi
Eigs=eigen(Hi);
@show μs=Eigs.values
Ai=Eigs.vectors  ./ Snorm ./ Eigs.values;
end
plot(solset[1])
plot!.(solset[2:end-1])
plot!(solset[end])

#plot(1:eigvecnum,1:eigvecnum,real.(Ai))

#------------------------------


Ntime=100
taumax=maximum(prob.constant_lags)
timelocations=-taumax  :taumax/(Ntime-1) :0

Vi=[solset[solind](t) for t in timelocations, solind in 1:eigvecnum];

plot(solset[1].t)

















    # figure(9)
    # clf()
    # surf(abs.(Hi))

    Si=Vi*Eigs.vectors
    for iloc in 1:eigvecnum
        sinorminv=1/norm(Si[:,iloc])
        Si[:,iloc] *= sinorminv
    end
    push!(eigsol_inter,Eigs.values)
println("---------------------------------")




solset = [solve(prob,alg) for k in 1:eigvecnum]

plot(solset[1])
plot!.(solset[2:end-1])
plot!(solset[end])

size(solset)


@show Si=[sol(t+T) for t in timelocations];

solset[1](2)



@show Si=[sol(t+T) for t in timelocations];


sum([Si'[i]*Si[i] for i in 1:Ntime])
Si'*Si
plot(timelocations,Si)
Si'[end]*Si[end]
# initialization of the Eigen vector
N = 50
eigvecnum = 6
# δ, ϵ, b, κ, τ = [2.0,1.0,0.02,0.0,2*pi] # [1.0,2.0,1,0.1]
# p=[δ, ϵ, b, κ, τ]
p=[δ,1.0,b,0.0,2*pi]
#
Si = rand(N, eigvecnum)
eigsol_inter=Array{Array{ComplexF64}}(undef,0)
#@time begin
for k in 1:6
    #global(Si)
    Vi=mappinggeneration(Si,p)
    Hi=Si\Vi
    # figure(9)
    # clf()
    # surf(abs.(Hi))
    Eigs=eigen(Hi)
    Si=Vi*Eigs.vectors
    for iloc in 1:eigvecnum
        sinorminv=1/norm(Si[:,iloc])
        Si[:,iloc] *= sinorminv
    end
    push!(eigsol_inter,Eigs.values)

    # figure(14)
    # plot(k*ones(size(Eigs.values)),abs.(Eigs.values), linestyle="none", marker="o",color="green", markerfacecolor="blue", markersize=4)

    #figure(4)
    # clf()
    #plot(real(Eigs.values),imag(Eigs.values), linestyle="none", marker="o",color="green", markerfacecolor="blue", markersize=4)
    # fi=0:0.01:2*pi
    # plot(sin.(fi),cos.(fi))
    #pause(0.0001)
end



















function mappinggeneration(Si,p)
    # figure(3)
    # clf()
    Vi = zero(Si)
    τ=p[5]
    N=size(Si,1)
    for i in 1:size(Si,2)

        hist = historymaker(Si[:,i],τ,N)

        # p = [1.0,2.0,1,0.1]

        T = 2 * pi
        lags = [τ]
        tspan = (0.0, T)
        Δt = 1e-6
        u0 = [(hist(p, 0, idxs = 2) - hist(p, -Δt, idxs = 2)) / Δt, hist(p, 0, idxs = 2)]
        prob = DDEProblem(delay_mathieu_model, u0, hist, tspan, p; constant_lags = lags)

        alg = MethodOfSteps(Tsit5())
        sol = solve(prob, alg,reltol=1e-3,abstol=1e-3)#,reltol=1e-8,abstol=1e-8)


        tplot = LinRange(T - τ, T, N)
        solplot = [sol(ttt)[2] for ttt in tplot]
        Vi[:,i] .= solplot

        # figure(3)
        # # clf()
        # tplot = LinRange(-τ, 0, N)
        # solplot = [hist(p, ttt, idxs = 1) for ttt in tplot]
        # plot(tplot, real.(solplot), marker="o",color="green", markerfacecolor="blue", markersize=4)
        #
        # # tplot = LinRange(0, T - τ, N)
        # # solplot = [sol(ttt)[1] for ttt in tplot]
        # # plot(tplot, real.(solplot), marker="o",color="green", markerfacecolor="blue", markersize=4)
        #
        # tplot = LinRange(T - τ, T, N)
        # solplot = [sol(ttt)[2] for ttt in tplot]
        # plot(tplot, real.(solplot), marker="o",color="green", markerfacecolor="blue", markersize=4)

    end


    return Vi

end

function delayoscillmapp(δ, b)

# initialization of the Eigen vector
N = 50
eigvecnum = 6
# δ, ϵ, b, κ, τ = [2.0,1.0,0.02,0.0,2*pi] # [1.0,2.0,1,0.1]
# p=[δ, ϵ, b, κ, τ]
p=[δ,1.0,b,0.0,2*pi]
#
Si = rand(N, eigvecnum)
eigsol_inter=Array{Array{ComplexF64}}(undef,0)
#@time begin
for k in 1:6
    #global(Si)
    Vi=mappinggeneration(Si,p)
    Hi=Si\Vi
    # figure(9)
    # clf()
    # surf(abs.(Hi))
    Eigs=eigen(Hi)
    Si=Vi*Eigs.vectors
    for iloc in 1:eigvecnum
        sinorminv=1/norm(Si[:,iloc])
        Si[:,iloc] *= sinorminv
    end
    push!(eigsol_inter,Eigs.values)

    # figure(14)
    # plot(k*ones(size(Eigs.values)),abs.(Eigs.values), linestyle="none", marker="o",color="green", markerfacecolor="blue", markersize=4)

    #figure(4)
    # clf()
    #plot(real(Eigs.values),imag(Eigs.values), linestyle="none", marker="o",color="green", markerfacecolor="blue", markersize=4)
    # fi=0:0.01:2*pi
    # plot(sin.(fi),cos.(fi))
    #pause(0.0001)
end

#end

# figure(16)
# clf()
# for keig in 1:min(10,eigvecnum)#1:2
#     scatter(1:length(eigsol_inter),log.(abs.([abs(e[keig])-abs(eigsol_inter[end][keig]) for e in eigsol_inter])) ./ log(10))
#     #scatter(1:length(eigsol_inter),log.(abs.([e[keig] for e in eigsol_inter])))
# end
#
#
# figure(4)
# clf()
# plot(real(eigsol_inter[end]),imag(eigsol_inter[end]), linestyle="none", marker="o",color="green", markerfacecolor="blue", markersize=4)
# fi=0:0.01:2*pi
# plot(sin.(fi),cos.(fi))
# pause(0.1)

return abs(eigsol_inter[end][1])-1

end





using MDBM

using PyPlot;
pygui(true);

mymdbm=MDBM_Problem(delayoscillmapp,[-1:0.25:9,-1.5:0.25:1.5])

for kinter in 1:4
MDBM.solve!(mymdbm,1)


x_eval,y_eval=getevaluatedpoints(mymdbm)

x_sol,y_sol=getinterpolatedsolution(mymdbm)

fig = figure(112);clf()
scatter(x_eval,y_eval,s=5)
scatter(x_sol,y_sol,s=5)
pause(0.01)
end
