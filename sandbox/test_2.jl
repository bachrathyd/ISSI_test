
import Pkg


using DifferentialEquations
using LinearAlgebra
using Random
using Plots


using Statistics
using Interpolations



using Revise
include("SimulationBasedStability.jl")
include("structures.jl")


function simpleONE(du,u,h,p,t)
    n,d,A,tau = p  
    du[1] = 0*n*A * exp(1.0im*n*t+d) - 0.0 * u[1]+ 0.2*h(p, t-tau)[1]
end


h(p, t) = [0.0] .* (0.0+0.0im)
n=1.0
d=0.0
A=1.0
tau=2.0pi
T=2.0pi

p=( n,d,A,tau)
lags = [tau]
taumax=maximum(lags)
tspan = (0.0, T) # The end of the integration time considert to be the timeperiod of the system.

u0=[1.0+1.0im]

probi = DDEProblem(simpleONE, u0, h, tspan, p; constant_lags = lags, dtmax=T/20.0);
alg = MethodOfSteps(Tsit5());

egi_prblem=dynamic_problem(probi,alg,taumax,Historyresolution=400,eigN=4,p=p);


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#plot(egi_prblem.SolutionSet[1])
#plot!(egi_prblem.SolutionSet[2])
#plot!(egi_prblem.SolutionSet[3])
#plot!(egi_prblem.SolutionSet[4])


Nshift=0;
si,vi =SVi1real(egi_prblem,1)
plot(egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(si),linewidth=2)
plot!(egi_prblem.StateSmaplingTime.+ T*(Nshift),imag.(si),linewidth=1)
plot!(egi_prblem.StateSmaplingTime.+ T*(Nshift+1),real.(vi),linewidth=2)
plot!(egi_prblem.StateSmaplingTime.+ T*(Nshift+1),imag.(vi),linewidth=1)

plot!(legend = :bottomleft)


ei,sii,vii,aii=compute_eig!(egi_prblem);
#@show (abs.(ei));



#egi_prblem.StateCombinations[:] .= [1.0,-10.0,0.0,1.0]

egi_prblem.StateCombinations[1,1]=1.0
egi_prblem.StateCombinations[3,1]=3.0
egi_prblem.StateCombinations[2,2]=0.0
egi_prblem.StateCombinations[3,2]=-1.0
egi_prblem.StateCombinations[3,3]=0.0
egi_prblem.StateCombinations[4,4]=0.0

egi_prblem.StateCombinations[1,1]=1.0
egi_prblem.StateCombinations[3,1]=0.0
egi_prblem.StateCombinations[2,2]=1.0
egi_prblem.StateCombinations[3,2]=0.0
egi_prblem.StateCombinations[3,3]=1.0
egi_prblem.StateCombinations[4,4]=1.0
A =egi_prblem.StateCombinations
plot!(egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(vi*A),linewidth=5)

Nshift += 1;

iterate!(egi_prblem);
si,vi =SVi1real(egi_prblem,1); #Si2,Vi2
plot(egi_prblem.StateSmaplingTime.+ T*(Nshift-1),real.(si),linewidth=3)
plot!(egi_prblem.StateSmaplingTime.+ T*(Nshift),real.(vi),linewidth=2)

#Sj=reshape(S[:,1:n+1,:],(2*(n+1),mult))
#Vj=reshape(S[:,end-n:end,:],(2*(n+1),mult))
#Hj=Sj\Vj
#(λj,Gj)=eigen(Hj)

#S[:,1:n+1,:] = (Vj*Gj)[:]



μs,Ai=eigen(A, sortby = x -> -abs(x));
Ai*diagm(μs)*inv(Ai)-A#ez az a felbontás, amit a DZ is csinált, ez rendben van!

vv=vi'*vi
sv=si'*vi
ss=si'*si

vv=vii'*vii
sv=sii'*vii
ss=sii'*sii
plot(real.(vi[:,2]))

Hi=ss\sv
Hi-si\vi


μs,Ai=eigen(Hi)

using GenericSchur
S = schur(Hi .+ 0im)
VR = eigvecs(S)
VL = eigvecs(S, left=true)
inv(VR)
μs=eigvals(Hi)

VL*diagm(μs)*VR
A' * diagm(μs) * A
-Hi

Snorm=[norm(si[:,k]) for k in 1:DelayMathieu_egi_prblem.eigN ];

