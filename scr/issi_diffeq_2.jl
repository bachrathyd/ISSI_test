using DifferentialEquations
using Interpolations
using LinearAlgebra
using PyPlot
pygui(true)


function historymaker(si,τ,N)
    #itph = interpolate(si, BSpline(Constant()))#Constant, Linear, Quadratic, and Cubic corresponding to B-splines of degree 0, 1, 2, and 3 respectively.
    itph = interpolate(si, BSpline(Linear()))#Constant, Linear, Quadratic, and Cubic corresponding to B-splines of degree 0, 1, 2, and 3 respectively.
    #itph = interpolate(si, BSpline(Quadratic()))#Constant, Linear, Quadratic, and Cubic corresponding to B-splines of degree 0, 1, 2, and 3 respectively.
    #itph = interpolate(si, BSpline(Cubic(Line(OnGrid()))))#Constant, Linear, Quadratic, and Cubic corresponding to B-splines of degree 0, 1, 2, and 3 respectively.
# sitp = scale(itp, A_x)
  h(p, t; idxs) =itph((1.0-(-t/τ))*(N-1)+1.0)
  # h([1],-0.1,1)
  # t=-6:0.01:0
  # (tt->h(1,tt)).(t)
end

function delay_mathieu_model(du,u,h,p,t)
  δ, ϵ, b, κ, τ = p
  # println(h(p, t-τ))
  # println(u[:,1])
  # println(u[:,2])
  # println(du[:,1])
  # println(du[:,2])
  # du[:,1] .= - κ .* u[:,1] - (δ + ϵ * cos(t)) .* u[:,2] - b .* h(p, t-τ)
  # du[:,2] .= u[:,1]
  du[1] = - κ * u[1] - (δ + ϵ * cos(t)) * u[2] + b * h(p, t-τ, idxs=2)
  du[2] = u[1]
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
delayoscillmapp(2.0, 0.02)



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
