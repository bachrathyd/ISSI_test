using LinearAlgebra
using BenchmarkTools
using Statistics
using PyPlot
# pygui(true);


function mathieu(z, zτ, t, κ, δ, ϵ, b, τ, ω)
        dzdt=zeros(ComplexF64,size(z))
                dzdt[1]=z[2]
                dzdt[2]=-(δ+(ϵ * cos(ω * t)))*z[1] - κ * z[2] + b * zτ[1]
                return dzdt
end


#integration

function Sinterp(Sloc::Array{ComplexF64,2},indf::Float64)
        i_sub=Int64(floor(indf))::Int64
        remi=(indf-i_sub)::Float64
        # if abs(remi)<(eps(typeof(remi))*5.0)
        #         return S[:,i_sub] .* (1-remi)
        # else
        return Sloc[:,i_sub] .* (1.0-remi) .+ Sloc[:,i_sub+1]*remi
        # end
end

function Sinterp(S::Array{ComplexF64,2},i::Integer)
        return S[:,i]
end

# S=rand(ComplexF64,1,n+1+m)
# @time Sinterp(S,4.3)
# @time Sinterp(S,4)
# @time S[:,4]
# figure(1)
# plot(1:n+1+m,S[1,:])
# v=1:0.1:n+1+m
# scatter(v,[Sinterp(S[:,:,1],i) for i in v])


function mathiestab( κ, δ, ϵ, b, τ, ω, n,m,mult,iter,RKtype::String)

dt=τ/n
S=rand(ComplexF64,2,n+1+m,mult)
λj=Array{ComplexF64}(undef,mult)
for kinter in 1:iter
        #Threads.@threads
        for mi in 1:mult
                for kt in n+1:n+m
                if RKtype=="EE"
                        k1=mathieu(S[:,kt,mi], Sinterp(S[:,:,mi],kt-n),(kt-(n+1))*dt, κ, δ, ϵ, b, τ, ω)
                        k=k1 ./ 1.0
                elseif RKtype=="RK2-fast"
                         k1=mathieu(S[:,kt,mi]      ,S[:,kt-n  ,mi],(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                         k2=mathieu(S[:,kt,mi]+k1*dt,S[:,kt-n+1,mi],(kt-(n+1)+1)*dt, κ, δ, ϵ, b, τ, ω)
                         k=(k1+k2) ./ 2.0
                elseif RKtype=="RK2"
                        # RK2
                        k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        k2=mathieu(S[:,kt,mi]+k1*dt,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+k2) ./ 2.0
                elseif RKtype=="RK3"
                        # RK3
                        k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        k2=mathieu(S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                        k3=mathieu(S[:,kt,mi]-k1*dt*1.0+k2*dt*2.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+4*k2+k3) ./ 6.0
                elseif RKtype=="RK4"
                         # RK4
                         k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                         k2=mathieu(S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                         k3=mathieu(S[:,kt,mi]+k2*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                         k4=mathieu(S[:,kt,mi]+k3*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                         k=(k1+2*k2+2*k3+k4) ./ 6.0
                else
                        error("nincs ilyen integrator")
                end
                S[:,kt+1,mi]=S[:,kt,mi]+dt*k;
                end
        end

        Sj=reshape(S[:,1:n+1,:],(2*(n+1),mult))
        Vj=reshape(S[:,end-n:end,:],(2*(n+1),mult))
        Hj=Sj\Vj
        (λj,Gj)=eigen(Hj)

        S[:,1:n+1,:] = (Vj*Gj)[:]

        # figure(3)
        # clf()
        # plot((1:n+1+m).*dt,S[2,1:n+1+m,:])
        # plot((1:n+1).*dt,S[2,1:n+1,:])
        # pause(0.2)
        #
        # figure(4)
        # scatter(kinter .* ones(mult),abs.(λj))
end


return λj
end




println("----------=======<<<<<<<<<<<<<################>>>>>>>>>>>>>>========----------")
κ=0.05
δ=5
ϵ=1.0
b=0.1
τ=2pi
ω=1
#
mult=20 # number of computed eigenvalues
iter=30 # number of iteration in issi


# @debug
# @elapsed mu_ref=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω,10^6,m,mult,iter,"RK4")))
@time  mu_ref=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω, 2^4, Int64(2^4*(1 / ω)),mult,iter,"RK4")))
@time  mu_ref=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω, 2^4, Int64(2^4*(1 / ω)),mult,iter,"RK4")))


n=2^7#16
m=Int64(n*(1 / ω))
@time mu_ref=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω, n,m,mult,iter,"RK4")))

integrators=["EE","RK2","RK3","RK4"]
ni=2 .^ (3:6)
mumax=zeros(length(ni),length(integrators))
timelength=zeros(length(ni),length(integrators))
timestd=zeros(length(ni),length(integrators))


figure(1)
# clf()
ax=gca()
ax.set_xscale("log")
ax.set_yscale("log")
grid("on")
title("n / error")
xlabel("n")
ylabel("error")
#
figure(2)
# clf()
ax=gca()
ax.set_xscale("log")
ax.set_yscale("log")
grid("on")
title("n / time")
xlabel("n")
ylabel("time")
#
figure(3)
# clf()
ax=gca()
ax.set_xscale("log")
ax.set_yscale("log")
grid("on")
title("time / error")
xlabel("time")
ylabel("error")

Nsample=10

for kn in 1:length(ni)


        Threads.@threads for ksample in 1:Nsample
                println("<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>")
                println([ksample,kn,ni[kn]])
                for kRKi in 1:length(integrators)
                        integratortype=integrators[kRKi]
                        print(integratortype)
                        print("  -  ")
        n=ni[kn]
        m=Int64(n*(1 / ω))
        # timeing = @benchmark (mumax[kn]=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω,n,m,mult,iter,"RK4"))))
        # timelength[kn]=mean(timeing.times)
        # timestd[kn]=std(timeing.times)
        timeing = @elapsed (mumax[kn,kRKi]=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω,n,m,mult,iter,integratortype))))
        timelength[kn,kRKi]+=timeing
        println(timeing)


        figure(1)
        scatter(n,abs(mumax[kn,kRKi] - mu_ref))
        figure(2)
        scatter(n,timeing)
        figure(3)
        scatter(timeing,abs(mumax[kn,kRKi] - mu_ref))
                end
        end

figure(1)
plot(ni[1:kn,:],abs.(mumax[1:kn,:] .- mu_ref))
savefig("1_nt")
# savefig("1_t")

figure(2)
plot(ni[1:kn,:],timelength[1:kn,:]./Nsample)
savefig("2_nt")
# savefig("2_t")

figure(3)
plot(timelength[1:kn,:]./Nsample,abs.(mumax[1:kn,:] .- mu_ref))
savefig("3_nt")
# savefig("3_t")
end




println("-------------------")
                # figure(1)
                # plot(log.(ni)./log(2),log.(abs.(mumax .- mu_ref))./log(2))
                # figure(2)
                # plot(log.(ni)./log(2),log.(timelength)./log(2))
                # # plot(log.(timelength .+ timestd))
                # # plot(log.(timelength .- timestd))


# #----------------------------
# using MDBM
#
# foo=(x,y)->(maximum(abs.(mathiestab( κ, x, ϵ, y, τ, ω,n,m,mult,iter,"RK4")))-1)
#
# δv=-0.1:0.2:10
# # ϵv=-3:0.1:3
# bv=-1:0.2:1
#
#
# mymdbm=MDBM_Problem(foo,Axis.([δv,bv]))
# MDBM.solve!(mymdbm,3)
#
# x_eval,y_eval=getevaluatedpoints(mymdbm)
# x_sol,y_sol=getinterpolatedsolution(mymdbm)
#
# fig = figure(11);clf()
# PyPlot.scatter(x_eval,y_eval,s=5)
# PyPlot.scatter(x_sol,y_sol)
#
# println("cuccma")
