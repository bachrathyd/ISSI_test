using NumericalIntegration

using LinearAlgebra
using BenchmarkTools
using Statistics
using PyPlot
pygui(true);


function mathieu(z, zτ, t, κ, δ, ϵ, b, τ, ω)
        dzdt=zeros(ComplexF64,size(z))
                dzdt[1]=z[2]
                dzdt[2]=-(δ+(ϵ * cos(ω * t)))*z[1] - κ * z[2] + b * zτ[1]
                return dzdt
end

function mathieu!(dzdt, z, zτ, t, κ, δ, ϵ, b, τ, ω)
        dzdt=zeros(ComplexF64,size(z))
                dzdt[1]=z[2]
                dzdt[2]=-(δ+(ϵ * cos(ω * t)))*z[1] - κ * z[2] + b * zτ[1]
                # return dzdt
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
function Sinterp(Sloc,indf::Float64)
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
function Sinterp(S,i::Integer)
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

# n=2^5
# m=Int64(n*(1 / ω))
# RKtype="RK4"
function mathiestab( κ, δ, ϵ, b, τ, ω, n,m,mult,iter,RKtype::String)::Array{ComplexF64,1}

dt=τ/n
S=rand(ComplexF64,2,n+1+m,mult)
λj=Array{ComplexF64}(undef,mult)
#
# k=Array{ComplexF64}(undef,2)
# k1=Array{ComplexF64}(undef,2)
# k2=Array{ComplexF64}(undef,2)
# k3=Array{ComplexF64}(undef,2)
# k4=Array{ComplexF64}(undef,2)

for kinter in 1:iter
        #Threads.@threads
        for mi in 1:mult
                for kt in n+1:n+m
                if RKtype=="EE"
                        @views k1=mathieu(S[:,kt,mi], Sinterp(S[:,:,mi],kt-n),(kt-(n+1))*dt, κ, δ, ϵ, b, τ, ω)
                        k=k1 ./ 1.0
                # elseif RKtype=="RK2-fast"
                #         @views k1=mathieu(S[:,kt,mi]      ,S[:,kt-n  ,mi],(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                #         @views k2=mathieu(S[:,kt,mi]+k1*dt,S[:,kt-n+1,mi],(kt-(n+1)+1)*dt, κ, δ, ϵ, b, τ, ω)
                #          k=(k1+k2) ./ 2.0
                elseif RKtype=="RK2-Heun"
                        @views k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        @views k2=mathieu(S[:,kt,mi]+k1*dt,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+k2) ./ 2.0
                elseif RKtype=="RK2-Ralston"
                        @views k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        @views k2=mathieu(S[:,kt,mi]+k1*dt*2.0/3.0,Sinterp(S[:,:,mi],kt-n+2.0/3.0),(kt-(n+1)+2.0/3.0)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+ 3* k2) ./ 4.0
                elseif RKtype=="RK3"
                        @views k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        @views k2=mathieu(S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                        @views k3=mathieu(S[:,kt,mi]-k1*dt*1.0+k2*dt*2.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+4*k2+k3) ./ 6.0
                elseif RKtype=="RK3-strong"
                        @views k1=mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)  ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                        @views k2=mathieu(S[:,kt,mi]+k1*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                        @views k3=mathieu(S[:,kt,mi]+k1*dt*0.25+k2*dt*0.25,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                        k=(k1+k2+4*k3) ./ 6.0
                elseif RKtype=="RK4"
                         @views k1 = mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                         @views k2 = mathieu(S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                         @views k3 = mathieu(S[:,kt,mi]+k2*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                         @views k4 = mathieu(S[:,kt,mi]+k3*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                          # mathieu!(k1,S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                          # mathieu!(k2,S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                          # mathieu!(k3,S[:,kt,mi]+k2*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                          # mathieu!(k4,S[:,kt,mi]+k3*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                         k = (k1+2*k2+2*k3+k4) ./ 6.0
                 elseif RKtype=="RK5"
                          @views k1 = mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                          @views k2 = mathieu(S[:,kt,mi]+k1*dt*0.25,Sinterp(S[:,:,mi],kt-n+0.25),(kt-(n+1)+0.25)*dt, κ, δ, ϵ, b, τ, ω)
                          @views k3 = mathieu(S[:,kt,mi]+k1*dt*0.125+k2*dt*0.125,Sinterp(S[:,:,mi],kt-n+0.25),(kt-(n+1)+0.25)*dt, κ, δ, ϵ, b, τ, ω)
                          @views k4 = mathieu(S[:,kt,mi]-k2*dt*0.5+k3*dt,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                          @views k5 = mathieu(S[:,kt,mi]+k1*dt*3.0/16.0+k4*dt*9.0/16.0,Sinterp(S[:,:,mi],kt-n+0.75),(kt-(n+1)+0.75)*dt, κ, δ, ϵ, b, τ, ω)
                          @views k6 = mathieu(S[:,kt,mi]-k1*dt*3.0/7.0+k2*dt*2.0/7.0+k3*dt*12.0/7.0-k4*dt*12.0/7.0+k5*dt*8.0/7.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                           # mathieu!(k1,S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                           # mathieu!(k2,S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                           # mathieu!(k3,S[:,kt,mi]+k2*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                           # mathieu!(k4,S[:,kt,mi]+k3*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                          k = (7*k1    +32*k3+12*k4+32*k5+7*k6) ./ 90.0
                  elseif RKtype=="RK5-2"
                           @views k1 = mathieu(S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                           @views k2 = mathieu(S[:,kt,mi]+k1*dt*1/5,Sinterp(S[:,:,mi],kt-n+1/5),(kt-(n+1)+1/5)*dt, κ, δ, ϵ, b, τ, ω)
                           @views k3 = mathieu(S[:,kt,mi]+k2*dt*2/5,Sinterp(S[:,:,mi],kt-n+2/5),(kt-(n+1)+2/5)*dt, κ, δ, ϵ, b, τ, ω)
                           @views k4 = mathieu(S[:,kt,mi]+k1*dt*9/4-k2*dt*5+k3*dt*15/4,Sinterp(S[:,:,mi],kt-n+1),(kt-(n+1)+1)*dt, κ, δ, ϵ, b, τ, ω)
                           @views k5 = mathieu(S[:,kt,mi]-k1*dt*63/100+k2*dt*9/5-k3*dt*13/20+k4*dt*2/25,Sinterp(S[:,:,mi],kt-n+3/5),(kt-(n+1)+3/5)*dt, κ, δ, ϵ, b, τ, ω)
                           @views k6 = mathieu(S[:,kt,mi]-k1*dt*6/25+k2*dt*4/5+k3*dt*2/15+k4*dt*8/75,Sinterp(S[:,:,mi],kt-n+4/5),(kt-(n+1)+4/5)*dt, κ, δ, ϵ, b, τ, ω)
                            # mathieu!(k1,S[:,kt,mi]      ,Sinterp(S[:,:,mi],kt-n)        ,(kt-(n+1)  )*dt, κ, δ, ϵ, b, τ, ω)
                            # mathieu!(k2,S[:,kt,mi]+k1*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                            # mathieu!(k3,S[:,kt,mi]+k2*dt*0.5,Sinterp(S[:,:,mi],kt-n+0.5),(kt-(n+1)+0.5)*dt, κ, δ, ϵ, b, τ, ω)
                            # mathieu!(k4,S[:,kt,mi]+k3*dt*1.0,Sinterp(S[:,:,mi],kt-n+1.0),(kt-(n+1)+1.0)*dt, κ, δ, ϵ, b, τ, ω)
                           k = (17/144*k1+25/36*k3+1/72*k4-25/72*k5+25/48*k6)
                else
                        println("RKtype")
                        error("nincs ilyen integrator")
                end
                 S[:,kt+1,mi] = S[:,kt,mi]+dt*k;
                end
        end



                Sj=reshape(S[:,1:n+1,:],(2*(n+1),mult))
                Vj=reshape(S[:,end-n:end,:],(2*(n+1),mult))
                        # Sj=S[:,1:n+1,1]
                        # Vj=S[:,end-n:end,1]

        # SS=[integrate(1:size(Sj,1), Sj[:,m1] .* Sj[:,m2], SimpsonEvenFast()) for m1 in 1:mult, m2 in 1:mult]
        # SV=[integrate(1:size(Sj,1), Sj[:,m1] .* Vj[:,m2], SimpsonEvenFast()) for m1 in 1:mult, m2 in 1:mult]

        # SS=[integrate(1:size(Sj,1), Sj[:,m1] .* Sj[:,m2], TrapezoidalEvenFast()) for m1 in 1:mult, m2 in 1:mult]
        # SV=[integrate(1:size(Sj,1), Sj[:,m1] .* Vj[:,m2], TrapezoidalEvenFast()) for m1 in 1:mult, m2 in 1:mult]
                # println(n)
                # println(size(Sj,1))
                # SS=[integrate(1:size(Sj,1), Sj[:,m1] .* Sj[:,m2], RombergEven()) for m1 in 1:mult, m2 in 1:mult]
                # SV=[integrate(1:size(Sj,1), Sj[:,m1] .* Vj[:,m2], RombergEven()) for m1 in 1:mult, m2 in 1:mult]

                        SS=(transpose(Sj)*Sj)
                        SV=(transpose(Sj)*Vj)



        #Hj=(Sj\Vj)
        Hj=(SS\SV)
        # Hj=((transpose(Sj)*Sj)\(transpose(Sj)*Vj))
        (λj,Gj)=eigen(Hj)
                # end


        S[:,1:n+1,:] = (Vj*(Gj/diagm(λj)))[:]
        # S[:,1:n+1,1] = (Vj*(Gj/diagm(λj)))[:]
        # S[:,1:n+1,:] = (Vj*(Gj))[:]


        # figure(3)
        # clf()
        # plot((1:n+1+m).*dt,S[2,1:n+1+m,:])
        # plot((1:n+1).*dt,S[2,1:n+1,:])
        # pause(0.01)
        #
        # figure(4)
        # scatter(kinter .* ones(mult) .+ 0.1,abs.(λj))
        # pause(0.01)
        #
        # figure(8)
        # clf()
        # surf(abs.(Gj))
end

return λj
end


println("----------=======<<<<<<<<<<<<<#######   START   #######>>>>>>>>>>>>>>========----------")
 const κ=0.05
 const δ=5
 const ϵ=1.0
 const b=0.1
 const τ=2pi
 const ω=0.25

 mult=8 # number of computed eigenvalues
 iter=8 # number of iteration in issi




@time mathiestab( κ, δ, ϵ, b, τ, ω, 2^5, Int64(2^5*(1 / ω)),mult,iter,"RK4")
 #
 # for NN in 5:12
 # @time mathiestab( κ, δ, ϵ, b, τ, ω, 2^NN, Int64(2^NN*(1 / ω)),mult,iter,"RK4")
 # end
 #
 # @btime mathiestab( κ, δ, ϵ, b, τ, ω, 2^6, Int64(2^6*(1 / ω)),mult,iter,"RK4") #First run "time" is too long



# figure(4)
# clf()
@elapsed mathiestab( κ, δ, ϵ, b, τ, ω, 2^5, Int64(2^5*(1 / ω)),mult,iter,"RK2-Ralston")
@debug mathiestab( κ, δ, ϵ, b, τ, ω, 2^5, Int64(2^5*(1 / ω)),mult,iter,"RK2-Ralston") #First run "time" is too long
# @elapsed mu_ref=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω,10^6,m,mult,iter,"RK4")))
# @code_warntype mathiestab( κ, δ, ϵ, b, τ, ω, 2^4, Int64(2^4*(1 / ω)),mult,iter,"RK4") #First run "time" is too long



println("Computing the reference")
# n=2^19 # 16 sec
# n=2^19#16
n=2^15#16
m=Int64(n*(1 / ω))
@time mu_ref=mathiestab( κ, δ, ϵ, b, τ, ω, n,m,mult,iter,"RK4")
if imag(mu_ref[end])<0
        conj!(mu_ref)
end

println("--------------------===================------------------------------")
integrators1=["EE"]
integrators2=["RK2-Ralston","RK2-Heun"]
integrators3=["RK3","RK3-strong"]
integrators4=["RK4"]
integrators5=["RK5","RK5-2"]
integrators=[integrators1...,integrators2...,integrators3...,integrators4...,integrators5...]


# ni=2 .^ (4:19)
ni=2 .^ (2:13)
mumax=zeros(ComplexF64,length(ni),length(integrators),mult)
timelength=zeros(length(ni),length(integrators))
timestd=zeros(length(ni),length(integrators))

# figure(1)
# clf()
# figure(2)
# clf()
# figure(3)
# clf()

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

Nsample=1

for kn in 1:length(ni)
        n=ni[kn]
        m=Int64(n*(1 / ω))

        for ksample in 1:Nsample#maximum([Nsample-Int(log(ni[kn])/log(2)),1])
                println("<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>")
                println([ksample,kn,ni[kn]])
                for kRKi in 1:length(integrators)
                        integratortype=integrators[kRKi]
                        print(integratortype)
                        print("  -  ")

        # timeing = @benchmark (mumax[kn]=maximum(abs.(mathiestab( κ, δ, ϵ, b, τ, ω,n,m,mult,iter,"RK4"))))
        # timelength[kn]=mean(timeing.times)
        # timestd[kn]=std(timeing.times)
        timeing = @elapsed (mui_set=mathiestab( κ, δ, ϵ, b, τ, ω,n,m,mult,iter,integratortype))
        if imag(mui_set[end])<0
                conj!(mui_set)
        end
        mumax[kn,kRKi,:] .= mui_set
        timelength[kn,kRKi]+=timeing
        println(timeing)


        figure(1)
        scatter(n,abs(mumax[kn,kRKi,end] - mu_ref[end]))
        # scatter(n,abs(mumax[kn,kRKi,end]) - abs(mu_ref[end]))
        figure(2)
        scatter(n,timeing)
        figure(3)
        scatter(timeing,abs(mumax[kn,kRKi,end] - mu_ref[end]))
        # scatter(timeing,abs(mumax[kn,kRKi,end]) - abs(mu_ref[end]))

        # figure(5)
        # scatter(abs(mumax[kn,kRKi,end] - mu_ref[end]),abs(mumax[kn,kRKi,end]) - abs(mu_ref[end]))
                end
        end


        figure(1)
        plot(ni[1:kn,:],abs.(mumax[1:kn,:,end] .- mu_ref[end]))
        savefig("1_v_nt")
        # savefig("1_t")

        figure(2)
        plot(ni[1:kn,:],timelength[1:kn,:]./Nsample)
        savefig("2_v_nt")
        # savefig("2_t")

        figure(3)
        plot(timelength[1:kn,:]./Nsample,abs.(mumax[1:kn,:,end] .- mu_ref[end]))
        savefig("3_v_nt")
        # savefig("3_t")

        # figure(5)
        # # clf()
        # ax=gca()
        # ax.set_xscale("log")
        # ax.set_yscale("log")
        # grid("on")
        # savefig("5_nt")
end


figure(1)
plot(ni,abs.(mumax[:,:,end] .- mu_ref[end]))
savefig("1_v_nt")
# savefig("1_t")

figure(2)
plot(ni,timelength./Nsample)
savefig("2_v_nt")
# savefig("2_t")

figure(3)
plot(timelength./Nsample,abs.(mumax[:,:,end] .- mu_ref[end]))
savefig("3_v_nt")
# savefig("3_t")

# figure(5)
# savefig("5_nt")
println("---------   <<<<<<<<  END >>>>>>>>>   ----------")
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




figure(11)
 clf()
ax=gca()
ax.set_xscale("log")
ax.set_yscale("log")
grid("on")
title("n / error")
xlabel("n")
ylabel("error")

A2=hcat([-1,1,-1,1],[1,1,1,1],[3,-2,1,0],[3,2,1,0])'
A1=hcat([-1,1],[1,1])'

w2=[0,2/3,0,2]'*inv(A2)
w1=[0,2]'*inv(A1)


# t=0:0.1:pi/2
# y=sin.(try)
# dydy=cos.(t)
for kk in 3:20
        dt=pi/(2.0 ^ (kk))

        t=0:dt:pi/2
        # y=t
        # dydt=1 .+ 0 .* (t)
        y=sin.(t)
        dydt=cos.(t)

        # # y=t
        # # dydt=1 .+ 0 .* (t)
        # y=t.^7.0
        # dydt=0.0 .+ 7.0 .* t .^6.0
        I=0.0
        for l in 1:(length(t)-1)
                # global I
                # I += dt*y[l]
                # I += dt*y[l+1]
                 # I += dt/2*(y[l]+y[l+1])
                  # I += (dt/2)*w1*[y[l],y[l+1]]
                I += (dt/2)*w2*[y[l],y[l+1],dydt[l]*(dt/2),dydt[l+1]*(dt/2)]
        end

        figure(11)
        scatter(length(t),abs(I-1.0))
        println(I)

        IS=integrate(t, y, SimpsonEvenFast())
        IT=integrate(t, y, TrapezoidalEvenFast())
        IR=integrate(t, y, RombergEven())
        scatter(length(t),abs(sum(y)*dt-1.0), color=:green)
        scatter(length(t),abs(IS-1.0), color=:red)
        scatter(length(t),abs(IT-1.0), color=:blue)
        scatter(length(t),abs(IR-1.0), color=:black)
end

length(dt)
