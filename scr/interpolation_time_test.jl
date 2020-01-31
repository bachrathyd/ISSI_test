using PProf

using NumericalIntegration

using LinearAlgebra
using BenchmarkTools
using Statistics
using PyPlot
pygui(true);


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
