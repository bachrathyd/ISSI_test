
using DifferentialEquations

function parameterized_lorenz(du,u,p,t)
  x,y,z = u
  σ,ρ,β = p
  du[1] = dx = σ*(y-x)
  du[2] = dy = x*(ρ-z) - y
  du[3] = dz = x*y - β*z
end

using ParameterizedFunctions
g = @ode_def begin
  dx = σ*(y-x)
  dy = x*(ρ-z) - y
  dz = x*y - β*z
end σ ρ β

# using LinearAlgebra
# f(du,u,p,t) = mul!(du,A,u)

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(g,u0,tspan,p)


using StaticArrays, DifferentialEquations
A  = @SMatrix [ 1.0  0.0 0.0 -5.0
                4.0 -2.0 4.0 -3.0
               -4.0  0.0 0.0  1.0
                5.0 -2.0 2.0  3.0]
u0 = @SMatrix rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
using Plots;
plot(sol)


# const κ = 1; const tau = 2 * pi
const out = zeros(2) # Define a cache variable
function delay_mathieu_model(du,u,h,p,t)
  δ, ϵ, b, κ = p
  h(out, p, t-tau)
  du[1] = - κ *u[1] - (δ+ϵ*cos(t)) *u[2] -b * out[2]
  du[2] = u[1]
end

p = [1.0,2.0,1,0.1]
tau=2*pi
T=2*pi
lags = [tau]
h(out, p, t) = (out.=1.0)
# h(out, p, t, deriv::Type{Val{i}})
tspan = (0.0,T)
u0 = [1.0,1.0]
prob = DDEProblem(delay_mathieu_model,h,tspan,p; constant_lags=lags)


alg = MethodOfSteps(Tsit5())
# alg = MethodOfSteps(BS3())
# alg = MethodOfSteps(Vern6())
sol = solve(prob,alg)#,reltol=1e-8,abstol=1e-8)

using Plots
plot(sol,vars=(0,2))

# plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
#      xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
# plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
sol[2][3]
plot(sol,vars=(0,2))
