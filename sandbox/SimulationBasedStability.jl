module SimulationBasedStability

using Reexport
@reexport using LinearAlgebra
@reexport using DifferentialEquations
@reexport using StaticArrays
@reexport using Random
@reexport using Statistics
@reexport using Interpolations
@reexport using Arpack
#using QuadGK
#using Lazy: iterated, take

include("structures.jl")

export SemiDiscretization,
dynamic_problem,
iterate!,
compute_eig!,
spectralRadiusOfMapping,
histopryremap,
getvalues,
SVi1real,
CompRand,
ComplexinterpFun,
LinearHistorymapping



end # module
