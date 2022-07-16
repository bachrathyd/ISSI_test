module SimulationBasedStability

using Reexport
@reexport using LinearAlgebra
#@reexport using SparseArrays
#@reexport using StaticArrays
#@reexport using Arpack
@reexport using DifferentialEquations
@reexport using Random
#using QuadGK
#using Lazy: iterated, take

include("structures.jl")

#include("functions.jl")
#include("functions_discretization.jl")
#include("functions_method.jl")

export 
dynamic_problem,
histopryremapp,
CompRand,
compute_eig!
#SemiDiscretization, NumericSD, 
#ProportionalMX,
#Delay,DelayMX,
#Additive,
#LDDEProblem,
#DiscreteMapping, DiscreteMapping_1step,
#fixPointOfMapping, spectralRadiusOfMapping,
#DiscreteMapping_LR,spectralRadiusOfMappingLR

end # module
