module CombinedUncertainDiffEq

# Imported functions
using ComponentArrays: ComponentArray, getaxes, getdata
using DifferentialEquations: ODEProblem, EnsembleProblem, solve
using DiffEqUncertainty: DiffEqUncertainty, Koopman, MonteCarlo, expectation, _rand
using Distributions: Distribution, Uniform
using IntervalArithmetic: Interval, ±, (..)
using Setfield: @set!


# Exported functions
export ±, (..), Interval, Uniform

include("conversions.jl")

include("expectation.jl")
export koopman_expectation, mc_solve

end
