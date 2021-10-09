module CombinedUncertainDiffEq

# Imported functions
import ReachabilityAnalysis
const RA = ReachabilityAnalysis
const ReachInterval = RA.Interval

using ComponentArrays: ComponentArrays, ComponentArray, getaxes, getdata
using ConcreteStructs: @concrete
using DifferentialEquations: DifferentialEquations, ODEProblem, EnsembleProblem, solve, remake, isinplace
using DiffEqUncertainty: DiffEqUncertainty, Koopman, MonteCarlo, expectation, _rand
using Distributions: Distribution, Uniform
using IntervalArithmetic: Interval, ±, (..), inf, sup
using ReachabilityAnalysis: HalfSpace, HPolyhedron, Hyperrectangle,
        ×, cartesian_product, concretize
using Setfield: @set!
using Statistics: mean
using Symbolics: @variables, Num
using UnPack: @unpack


# Exported functions
export ±, (..), Interval, Uniform

include("soft_functions.jl")

include("sim_utils.jl")
export Saturated, SoftSaturated, WithControls, reachable_constraint

include("conversions.jl")

include("solve.jl")
export mc_solve, mc_expectation, koopman_expectation, reachability_problem, reachability_solve

end
