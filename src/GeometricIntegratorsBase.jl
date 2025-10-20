module GeometricIntegratorsBase

using GeometricBase
using GeometricEquations
using GeometricSolutions
using LinearAlgebra
using OffsetArrays
using SimpleSolvers

import Base: Callable
using Unicode: normalize

import GeometricBase: equations, initialguess, nconstraints, parameters, tableau
import GeometricBase: timestep, timespan
import GeometricBase: reset!, solutionstep!, update!
import GeometricBase: integrate, integrate!
import GeometricBase: periodic, verifyrange
import GeometricBase: AbstractVariable, AbstractScalarVariable, AbstractStateVariable


# compat workaround
Base.ndims(prob::EquationProblem) = length(vec(prob.ics.q))


export update!, reset!


export Extrapolation
export EulerExtrapolation,
    MidpointExtrapolation,
    HermiteExtrapolation
export extrapolate!, solutionstep!

include("extrapolation/extrapolation.jl")
include("extrapolation/aitken_neville.jl")
include("extrapolation/euler.jl")
include("extrapolation/hermite.jl")
include("extrapolation/midpoint.jl")


export InitialGuess, NoInitialGuess
export initialguess!

include("initialguess.jl")


export GeometricMethod
export default_solver, default_iguess, default_projection
export initmethod, internal_variables, nullvector, tableau
export isexplicit, isimplicit, issymmetric, issymplectic, isenergypreserving, isstifflyaccurate
export isodemethod, ispodemethod, ishodemethod, isiodemethod, islodemethod, issodemethod
export isdaemethod, ispdaemethod, ishdaemethod, isidaemethod, isldaemethod

include("method.jl")


export SolutionStep
export parameters, solution, vectorfield
export current, previous, history

include("solutionstep.jl")
include("solution.jl")


export Cache, CacheDict, NoCache
export cache, nlsolution

include("cache.jl")


export NoSolver

include("solvers.jl")


export GeometricIntegrator
export integrate, integrate!, integrate_step!

include("integrate.jl")
include("integrator.jl")

export NoProjection, projection

include("projection.jl")

export ExplicitEuler, ImplicitEuler

include("euler/explicit_euler.jl")
include("euler/implicit_euler.jl")

end
