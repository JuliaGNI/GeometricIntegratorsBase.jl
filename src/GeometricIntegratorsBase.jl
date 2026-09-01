module GeometricIntegratorsBase

using GeometricBase
using GeometricEquations
using GeometricSolutions
using LinearAlgebra
using OffsetArrays
using SimpleSolvers

using Unicode: normalize

import Base: Callable
import GeometricBase: equations, nconstraints, parameters, tableau
import GeometricBase: initialguess, initialstate, initialtime
import GeometricBase: timestep, timespan
import GeometricBase: reset!, solutionstep!, update!
import GeometricBase: solution, solutionkeys, state, vectorfield
import GeometricBase: integrate, integrate!
import GeometricBase: periodic, verifyrange
import GeometricBase: AbstractVariable, AbstractScalarVariable, AbstractStateVariable,
                      TimeVariable
import GeometricBase: NoSolver
import SimpleSolvers: NonlinearSolverException, NonlinearSolverMethod

# The problem unions of GeometricEquations also cover the constrained variants, that is
# `AbstractProblemODE` includes `DAEProblem`, `AbstractProblemPODE` includes `PDAEProblem` and
# `HDAEProblem`, and `AbstractProblemIODE` includes `IDAEProblem` and `LDAEProblem`. An
# integrator that does not enforce a constraint must not accept those, as it would otherwise
# silently return a solution that ignores it, so the integrators in this package dispatch on
# the unconstrained unions below. A combination that is left unimplemented is rejected with an
# `ArgumentError` naming method and problem, by `initsolver` for a method that solves and by
# the `integrate_step!` fallback in `src/integrator.jl` for one that does not.
const ProblemODE{DT, TT} = Union{ODEProblem{DT, TT}, SubstepProblem{DT, TT}}
const ProblemIODE{DT, TT} = Union{IODEProblem{DT, TT}, LODEProblem{DT, TT}}
const ProblemPODE{DT, TT} = Union{PODEProblem{DT, TT}, HODEProblem{DT, TT}}

export update!, reset!
export State
export NoSolver

export InitialGuess, NoInitialGuess
export initialguess!

include("initialguess.jl")

export Extrapolation, NoExtrapolation
export EulerExtrapolation,
       MidpointExtrapolation,
       HermiteExtrapolation,
       NormalizedHermiteExtrapolation
export extrapolate!, solutionstep!
export default_extrapolation

include("extrapolation/extrapolation.jl")
include("extrapolation/aitken_neville.jl")
include("extrapolation/euler.jl")
include("extrapolation/hermite.jl")
include("extrapolation/midpoint.jl")

export GeometricMethod
export default_solver, default_iguess, default_projection
export initmethod, internal_variables, nullvector, tableau
export isexplicit, isimplicit, issymmetric, issymplectic, isenergypreserving,
       isstifflyaccurate
export isodemethod, ispodemethod, ishodemethod, isiodemethod, islodemethod, issodemethod
export isdaemethod, ispdaemethod, ishdaemethod, isidaemethod, isldaemethod

include("method.jl")

export SolutionStep
export parameters, solution, state, vectorfield
export current, history, previous

include("solutionstep.jl")
include("solution.jl")

export Cache, CacheDict, NoCache
export cache, nlsolution

include("cache.jl")

export check_solver_status

include("solvers.jl")

export GeometricIntegrator
export integrate, integrate!, integrate_step!

include("integrate.jl")
include("integrator.jl")

export NoProjection, projection

include("projection.jl")

export ExplicitEuler, ImplicitEuler
export SymplecticEulerMethod, SymplecticEulerA, SymplecticEulerB
export ImplicitMidpoint, CrankNicolson

include("integrators/explicit_euler.jl")
include("integrators/implicit_euler.jl")
include("integrators/symplectic_euler.jl")
include("integrators/implicit_midpoint.jl")
include("integrators/crank_nicolson.jl")

end
