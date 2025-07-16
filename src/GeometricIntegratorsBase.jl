module GeometricIntegratorsBase

using GeometricBase
using GeometricEquations
using GeometricSolutions
using OffsetArrays
using SimpleSolvers

using Unicode: normalize

import GeometricBase: parameters, tableau
import GeometricBase: timestep, timespan
import GeometricBase: reset!, update!
import GeometricBase: periodic, verifyrange
import GeometricBase: AbstractVariable, AbstractScalarVariable, AbstractStateVariable


export update!, reset!


export Extrapolation

include("extrapolation.jl")


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

end
