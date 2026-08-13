using GeometricIntegratorsBase
using Test

import GeometricIntegratorsBase: DEFAULT_F_STALL_WINDOW
import GeometricIntegratorsBase: default_linesearch, default_options, initsolver

using SimpleSolvers: Backtracking, F_STALL_REPORT_MINIMUM

using ..HarmonicOscillator


struct TestMethod <: GeometricMethod end

prob = odeproblem()

caches = CacheDict(prob, TestMethod())


@test initsolver(NoSolver(), TestMethod(), caches) == NoSolver()


# the untyped `default_linesearch` is the downstream hook and must keep working unchanged,
# while the typed one builds the method at the working precision directly
@test default_linesearch() isa Backtracking{Float64}
@test default_linesearch(TestMethod()) isa Backtracking{Float64}
@test default_linesearch(Float32) isa Backtracking{Float32}
@test default_linesearch(Float32, TestMethod()) isa Backtracking{Float32}

@test default_linesearch(BigFloat) isa Backtracking{BigFloat}

# the typed method must not intercept a call meant for the untyped hook: passing a method *type*
# instead of an instance still resolves to the untyped one, as it did before the typed one existed
@test default_linesearch(TestMethod) isa Backtracking{Float64}
@test which(default_linesearch, Tuple{Type{TestMethod}}) !== which(default_linesearch, Tuple{Type{Float32}})

# ... while the constraint stays wide enough for a real type that is not an `AbstractFloat`: the
# `ForwardDiff.Dual` a caller differentiating through an integration arrives with has to reach the
# typed method, which is what `T<:Real` rather than `T<:AbstractFloat` buys
@test which(default_linesearch, Tuple{Type{Rational{Int}}}) === which(default_linesearch, Tuple{Type{Float32}})

# the expansion phase of SimpleSolvers 0.11 is deliberately left off, since every method here
# solves with `Newton()`, whose direction is already scaled like a Newton step
@test default_linesearch().expand == false
@test default_linesearch(Float32).expand == false


# a solve whose residual goes nowhere is bounded per time step rather than spending
# `max_iterations` on every one of them
@test default_options(TestMethod(), prob).f_stall_window == DEFAULT_F_STALL_WINDOW
@test default_options(TestMethod(), prob).min_iterations == 1

# the window has to stay above the count below which SimpleSolvers holds the no-progress
# proportion to be no evidence of anything: giving up on less than it will report on is incoherent
@test DEFAULT_F_STALL_WINDOW > F_STALL_REPORT_MINIMUM
