using GeometricIntegratorsBase
using Test

import GeometricIntegratorsBase: default_linesearch, default_options, initsolver

using SimpleSolvers: Backtracking

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

# the expansion phase of SimpleSolvers 0.11 is deliberately left off, since every method here
# solves with `Newton()`, whose direction is already scaled like a Newton step
@test default_linesearch().expand == false
@test default_linesearch(Float32).expand == false


# a solve whose residual goes nowhere is bounded per time step rather than spending
# `max_iterations` on every one of them
@test default_options(TestMethod(), prob).f_stall_window == 50
@test default_options(TestMethod(), prob).min_iterations == 1
