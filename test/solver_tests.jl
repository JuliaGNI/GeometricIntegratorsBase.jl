using GeometricIntegratorsBase
using Test

import GeometricIntegratorsBase: initsolver

using ..HarmonicOscillator


struct TestMethod <: GeometricMethod end

prob = odeproblem()

caches = CacheDict(prob, TestMethod())


@test initsolver(NoSolver(), TestMethod(), caches) == NoSolver()
