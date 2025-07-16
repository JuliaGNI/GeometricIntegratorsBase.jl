using GeometricBase
using GeometricEquations
using GeometricIntegratorsBase
using Test

import GeometricIntegratorsBase: IntegratorCache, CacheType

using ..HarmonicOscillator


prob = odeproblem()


struct TestMethod <: GeometricMethod end

struct TestCache{DT,D} <: IntegratorCache{DT,D}
    x::Vector{DT}
    q::Vector{DT}

    function TestCache{DT}(q) where {DT}
        x = zeros(DT, 2 * length(axes(q)))
        q = zero(similar(q, DT))
        new{DT,length(x)}(x, q)
    end
end

@test Cache{Float64}(prob, TestMethod()) == NoCache{Float64}()
@test CacheType(Float64, prob, TestMethod()) == NoCache{Float64}


function Cache{ST}(problem::AbstractProblem, ::TestMethod; kwargs...) where {ST}
    TestCache{ST}(initial_conditions(problem).q; kwargs...)
end

CacheType(ST, ::AbstractProblem, ::TestMethod) = TestCache{ST}

@test typeof(Cache(prob, TestMethod())) <: TestCache{datatype(prob)}
@test typeof(Cache{datatype(prob)}(prob, TestMethod())) <: TestCache{datatype(prob)}

for ST ∈ (Float32, Float64, Rational{Int64})
    @test CacheType(ST, prob, TestMethod()) == TestCache{ST}
end


tcache = Cache(prob, TestMethod())
caches = CacheDict(prob, TestMethod())

@test cache(caches) == caches[datatype(prob)]
@test ismissing(nlsolution(caches))

@test caches[Float32].x == zeros(Float32, axes(tcache.x))
@test caches[Float32].q == zeros(Float32, axes(tcache.q))
