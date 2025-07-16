
abstract type Cache{DT} end

struct NoCache{DT} <: Cache{DT} end

Cache{DT}(::AbstractProblem, ::GeometricMethod) where {DT} = NoCache{DT}()
Cache(problem::AbstractProblem, method::GeometricMethod) = Cache{datatype(problem)}(problem, method)

CacheType(T, problem::AbstractProblem, method::GeometricMethod) = NoCache{T}

reset!(::Cache) = nothing


abstract type IntegratorCache{D,DT} <: Cache{DT} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DELEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end

# cache(::AbstractIntegrator, DT) = nothing
# cache(::AbstractIntegrator) = nothing

copy_internal_variables!(::SolutionStep, ::Union{IntegratorCache,NoCache}) = nothing

nlsolution(::IntegratorCache) = missing


struct CacheDict{PT,MT}
    problem::PT
    method::MT
    caches::Dict{UInt64,Cache}

    function CacheDict(prob::AbstractProblem, method::GeometricMethod)
        new{typeof(prob),typeof(method)}(prob, method, Dict{UInt64,Cache}())
    end
end

@inline function Base.getindex(c::CacheDict, ST::DataType)
    key = hash(Threads.threadid(), hash(ST))
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = Cache{ST}(c.problem, c.method)
    end::CacheType(ST, c.problem, c.method)
end

function updateall!(c::CacheDict, data)
    for (key, cache) in c.caches
        update!(c.caches[key], data...)
    end
end

cache(c::CacheDict) = c[datatype(c.problem)]
nlsolution(c::CacheDict) = nlsolution(cache(c))
