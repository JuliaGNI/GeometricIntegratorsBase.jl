"""
Explicit Euler Method.

"""
struct ExplicitEuler <: ODEMethod end
# $(reference(Val(:ExplicitEuler)))

isexplicit(method::ExplicitEuler) = true
isimplicit(method::ExplicitEuler) = false
issymmetric(method::ExplicitEuler) = false
issymplectic(method::ExplicitEuler) = false

@doc raw"""
Explicit Euler integrator cache.
"""
struct ExplicitEulerCache{DT} <: ODEIntegratorCache{DT}
    v::Vector{DT}

    function ExplicitEulerCache{DT}(ics) where {DT}
        v = zeros(DT, axes(ics.q))
        new(v)
    end
end

function Cache{ST}(problem::AbstractProblem, method::ExplicitEuler; kwargs...) where {ST}
    ExplicitEulerCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblem, ::ExplicitEuler) = ExplicitEulerCache{ST}


function update!(sol, params, _, int::GeometricIntegrator{<:ExplicitEuler})
    # compute final update
    sol.q .+= timestep(int) .* cache(int).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ExplicitEuler,<:AbstractProblemODE})
    # compute vector field
    equations(int).v(cache(int).v, sol.t, sol.q, params)

    # compute final update
    update!(sol, params, nothing, int)
end
