@doc raw"""
Explicit Euler Method.

```math
q_{n+1} = q_{n} + h \, v (t_{n}, q_{n})
```
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

function Cache{ST}(problem::ProblemODE, method::ExplicitEuler; kwargs...) where {ST}
    ExplicitEulerCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemODE, ::ExplicitEuler) = ExplicitEulerCache{ST}

function update!(sol, params, _, int::GeometricIntegrator{<:ExplicitEuler})
    # compute final update
    sol.q .+= timestep(int) .* cache(int).v
end

function integrate_step!(sol, history, params, int::GeometricIntegrator{
        <:ExplicitEuler, <:ProblemODE})
    # on entry, sol.t holds tₙ₊₁ while sol.q still holds qₙ
    t̄ = sol.t - timestep(int)

    # compute vector field v = v(t̄, q̄)
    equations(int).v(cache(int).v, t̄, sol.q, params)

    # compute final update
    update!(sol, params, nothing, int)
end
