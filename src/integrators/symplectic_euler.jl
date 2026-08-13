@doc raw"""
Abstract supertype of the symplectic Euler methods [`SymplecticEulerA`](@ref) and
[`SymplecticEulerB`](@ref).

Both methods are implemented for a *separable* Hamiltonian, that is
```math
H (t, q, p) = T (t, p) + V (t, q) ,
```
so that the vector fields satisfy ``v = v(t,p)`` and ``f = f(t,q)``. Under this assumption the
otherwise implicit partitioned scheme decouples into two explicit substeps and no nonlinear
solver is required. Separability cannot be checked at runtime, so applying these methods to a
non-separable Hamiltonian silently computes something that is neither symplectic nor
consistent with the symplectic Euler method.
"""
abstract type SymplecticEulerMethod <: PODEMethod end

isexplicit(method::SymplecticEulerMethod) = true
isimplicit(method::SymplecticEulerMethod) = false
issymmetric(method::SymplecticEulerMethod) = false
issymplectic(method::SymplecticEulerMethod) = true


@doc raw"""
Symplectic Euler-A Method for separable Hamiltonians.

The momentum is updated explicitly and the position is updated with the new momentum,
```math
\begin{aligned}
p_{n+1} &= p_{n} + h \, f (t_{n}, q_{n}) , \\
q_{n+1} &= q_{n} + h \, v (t_{n+1}, p_{n+1}) .
\end{aligned}
```
This is the adjoint of [`SymplecticEulerB`](@ref). See [`SymplecticEulerMethod`](@ref) for the
separability assumption.
"""
struct SymplecticEulerA <: SymplecticEulerMethod end

@doc raw"""
Symplectic Euler-B Method for separable Hamiltonians.

The position is updated explicitly and the momentum is updated with the new position,
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \, v (t_{n}, p_{n}) , \\
p_{n+1} &= p_{n} + h \, f (t_{n+1}, q_{n+1}) .
\end{aligned}
```
This is the adjoint of [`SymplecticEulerA`](@ref). See [`SymplecticEulerMethod`](@ref) for the
separability assumption.
"""
struct SymplecticEulerB <: SymplecticEulerMethod end


@doc raw"""
Symplectic Euler integrator cache.
"""
struct SymplecticEulerCache{DT} <: PODEIntegratorCache{DT}
    v::Vector{DT}
    f::Vector{DT}

    function SymplecticEulerCache{DT}(ics) where {DT}
        v = zeros(DT, axes(ics.q))
        f = zeros(DT, axes(ics.p))
        new(v, f)
    end
end

function Cache{ST}(problem::AbstractProblemPODE, method::SymplecticEulerMethod; kwargs...) where {ST}
    SymplecticEulerCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblemPODE, ::SymplecticEulerMethod) = SymplecticEulerCache{ST}


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:SymplecticEulerA,<:AbstractProblemPODE})
    # on entry, sol.t holds tₙ₊₁ while sol.q and sol.p still hold qₙ and pₙ
    t̄ = sol.t - timestep(int)

    # compute f = f(t̄, q̄) and update p (separable: f does not depend on p)
    equations(int).f(cache(int).f, t̄, sol.q, sol.p, params)
    sol.p .+= timestep(int) .* cache(int).f

    # compute v = v(t, p) at the new momentum and update q (separable: v does not depend on q)
    equations(int).v(cache(int).v, sol.t, sol.q, sol.p, params)
    sol.q .+= timestep(int) .* cache(int).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:SymplecticEulerB,<:AbstractProblemPODE})
    # on entry, sol.t holds tₙ₊₁ while sol.q and sol.p still hold qₙ and pₙ
    t̄ = sol.t - timestep(int)

    # compute v = v(t̄, p̄) and update q (separable: v does not depend on q)
    equations(int).v(cache(int).v, t̄, sol.q, sol.p, params)
    sol.q .+= timestep(int) .* cache(int).v

    # compute f = f(t, q) at the new position and update p (separable: f does not depend on p)
    equations(int).f(cache(int).f, sol.t, sol.q, sol.p, params)
    sol.p .+= timestep(int) .* cache(int).f
end
