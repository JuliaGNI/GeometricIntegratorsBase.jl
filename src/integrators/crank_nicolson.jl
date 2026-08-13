@doc raw"""
Crank-Nicolson Method, also known as the trapezoidal rule.

```math
q_{n+1} = q_{n} + \frac{h}{2} \, \big[ v (t_{n}, q_{n}) + v (t_{n+1}, q_{n+1}) \big]
```
The nonlinear solver solves for the vector field at the new time step,
``V = v(t_{n+1}, q_{n+1})``, while ``\bar{V} = v(t_{n}, q_{n})`` is computed once per time step.
The method is symmetric and second order, but not symplectic (it is conjugate to a symplectic
method).
"""
struct CrankNicolson <: ODEMethod end

isexplicit(method::CrankNicolson) = false
isimplicit(method::CrankNicolson) = true
issymmetric(method::CrankNicolson) = true
issymplectic(method::CrankNicolson) = false


@doc raw"""
Crank-Nicolson integrator cache.

The field `v̄` holds the vector field at the previous time step, ``v(t_{n}, q_{n})``, which is
constant throughout the nonlinear solve.
"""
struct CrankNicolsonCache{DT} <: ODEIntegratorCache{DT}
    x::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}
    v̄::Vector{DT}

    function CrankNicolsonCache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)))
        q = zeros(DT, axes(ics.q))
        v = zeros(DT, axes(ics.q))
        v̄ = zeros(DT, axes(ics.q))
        new(x, q, v, v̄)
    end
end

nlsolution(cache::CrankNicolsonCache) = cache.x

function Cache{ST}(problem::AbstractProblem, method::CrankNicolson; kwargs...) where {ST}
    CrankNicolsonCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblem, ::CrankNicolson) = CrankNicolsonCache{ST}


solversize(::CrankNicolson, problem::AbstractProblemODE) = length(vec(initial_conditions(problem).q))

default_solver(::CrankNicolson) = Newton()
default_iguess(::CrankNicolson) = HermiteExtrapolation()

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson})
    # temporary solution
    ig = (t=sol.t, q=cache(int).q, q̇=cache(int).v)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v

    # v̄ = v(t̄, q̄) is constant during the solve, so it is read from the cache at working
    # precision, which also holds it during the automatic differentiation of the residual
    v̄ = cache(int).v̄

    # compute q = q̄ + Δt/2 * (v̄ + x) (v = x)
    q .= sol.q .+ (timestep(int) / 2) .* (v̄ .+ x)

    # compute v = v(q)
    equations(int).v(v, sol.t, q, params)
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:CrankNicolson}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v

    # compute b = - (v-x)
    b .= v .- x
end


# Compute stages of Crank-Nicolson methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CrankNicolson}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= (timestep(int) / 2) .* (cache(int).v̄ .+ cache(int, DT).v)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:AbstractProblemODE})
    # compute vector field at the previous time step
    equations(int).v(cache(int).v̄, sol.t - timestep(int), sol.q, params)

    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
