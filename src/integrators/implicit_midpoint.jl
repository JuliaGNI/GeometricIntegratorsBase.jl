@doc raw"""
Implicit Midpoint Method.

```math
q_{n+1} = q_{n} + h \, v \bigg( t_{n} + \frac{h}{2} , \frac{q_{n} + q_{n+1}}{2} \bigg)
```
The nonlinear solver solves for the stage vector field ``V = v(t_{n} + h/2, q_{n} + h V / 2)``,
so that the update reads ``q_{n+1} = q_{n} + h V``.
"""
struct ImplicitMidpoint <: ODEMethod end

isexplicit(method::ImplicitMidpoint) = false
isimplicit(method::ImplicitMidpoint) = true
issymmetric(method::ImplicitMidpoint) = true
issymplectic(method::ImplicitMidpoint) = true


@doc raw"""
Implicit midpoint integrator cache.
"""
struct ImplicitMidpointCache{DT} <: ODEIntegratorCache{DT}
    x::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}

    function ImplicitMidpointCache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)))
        q = zeros(DT, axes(ics.q))
        v = zeros(DT, axes(ics.q))
        new(x, q, v)
    end
end

nlsolution(cache::ImplicitMidpointCache) = cache.x

function Cache{ST}(problem::AbstractProblem, method::ImplicitMidpoint; kwargs...) where {ST}
    ImplicitMidpointCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblem, ::ImplicitMidpoint) = ImplicitMidpointCache{ST}


solversize(::ImplicitMidpoint, problem::AbstractProblemODE) = length(vec(initial_conditions(problem).q))

default_solver(::ImplicitMidpoint) = Newton()
default_iguess(::ImplicitMidpoint) = HermiteExtrapolation()

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint})
    # temporary solution, extrapolated to the midpoint of the time step
    ig = (t=sol.t - timestep(int) / 2, q=cache(int).q, q̇=cache(int).v)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v

    # compute midpoint q = q̄ + Δt/2 * x (v = x)
    q .= sol.q .+ (timestep(int) / 2) .* x

    # compute v = v(q) at the midpoint
    equations(int).v(v, sol.t - timestep(int) / 2, q, params)
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:ImplicitMidpoint}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v

    # compute residual b = v - x
    b .= v .- x
end


# Compute stages of implicit midpoint methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitMidpoint}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:AbstractProblemODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
