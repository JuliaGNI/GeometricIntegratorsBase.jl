@doc raw"""
Implicit Midpoint Method.

For an ordinary differential equation ``\dot{q} = v(t,q)``, that is an `ODEProblem`,
the method reads
```math
q_{n+1} = q_{n} + h \, v \bigg( t_{n} + \frac{h}{2} , \frac{q_{n} + q_{n+1}}{2} \bigg)
```
The nonlinear solver solves for the stage vector field ``V = v(t_{n} + h/2, q_{n} + h V / 2)``,
so that the update reads ``q_{n+1} = q_{n} + h V``.

For a partitioned differential equation, that is a `PODEProblem` or `HODEProblem`,
```math
\dot{q} = v (t, q, p) ,
\qquad
\dot{p} = f (t, q, p) ,
```
the same quadrature is applied to both components. With the stage time ``\tilde{t} = t_{n} + h/2``
and the midpoints ``Q = q_{n} + h V / 2`` and ``P = p_{n} + h F / 2``, the nonlinear solver solves
```math
V = v (\tilde{t}, Q, P) ,
\qquad
F = f (\tilde{t}, Q, P) ,
```
for the stage vector fields ``V`` and ``F``, and the updates read
```math
q_{n+1} = q_{n} + h \, V ,
\qquad
p_{n+1} = p_{n} + h \, F .
```
Both components are coupled through ``v`` and ``f``, so the solver solution vector holds both
stage vector fields and the nonlinear system is twice the size of the one for an ordinary
differential equation. Even for a separable Hamiltonian the method does not decouple into two
explicit substeps the way the symplectic Euler methods do: it is the Gauss method with a single
stage, applied to the partitioned system.

For an implicit differential equation, that is an `IODEProblem` or `LODEProblem`,
```math
\begin{aligned}
p &= \vartheta (t, q, v) , &
\dot{p} &= f (t, q, v) , &
\dot{q} &= v ,
\end{aligned}
```
the same quadrature is applied to the momentum map and to the force, which amounts to the Gauss
method with a single stage. With the stage time ``\tilde{t} = t_{n} + h/2`` and the midpoint
``Q = q_{n} + h V / 2``, the nonlinear solver solves
```math
\vartheta (\tilde{t}, Q, V) = p_{n} + \frac{h}{2} \, f (\tilde{t}, Q, V)
```
for the stage velocity ``V``, and the updates read
```math
\begin{aligned}
q_{n+1} &= q_{n} + h \, V , &
p_{n+1} &= p_{n} + h \, f (\tilde{t}, Q, V) .
\end{aligned}
```
Whenever ``\vartheta`` is regular, so that the implicit equation is equivalent to a partitioned
ordinary differential equation, this is the same map as the implicit midpoint method applied to
that equation. The solver solution vector holds ``V`` in both cases, so the nonlinear system has
the same size as for an ordinary differential equation.
"""
struct ImplicitMidpoint <: ODEMethod end

isexplicit(method::ImplicitMidpoint) = false
isimplicit(method::ImplicitMidpoint) = true
issymmetric(method::ImplicitMidpoint) = true
issymplectic(method::ImplicitMidpoint) = true

# besides ordinary differential equations the method is implemented for partitioned and implicit
# ones, so the corresponding traits are set explicitly, as they cannot follow from the supertype
ispodemethod(::Union{ImplicitMidpoint,Type{<:ImplicitMidpoint}}) = true
ishodemethod(::Union{ImplicitMidpoint,Type{<:ImplicitMidpoint}}) = true
isiodemethod(::Union{ImplicitMidpoint,Type{<:ImplicitMidpoint}}) = true
islodemethod(::Union{ImplicitMidpoint,Type{<:ImplicitMidpoint}}) = true


@doc raw"""
Implicit midpoint integrator cache.

### Fields

* `x`: nonlinear solver solution vector, holding the stage vector field ``V``
* `q`: midpoint of the time step, ``Q = q_{n} + h V / 2``
* `v`: stage vector field ``V``
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

function Cache{ST}(problem::ProblemODE, method::ImplicitMidpoint; kwargs...) where {ST}
    ImplicitMidpointCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemODE, ::ImplicitMidpoint) = ImplicitMidpointCache{ST}


@doc raw"""
Implicit midpoint integrator cache for partitioned differential equations.

### Fields

* `x`: nonlinear solver solution vector, holding the stage vector fields ``V`` and ``F``
* `q`: midpoint of the time step, ``Q = q_{n} + h V / 2``
* `p`: midpoint of the time step, ``P = p_{n} + h F / 2``
* `v`: stage vector field ``V``
* `f`: stage vector field ``F``
"""
struct ImplicitMidpointPODECache{DT} <: PODEIntegratorCache{DT}
    x::Vector{DT}
    q::Vector{DT}
    p::Vector{DT}
    v::Vector{DT}
    f::Vector{DT}

    function ImplicitMidpointPODECache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)) + length(vec(ics.p)))
        q = zeros(DT, axes(ics.q))
        p = zeros(DT, axes(ics.p))
        v = zeros(DT, axes(ics.q))
        f = zeros(DT, axes(ics.p))
        new(x, q, p, v, f)
    end
end

nlsolution(cache::ImplicitMidpointPODECache) = cache.x

function Cache{ST}(problem::ProblemPODE, method::ImplicitMidpoint; kwargs...) where {ST}
    ImplicitMidpointPODECache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemPODE, ::ImplicitMidpoint) = ImplicitMidpointPODECache{ST}


@doc raw"""
Implicit midpoint integrator cache for implicit differential equations.

### Fields

* `x`: nonlinear solver solution vector, holding the stage velocity ``V``
* `q`: midpoint of the time step, ``Q = q_{n} + h V / 2``
* `v`: stage velocity ``V``
* `θ`: momentum map ``\vartheta`` evaluated at the stage
* `f`: force ``f`` evaluated at the stage
"""
struct ImplicitMidpointIODECache{DT} <: IODEIntegratorCache{DT}
    x::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    function ImplicitMidpointIODECache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)))
        q = zeros(DT, axes(ics.q))
        v = zeros(DT, axes(ics.q))
        θ = zeros(DT, axes(ics.q))
        f = zeros(DT, axes(ics.q))
        new(x, q, v, θ, f)
    end
end

nlsolution(cache::ImplicitMidpointIODECache) = cache.x

function Cache{ST}(problem::ProblemIODE, method::ImplicitMidpoint; kwargs...) where {ST}
    ImplicitMidpointIODECache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemIODE, ::ImplicitMidpoint) = ImplicitMidpointIODECache{ST}


# the solver solves for the stage vector field in the case of an ordinary and for the stage
# velocity in the case of an implicit differential equation, so the size is the same for both
solversize(::ImplicitMidpoint, problem::Union{ProblemODE,ProblemIODE}) = length(vec(initial_conditions(problem).q))

# for a partitioned differential equation both stage vector fields are solved for, so the
# nonlinear system is twice as large
solversize(::ImplicitMidpoint, problem::ProblemPODE) =
    length(vec(initial_conditions(problem).q)) + length(vec(initial_conditions(problem).p))

default_solver(::ImplicitMidpoint) = Newton()
default_iguess(::ImplicitMidpoint) = HermiteExtrapolation()


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE})
    # temporary solution, extrapolated to the midpoint of the time step
    ig = (t=sol.t - timestep(int) / 2, q=cache(int).q, q̇=cache(int).v)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v

    # compute midpoint q = q̄ + Δt/2 * x (v = x)
    q .= sol.q .+ (timestep(int) / 2) .* x

    # compute v = v(q) at the midpoint
    equations(int).v(v, sol.t - timestep(int) / 2, q, params)
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v

    # compute residual b = v - x
    b .= v .- x
end


# Compute stages of implicit midpoint methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE})
    local D = length(cache(int).q)
    local x = nlsolution(int)

    # temporary solution, extrapolated to the midpoint of the time step
    ig = (t=sol.t - timestep(int) / 2, q=cache(int).q, p=cache(int).p, q̇=cache(int).v, ṗ=cache(int).f)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    # in contrast to the implicit case, the solver variables are the vector fields of the solution
    # itself, so the extrapolated q̇ and ṗ are guesses for them as they are
    for k in 1:D
        x[k] = ig.q̇[k]
        x[D+k] = ig.ṗ[k]
    end
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE}) where {ST}
    q = cache(int, ST).q
    p = cache(int, ST).p
    v = cache(int, ST).v
    f = cache(int, ST).f

    local D = length(q)

    # stage time at the midpoint of the time step
    t̃ = sol.t - timestep(int) / 2

    # compute midpoints q = q̄ + Δt/2 * x[1:D] and p = p̄ + Δt/2 * x[D+1:2D],
    # as the solver solution vector holds the stage vector fields (v, f)
    for k in 1:D
        q[k] = sol.q[k] + (timestep(int) / 2) * x[k]
        p[k] = sol.p[k] + (timestep(int) / 2) * x[D+k]
    end

    # compute v = v(t̃, q, p) and f = f(t̃, q, p) at the midpoint
    equations(int).v(v, t̃, q, p, params)
    equations(int).f(f, t̃, q, p, params)
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v
    f = cache(int, ST).f

    local D = length(cache(int, ST).q)

    # compute residual b = (v, f) - x
    for k in 1:D
        b[k] = v[k] - x[k]
        b[D+k] = f[k] - x[D+k]
    end
end


# Compute stages of implicit midpoint methods for partitioned differential equations.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE}) where {DT}
    # compute vector fields at internal stages
    # this has to precede the update, as the stages are computed relative to sol.q and sol.p
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
    sol.p .+= timestep(int) .* cache(int, DT).f
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemPODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE})
    # temporary solution, extrapolated to the midpoint of the time step
    ig = (t=sol.t - timestep(int) / 2, q=cache(int).q, p=cache(int).θ, q̇=cache(int).v, ṗ=cache(int).f)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # the extrapolated q̇ is a guess for the velocity of the solution, whereas the solver
    # variable is the stage velocity determined by ϑ(t̃, Q, V) = P, so refine the guess with
    # the velocity the problem provides for exactly that purpose
    initialguess(problem(int)).v(ig.q̇, ig.t, ig.q, ig.p, params)

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v
    θ = cache(int, ST).θ
    f = cache(int, ST).f

    # stage time at the midpoint of the time step
    t̃ = sol.t - timestep(int) / 2

    # the nonlinear solver solution vector holds the stage velocity
    v .= x

    # compute midpoint q = q̄ + Δt/2 * v
    q .= sol.q .+ (timestep(int) / 2) .* v

    # compute θ = ϑ(t̃, q, v) and f = f(t̃, q, v) at the midpoint
    equations(int).ϑ(θ, t̃, q, v, params)
    equations(int).f(f, t̃, q, v, params)
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE}) where {ST}
    # get cache for internal stages
    θ = cache(int, ST).θ
    f = cache(int, ST).f

    # compute residual b = ϑ - p̄ - Δt/2 * f
    b .= θ .- sol.p .- (timestep(int) / 2) .* f
end


# Compute stages of implicit midpoint methods for implicit differential equations.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE}) where {DT}
    # compute vector fields at internal stages
    # this has to precede the update, as the stages are computed relative to sol.q and sol.p
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
    sol.p .+= timestep(int) .* cache(int, DT).f
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitMidpoint,<:ProblemIODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
