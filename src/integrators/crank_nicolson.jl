@doc raw"""
Crank-Nicolson Method, also known as the trapezoidal rule.

For an ordinary differential equation ``\dot{q} = v(t,q)``, that is an `ODEProblem`,
the method reads
```math
q_{n+1} = q_{n} + \frac{h}{2} \, \big[ v (t_{n}, q_{n}) + v (t_{n+1}, q_{n+1}) \big]
```
The nonlinear solver solves for the vector field at the new time step,
``V = v(t_{n+1}, q_{n+1})``, while ``\bar{V} = v(t_{n}, q_{n})`` is computed once per time step.
The method is symmetric and second order, but not symplectic (it is conjugate to a symplectic
method).

For a partitioned differential equation, that is a `PODEProblem` or `HODEProblem`,
```math
\dot{q} = v (t, q, p) ,
\qquad
\dot{p} = f (t, q, p) ,
```
the trapezoidal rule is applied to both components,
```math
\begin{aligned}
q_{n+1} &= q_{n} + \frac{h}{2} \, ( \bar{V} + V ) , &
p_{n+1} &= p_{n} + \frac{h}{2} \, ( \bar{F} + F ) ,
\end{aligned}
```
where ``\bar{V} = v(t_{n}, q_{n}, p_{n})`` and ``\bar{F} = f(t_{n}, q_{n}, p_{n})`` are computed
once per time step, as they are given by the solution at the beginning of the time step, while
``V = v(t_{n+1}, q_{n+1}, p_{n+1})`` and ``F = f(t_{n+1}, q_{n+1}, p_{n+1})`` are solved for. The
solver solution vector holds both of them, so the nonlinear system is twice the size of the one
for an ordinary differential equation.

For an implicit differential equation, that is an `IODEProblem` or `LODEProblem`,
```math
\begin{aligned}
p &= \vartheta (t, q, v) , &
\dot{p} &= f (t, q, v) , &
\dot{q} &= v ,
\end{aligned}
```
the trapezoidal rule is applied to the position and to the momentum alike, which amounts to the
Lobatto IIIA method with two stages,
```math
\begin{aligned}
q_{n+1} &= q_{n} + \frac{h}{2} \, ( \bar{V} + V ) , &
p_{n+1} &= p_{n} + \frac{h}{2} \, \big[ f (t_{n}, q_{n}, \bar{V}) + f (t_{n+1}, q_{n+1}, V) \big] .
\end{aligned}
```
In contrast to the explicit case, the velocity at the beginning of the time step is not given by
a function evaluation but implicitly by ``\vartheta (t_{n}, q_{n}, \bar{V}) = p_{n}``. Both
velocities are therefore solved for simultaneously, so that the nonlinear system reads
```math
\begin{aligned}
0 &= \vartheta (t_{n}, q_{n}, \bar{V}) - p_{n} , \\
0 &= \vartheta (t_{n+1}, q_{n+1}, V) - p_{n} - \frac{h}{2} \, \big[ f (t_{n}, q_{n}, \bar{V}) + f (t_{n+1}, q_{n+1}, V) \big] ,
\end{aligned}
```
and is twice the size of the one for an ordinary differential equation. Whenever ``\vartheta``
is regular, so that the implicit equation is equivalent to a partitioned ordinary differential
equation, this is the same map as the Crank-Nicolson method applied to that equation.
"""
struct CrankNicolson <: ODEMethod end

isexplicit(method::CrankNicolson) = false
isimplicit(method::CrankNicolson) = true
issymmetric(method::CrankNicolson) = true
issymplectic(method::CrankNicolson) = false

# besides ordinary differential equations the method is implemented for partitioned and implicit
# ones, so the corresponding traits are set explicitly, as they cannot follow from the supertype
ispodemethod(::Union{CrankNicolson,Type{<:CrankNicolson}}) = true
ishodemethod(::Union{CrankNicolson,Type{<:CrankNicolson}}) = true
isiodemethod(::Union{CrankNicolson,Type{<:CrankNicolson}}) = true
islodemethod(::Union{CrankNicolson,Type{<:CrankNicolson}}) = true


@doc raw"""
Crank-Nicolson integrator cache.

### Fields

* `x`: nonlinear solver solution vector, holding the vector field ``V = v(t_{n+1}, q_{n+1})``
* `q`: solution at the end of the time step
* `v`: vector field at the end of the time step
* `v̄`: vector field at the beginning of the time step, ``v(t_{n}, q_{n})``, which is constant
  throughout the nonlinear solve
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

function Cache{ST}(problem::ProblemODE, method::CrankNicolson; kwargs...) where {ST}
    CrankNicolsonCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemODE, ::CrankNicolson) = CrankNicolsonCache{ST}


@doc raw"""
Crank-Nicolson integrator cache for partitioned differential equations.

As in [`CrankNicolsonCache`](@ref), the vector fields at the beginning of the time step are
constant during the solve, so they are read from the cache at working precision.

### Fields

* `x`: nonlinear solver solution vector, holding ``V`` and ``F``
* `q`, `p`: solution at the end of the time step
* `v`, `f`: vector fields at the end of the time step
* `v̄`, `f̄`: vector fields at the beginning of the time step, ``v(t_{n}, q_{n}, p_{n})`` and
  ``f(t_{n}, q_{n}, p_{n})``, which are constant throughout the nonlinear solve
"""
struct CrankNicolsonPODECache{DT} <: PODEIntegratorCache{DT}
    x::Vector{DT}

    q::Vector{DT}
    p::Vector{DT}
    v::Vector{DT}
    f::Vector{DT}

    v̄::Vector{DT}
    f̄::Vector{DT}

    function CrankNicolsonPODECache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)) + length(vec(ics.p)))
        q = zeros(DT, axes(ics.q))
        p = zeros(DT, axes(ics.p))
        v = zeros(DT, axes(ics.q))
        f = zeros(DT, axes(ics.p))
        v̄ = zeros(DT, axes(ics.q))
        f̄ = zeros(DT, axes(ics.p))
        new(x, q, p, v, f, v̄, f̄)
    end
end

nlsolution(cache::CrankNicolsonPODECache) = cache.x

function Cache{ST}(problem::ProblemPODE, method::CrankNicolson; kwargs...) where {ST}
    CrankNicolsonPODECache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemPODE, ::CrankNicolson) = CrankNicolsonPODECache{ST}


@doc raw"""
Crank-Nicolson integrator cache for implicit differential equations.

In contrast to [`CrankNicolsonCache`](@ref), nothing is constant during the solve: the velocity
at the beginning of the time step is part of the nonlinear solver solution vector, so every field
is computed from it and has to be accessed through `cache(int, ST)`.

### Fields

* `x`: nonlinear solver solution vector, holding ``\bar{V}`` and ``V``
* `q`: solution at the end of the time step
* `v`, `θ`, `f`: velocity, momentum map and force at the end of the time step
* `v̄`, `θ̄`, `f̄`: velocity, momentum map and force at the beginning of the time step
"""
struct CrankNicolsonIODECache{DT} <: IODEIntegratorCache{DT}
    x::Vector{DT}

    q::Vector{DT}
    v::Vector{DT}
    θ::Vector{DT}
    f::Vector{DT}

    v̄::Vector{DT}
    θ̄::Vector{DT}
    f̄::Vector{DT}

    function CrankNicolsonIODECache{DT}(ics) where {DT}
        x = zeros(DT, 2 * length(vec(ics.q)))
        q = zeros(DT, axes(ics.q))
        v = zeros(DT, axes(ics.q))
        θ = zeros(DT, axes(ics.q))
        f = zeros(DT, axes(ics.q))
        v̄ = zeros(DT, axes(ics.q))
        θ̄ = zeros(DT, axes(ics.q))
        f̄ = zeros(DT, axes(ics.q))
        new(x, q, v, θ, f, v̄, θ̄, f̄)
    end
end

nlsolution(cache::CrankNicolsonIODECache) = cache.x

function Cache{ST}(problem::ProblemIODE, method::CrankNicolson; kwargs...) where {ST}
    CrankNicolsonIODECache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::ProblemIODE, ::CrankNicolson) = CrankNicolsonIODECache{ST}


solversize(::CrankNicolson, problem::ProblemODE) = length(vec(initial_conditions(problem).q))
solversize(::CrankNicolson, problem::ProblemIODE) = 2 * length(vec(initial_conditions(problem).q))
solversize(::CrankNicolson, problem::ProblemPODE) =
    length(vec(initial_conditions(problem).q)) + length(vec(initial_conditions(problem).p))

default_solver(::CrankNicolson) = Newton()
default_iguess(::CrankNicolson) = HermiteExtrapolation()


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE})
    # temporary solution
    ig = (t=sol.t, q=cache(int).q, q̇=cache(int).v)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

@doc raw"""
Compute the stages of the Crank-Nicolson method from the nonlinear solver solution `x`.

Requires `v̄` in the cache at working precision to hold ``v(t_{n}, q_{n})``, which
`integrate_step!` computes at the beginning of every time step. Calling this function before
that would silently use a stale or zero `v̄`.
"""
function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE}) where {ST}
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


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v

    # compute residual b = v - x
    b .= v .- x
end


# Compute stages of Crank-Nicolson methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE}) where {DT}
    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= (timestep(int) / 2) .* (cache(int).v̄ .+ cache(int, DT).v)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemODE})
    # compute vector field at the previous time step
    # this cannot be taken from the vector field stored in the solution step, as that is
    # computed with initialguess(problem).v, which may be a surrogate for equations(int).v
    equations(int).v(cache(int).v̄, sol.t - timestep(int), sol.q, params)

    # call nonlinear solver
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    check_solver_status(solverstatus, int)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE})
    local D = length(cache(int).q)
    local x = nlsolution(int)

    # temporary solution, extrapolated to the end of the time step
    ig = (t=sol.t, q=cache(int).q, p=cache(int).p, q̇=cache(int).v, ṗ=cache(int).f)

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

@doc raw"""
Compute the stages of the Crank-Nicolson method for a partitioned differential equation from the
nonlinear solver solution `x`.

Requires `v̄` and `f̄` in the cache at working precision to hold ``v(t_{n}, q_{n}, p_{n})`` and
``f(t_{n}, q_{n}, p_{n})``, which `integrate_step!` computes at the beginning of every time step.
Calling this function before that would silently use stale or zero values.
"""
function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE}) where {ST}
    q = cache(int, ST).q
    p = cache(int, ST).p
    v = cache(int, ST).v
    f = cache(int, ST).f

    # v̄ and f̄ are constant during the solve, so they are read from the cache at working
    # precision, which also holds them during the automatic differentiation of the residual
    v̄ = cache(int).v̄
    f̄ = cache(int).f̄

    local D = length(q)

    # compute q = q̄ + Δt/2 * (v̄ + x[1:D]) and p = p̄ + Δt/2 * (f̄ + x[D+1:2D]),
    # as the solver solution vector holds the vector fields (v, f) at the end of the time step
    for k in 1:D
        q[k] = sol.q[k] + (timestep(int) / 2) * (v̄[k] + x[k])
        p[k] = sol.p[k] + (timestep(int) / 2) * (f̄[k] + x[D+k])
    end

    # compute v = v(t, q, p) and f = f(t, q, p) at the end of the time step
    equations(int).v(v, sol.t, q, p, params)
    equations(int).f(f, sol.t, q, p, params)
end


function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE}) where {ST}
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


# Compute stages of Crank-Nicolson methods for partitioned differential equations.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, x, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE}) where {DT}
    # compute vector fields at internal stages
    # this has to precede the update, as the stages are computed relative to sol.q and sol.p
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= (timestep(int) / 2) .* (cache(int).v̄ .+ cache(int, DT).v)
    sol.p .+= (timestep(int) / 2) .* (cache(int).f̄ .+ cache(int, DT).f)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemPODE})
    # compute vector fields at the previous time step
    # this cannot be taken from the vector fields stored in the solution step, as those are
    # computed with initialguess(problem).v and .f, which may be surrogates for the vector
    # fields of the equation
    equations(int).v(cache(int).v̄, sol.t - timestep(int), sol.q, sol.p, params)
    equations(int).f(cache(int).f̄, sol.t - timestep(int), sol.q, sol.p, params)

    # call nonlinear solver
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    check_solver_status(solverstatus, int)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end


function initial_guess!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE})
    local D = length(cache(int).q)
    local x = nlsolution(int)

    # the solver variables are the stage velocities determined by ϑ(t, q, V) = p rather than the
    # velocity of the solution, so they are guessed with the velocity the problem provides for
    # exactly that purpose, evaluated at both ends of the time step

    # at the beginning of the time step the solution is known: `sol` still holds it, as it is
    # only updated at the end of `integrate_step!`
    v̄ = cache(int).v̄
    initialguess(problem(int)).v(v̄, sol.t - timestep(int), sol.q, sol.p, params)

    for k in 1:D
        x[k] = v̄[k]
    end

    # at the end of the time step the solution is extrapolated
    ig = (t=sol.t, q=cache(int).q, p=cache(int).θ, q̇=cache(int).v, ṗ=cache(int).f)
    solutionstep!(ig, history, problem(int), iguess(int))
    initialguess(problem(int)).v(ig.q̇, ig.t, ig.q, ig.p, params)

    for k in 1:D
        x[D+k] = ig.q̇[k]
    end
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v
    θ = cache(int, ST).θ
    f = cache(int, ST).f
    v̄ = cache(int, ST).v̄
    θ̄ = cache(int, ST).θ̄
    f̄ = cache(int, ST).f̄

    local D = length(q)
    local t̄ = sol.t - timestep(int)

    # the nonlinear solver solution vector holds the stage velocities at both ends of the step
    for k in 1:D
        v̄[k] = x[k]
        v[k] = x[D+k]
    end

    # compute q = q̄ + Δt/2 * (v̄ + v)
    q .= sol.q .+ (timestep(int) / 2) .* (v̄ .+ v)

    # compute θ̄ = ϑ(t̄, q̄, v̄) and f̄ = f(t̄, q̄, v̄) at the beginning of the time step
    equations(int).ϑ(θ̄, t̄, sol.q, v̄, params)
    equations(int).f(f̄, t̄, sol.q, v̄, params)

    # compute θ = ϑ(t, q, v) and f = f(t, q, v) at the end of the time step
    equations(int).ϑ(θ, sol.t, q, v, params)
    equations(int).f(f, sol.t, q, v, params)
end


function residual!(b::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE}) where {ST}
    # get cache for internal stages
    θ = cache(int, ST).θ
    f = cache(int, ST).f
    θ̄ = cache(int, ST).θ̄
    f̄ = cache(int, ST).f̄

    local D = length(cache(int, ST).q)

    for k in 1:D
        # the velocity at the beginning of the time step is determined by the momentum there
        b[k] = θ̄[k] - sol.p[k]

        # trapezoidal rule for the momentum
        b[D+k] = θ[k] - sol.p[k] - timestep(int) * (f̄[k] + f[k]) / 2
    end
end


# Compute stages of Crank-Nicolson methods for implicit differential equations.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE}) where {ST}
    axes(x) == axes(b) || throw(ArgumentError("x and b must have the same axes"))

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, sol, params, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE}) where {DT}
    # compute vector fields at internal stages
    # this has to precede the update, as the stages are computed relative to sol.q and sol.p
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= (timestep(int) / 2) .* (cache(int, DT).v̄ .+ cache(int, DT).v)
    sol.p .+= (timestep(int) / 2) .* (cache(int, DT).f̄ .+ cache(int, DT).f)
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:CrankNicolson,<:ProblemIODE})
    # call nonlinear solver
    # in contrast to the explicit case, nothing is precomputed at the beginning of the time
    # step: the velocity there is determined implicitly and is part of the solver solution
    solverstatus = solve_with_status!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))
    check_solver_status(solverstatus, int)

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
