"""
Implicit Euler Method.

"""
struct ImplicitEuler <: ODEMethod end
# $(reference(Val(:ImplicitEuler)))

isexplicit(method::ImplicitEuler) = false
isimplicit(method::ImplicitEuler) = true
issymmetric(method::ImplicitEuler) = false
issymplectic(method::ImplicitEuler) = false


@doc raw"""
Implicit Euler integrator cache.
"""
struct ImplicitEulerCache{DT} <: ODEIntegratorCache{DT}
    x::Vector{DT}
    q::Vector{DT}
    v::Vector{DT}
    v̄::Vector{DT}

    function ImplicitEulerCache{DT}(ics) where {DT}
        x = zeros(DT, length(vec(ics.q)))
        q = zeros(DT, axes(ics.q))
        v = zeros(DT, axes(ics.q))
        v̄ = zeros(DT, axes(ics.q))
        new(x, q, v, v̄)
    end
end

nlsolution(cache::ImplicitEulerCache) = cache.x

function Cache{ST}(problem::AbstractProblem, method::ImplicitEuler; kwargs...) where {ST}
    ImplicitEulerCache{ST}(initial_conditions(problem); kwargs...)
end

@inline CacheType(ST, ::AbstractProblem, ::ImplicitEuler) = ImplicitEulerCache{ST}


solversize(problem::AbstractProblemODE, ::ImplicitEuler) = length(vec(initial_conditions(problem).q))

default_solver(::ImplicitEuler) = Newton()
default_iguess(::ImplicitEuler) = HermiteExtrapolation()

function initial_guess!(sol, history, params, int::GeometricIntegrator{<:ImplicitEuler})
    # temporary solution
    ig = (t=sol.t, q=cache(int).q, q̇=cache(int).v)

    # compute initial guess
    solutionstep!(ig, history, problem(int), iguess(int))

    # assemble initial guess for nonlinear solver solution vector
    nlsolution(int) .= ig.q̇
end

function components!(x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    q = cache(int, ST).q
    v = cache(int, ST).v
    v̄ = cache(int, ST).v̄

    # compute q = q̄ + Δt * x (v = x)
    v̄ .= x
    q .= sol.q .+ timestep(int) .* v̄

    # compute v = v(q)
    equations(int).v(v, sol.t, q, params)
end


function residual!(b::AbstractVector{ST}, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    # get cache for internal stages
    v = cache(int, ST).v
    v̄ = cache(int, ST).v̄

    # compute b = - (v-v)
    b .= v .- v̄
end


# Compute stages of implicit Euler methods.
function residual!(b::AbstractVector{ST}, x::AbstractVector{ST}, sol, params, int::GeometricIntegrator{<:ImplicitEuler}) where {ST}
    @assert axes(x) == axes(b)

    # copy previous solution from solstep to cache
    # reset!(cache(int, ST), sol...)

    # compute stages from nonlinear solver solution x
    components!(x, sol, params, int)

    # compute residual vector
    residual!(b, int)
end


function update!(sol, params, x::AbstractVector{DT}, int::GeometricIntegrator{<:ImplicitEuler}) where {DT}
    # copy previous solution from solstep to cache
    # reset!(cache(int, DT), sol...)

    # compute vector field at internal stages
    components!(x, sol, params, int)

    # compute final update
    sol.q .+= timestep(int) .* cache(int, DT).v
end


function integrate_step!(sol, history, params, int::GeometricIntegrator{<:ImplicitEuler,<:AbstractProblemODE})
    # call nonlinear solver
    solve!(nlsolution(int), solver(int), solverstate(int), (sol, params, int))

    # print solver status
    # println(status(solver))

    # check if solution contains NaNs or error bounds are violated
    # println(meets_stopping_criteria(status(solver)))

    # compute final update
    update!(sol, params, nlsolution(int), int)
end
