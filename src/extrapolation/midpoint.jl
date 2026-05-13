@doc raw"""
Midpoint extrapolation method with arbitrary order p.

For an [`ODEProblem`](@ref), this solves the ordinary differential equation

```math
\begin{aligned}
\dot{x} &= v(t, x) , &
x(t_0) &= x_0 ,
\end{aligned}
```

for $x_1 = x(t_1)$, and is called with

```julia
extrapolate!(t₀, x₀, t₁, x₁, ::ODEProblem, MidpointExtrapolation(s))
```

where

* `t₀`: initial time
* `x₀`: initial value $x_0 = x(t_0)$
* `t₁`: final   time
* `x₁`: final   value $x_1 = x(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)


For a [`PODEProblem`](@ref) or [`HODEProblem`](@ref),
this solves the partitioned ordinary differential equation

```math
\begin{aligned}
\dot{q} &= v(t, q, p) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, p) , &
p(t_0) &= p_0 ,
\end{aligned}
```

for $q_1 = q(t_1)$ and $p_1 = p(t_1)$, and is called with

```julia
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::PODEProblem, MidpointExtrapolation(s))
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::HODEProblem, MidpointExtrapolation(s))
```

where

* `t₀`: initial time
* `q₀`: initial position $q_0 = q(t_0)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `t₁`: final   time
* `q₁`: final   position $q_1 = q(t_1)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)


Similarly, for a [`IODEProblem`](@ref) or [`LODEProblem`](@ref),
this solves the explicit dynamical equation

```math
\begin{aligned}
\dot{q} &= v(t, q) , &
q(t_0) &= q_0 , \\
\dot{p} &= f(t, q, v) , &
p(t_0) &= p_0 ,
\end{aligned}
```

corresponding to the implicit problem, for $q_1 = q(t_1)$ and $p_1 = p(t_1)$, and is called with

    ```julia
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::IODEProblem, MidpointExtrapolation(s))
extrapolate!(t₀, q₀, p₀, t₁, q₁, p₁, ::LODEProblem, MidpointExtrapolation(s))
```

where

* `t₀`: initial time
* `q₀`: initial position $q_0 = q(t_0)$
* `p₀`: initial momentum $p_0 = p(t_0)$
* `t₁`: final   time
* `q₁`: final   position $q_1 = q(t_1)$
* `p₁`: final   momentum $p_1 = p(t_1)$
* `s`:  number of interpolations (order $p=2s+2$)

"""
struct MidpointExtrapolation{TT} <: Extrapolation
    s::Int
    Δt::TT
    MidpointExtrapolation(s=default_extrapolation_stages, Δt=1.0) = new{typeof(Δt)}(s, Δt)
end

function extrapolate!(
    t₀::TT, x₀::AbstractArray{DT},
    t₁::TT, x₁::AbstractArray{DT},
    problem::Union{AbstractProblemODE,SODEProblem},
    extrap::MidpointExtrapolation) where {DT,TT}

    axes(x₀) == axes(x₁) || throw(ArgumentError("x₀ and x₁ must have the same axes"))

    F = [2i * one(TT) for i in 1:extrap.s+1]
    σ = (t₁ - t₀) ./ F
    σ² = σ .^ 2
    pts = [zero(x₀) for _ in 1:extrap.s+1]

    xᵢ₁ = zero(x₀)
    xᵢ₂ = zero(x₀)
    xᵢₜ = zero(x₀)
    vᵢ = zero(x₀)
    v₀ = zero(x₀)

    initialguess(problem).v(v₀, t₀, x₀, parameters(problem))

    for i in eachindex(pts)
        xᵢ₁ .= x₀
        xᵢ₂ .= x₀ .+ σ[i] .* v₀
        for k in 1:(F[i]-1)
            initialguess(problem).v(vᵢ, t₀ + k * σ[i], xᵢ₂, parameters(problem))
            xᵢₜ .= xᵢ₁ .+ 2σ[i] .* vᵢ
            xᵢ₁ .= xᵢ₂
            xᵢ₂ .= xᵢₜ
        end
        pts[i] .= xᵢ₂
    end

    aitken_neville!(x₁, zero(TT), σ², pts)

    return x₁
end

function _extrapolate!(newsol, oldsol, problem::Union{AbstractProblemODE,SODEProblem}, extrap::MidpointExtrapolation)
    extrapolate!(oldsol.t, oldsol.q, newsol.t, newsol.q, problem, extrap)
    return newsol
end

function solutionstep!(sol, history, problem::Union{AbstractProblemODE,SODEProblem}, extrap::MidpointExtrapolation)
    extrapolate!(history[1].t, history[1].q, sol.t, sol.q, problem, extrap)
    initialguess(problem).v(sol.q̇, sol.t, sol.q, parameters(problem))
    # update_vectorfields!(sol, problem)
    return sol
end


function extrapolate!(t₀::TT, q₀::AbstractVector{DT}, p₀::AbstractVector{DT},
    t₁::TT, q₁::AbstractVector{DT}, p₁::AbstractVector{DT},
    problem::AbstractProblemPODE,
    extrap::MidpointExtrapolation) where {DT,TT}

    axes(q₀) == axes(q₁) || throw(ArgumentError("q₀ and q₁ must have the same axes"))
    axes(p₀) == axes(p₁) || throw(ArgumentError("p₀ and p₁ must have the same axes"))

    F = [2i * one(TT) for i in 1:extrap.s+1]
    σ = (t₁ - t₀) ./ F
    σ2 = σ .^ 2

    qts = [zero(q₀) for _ in 1:extrap.s+1]
    pts = [zero(p₀) for _ in 1:extrap.s+1]

    qᵢ₁ = zero(q₀)
    qᵢ₂ = zero(q₀)
    qᵢₜ = zero(q₀)

    pᵢ₁ = zero(p₀)
    pᵢ₂ = zero(p₀)
    pᵢₜ = zero(p₀)

    v₀ = zero(q₀)
    vᵢ = zero(q₀)

    f₀ = zero(p₀)
    fᵢ = zero(p₀)

    initialguess(problem).v(v₀, t₀, q₀, p₀, parameters(problem))
    initialguess(problem).f(f₀, t₀, q₀, p₀, parameters(problem))

    for i in 1:extrap.s+1
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for k in 1:(F[i]-1)
            initialguess(problem).v(vᵢ, t₀ + k * σ[i], qᵢ₂, pᵢ₂, parameters(problem))
            initialguess(problem).f(fᵢ, t₀ + k * σ[i], qᵢ₂, pᵢ₂, parameters(problem))
            qᵢₜ .= qᵢ₁ .+ 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ .+ 2σ[i] .* fᵢ
            pᵢ₁ .= pᵢ₂
            pᵢ₂ .= pᵢₜ
        end
        qts[i] .= qᵢ₂
        pts[i] .= pᵢ₂
    end

    aitken_neville!(q₁, zero(TT), σ2, qts)
    aitken_neville!(p₁, zero(TT), σ2, pts)

    return (q₁, p₁)
end

function _extrapolate!(newsol, oldsol, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    extrapolate!(oldsol.t, oldsol.q, oldsol.p, newsol.t, newsol.q, newsol.p, problem, extrap)
    return newsol
end

# function solutionstep!(sol, history, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
#     extrapolate!(history.t[1], history.q[1], history.p[1], sol.t, sol.q, sol.p, problem, extrap)
#     initialguess(problem).v(sol.q̇, sol.t, sol.q, sol.p, parameters(problem))
#     initialguess(problem).f(sol.ṗ, sol.t, sol.q, sol.p, parameters(problem))
#     # update_vectorfields!(sol, problem)
#     return sol
# end

function solutionstep!(sol, history, problem::AbstractProblemPODE, extrap::MidpointExtrapolation)
    Δt = sign(sol.t - history[1].t) * extrap.Δt
    t = history[1].t

    tmpsol = copy(history[1])

    while abs(t + Δt) < abs(sol.t)
        prvsol = copy(tmpsol)
        tmpsol.t .= t + Δt
        _extrapolate!(tmpsol, prvsol, problem, extrap)
        t += Δt
    end

    _extrapolate!(sol, tmpsol, problem, extrap)

    initialguess(problem).v(sol.q̇, sol.t, sol.q, sol.p, parameters(problem))
    initialguess(problem).f(sol.ṗ, sol.t, sol.q, sol.p, parameters(problem))
    # update_vectorfields!(sol, problem)

    return sol
end


function extrapolate!(
    t₀::TT, q₀::AbstractArray{DT}, p₀::AbstractArray{DT},
    t₁::TT, q₁::AbstractArray{DT}, p₁::AbstractArray{DT},
    problem::AbstractProblemIODE,
    extrap::MidpointExtrapolation) where {DT,TT}

    axes(q₀) == axes(q₁) || throw(ArgumentError("q₀ and q₁ must have the same axes"))
    axes(p₀) == axes(p₁) || throw(ArgumentError("p₀ and p₁ must have the same axes"))

    F = [2i * one(TT) for i in 1:extrap.s+1]
    σ = (t₁ - t₀) ./ F
    σ2 = σ .^ 2

    qts = [zero(q₀) for _ in 1:extrap.s+1]
    pts = [zero(p₀) for _ in 1:extrap.s+1]

    qᵢ₁ = zero(q₀)
    qᵢ₂ = zero(q₀)
    qᵢₜ = zero(q₀)

    pᵢ₁ = zero(p₀)
    pᵢ₂ = zero(p₀)
    pᵢₜ = zero(p₀)

    v₀ = zero(q₀)
    vᵢ = zero(q₀)

    f₀ = zero(p₀)
    fᵢ = zero(p₀)

    initialguess(problem).v(v₀, t₀, q₀, p₀, parameters(problem))
    initialguess(problem).f(f₀, t₀, q₀, v₀, parameters(problem))

    for i in 1:extrap.s+1
        qᵢ₁ .= q₀
        qᵢ₂ .= q₀ .+ σ[i] .* v₀
        pᵢ₁ .= p₀
        pᵢ₂ .= p₀ .+ σ[i] .* f₀
        for k in 1:(F[i]-1)
            initialguess(problem).v(vᵢ, t₀ + k * σ[i], qᵢ₂, pᵢ₂, parameters(problem))
            initialguess(problem).f(fᵢ, t₀ + k * σ[i], qᵢ₂, vᵢ, parameters(problem))
            qᵢₜ .= qᵢ₁ .+ 2σ[i] .* vᵢ
            qᵢ₁ .= qᵢ₂
            qᵢ₂ .= qᵢₜ
            pᵢₜ .= pᵢ₁ .+ 2σ[i] .* fᵢ
            pᵢ₁ .= pᵢ₂
            pᵢ₂ .= pᵢₜ
        end
        qts[i] .= qᵢ₂
        pts[i] .= pᵢ₂
    end

    aitken_neville!(q₁, zero(TT), σ2, qts)
    aitken_neville!(p₁, zero(TT), σ2, pts)

    return (q₁, p₁)
end

function _extrapolate!(newsol, oldsol, problem::AbstractProblemIODE, extrap::MidpointExtrapolation)
    extrapolate!(oldsol.t, oldsol.q, oldsol.p, newsol.t, newsol.q, newsol.p, problem, extrap)
    return newsol
end

function solutionstep!(sol, history, problem::AbstractProblemIODE, extrap::MidpointExtrapolation)
    extrapolate!(history[1].t, history[1].q, history[1].p, sol.t, sol.q, sol.p, problem, extrap)
    initialguess(problem).v(sol.q̇, sol.t, sol.q, sol.p, parameters(problem))
    initialguess(problem).f(sol.ṗ, sol.t, sol.q, sol.q̇, parameters(problem))
    # update_vectorfields!(sol, problem)
    return sol
end


# TODO: Revise this function! Adapt interface or merge functionality into solutionstep! and remove.
function extrapolate!(newsol, oldsol, problem::GeometricProblem, extrap::MidpointExtrapolation)
    Δt = sign(newsol.t - oldsol.t) * extrap.Δt
    t = oldsol.t
    tmpsol = copy(oldsol)

    while abs(t + Δt) < abs(newsol.t)
        prvsol = copy(tmpsol)
        tmpsol.t .= t + Δt
        _extrapolate!(tmpsol, prvsol, problem, extrap)
        t += Δt
    end

    _extrapolate!(newsol, tmpsol, problem, extrap)

    return newsol
end
