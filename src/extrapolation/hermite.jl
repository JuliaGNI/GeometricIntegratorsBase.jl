@doc raw"""
# Hermite's Interpolating Polynomials

Implements a two point Hermite inter-/extrapolation function which passes
through the function and its first derivative for the interval ``[0,1]``.
The polynomial is determined by four constraint equations, matching the
function and its derivative at the points ``0`` and ``1``.

Call with one of the following methods
```julia
extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, t, x, HermiteExtrapolation())
extrapolate!(t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, t, x, ẋ, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, v, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, ẋ, v, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, problem, HermiteExtrapolation())
extrapolate!(t₀, x₀, t₁, x₁, t, x, ẋ, problem, HermiteExtrapolation())
```

where

* `t₀`: first  sample time $t_0$
* `x₀`: first  solution value $x_0 = x(t_0)$
* `ẋ₀`: first  vector field value $ẋ_0 = v(t_0, x(t_0))$
* `t₁`: second sample time $t_1$
* `x₁`: second solution value $x_1 = x(t_1)$
* `ẋ₁`: second vector field value $ẋ_1 = v(t_1, x(t_1))$
* `t`:  time $t$ to extrapolate
* `x`:  extrapolated solution value $x(t)$
* `ẋ`:  extrapolated vector field value $ẋ(t)$
* `v`:  function to compute vector field with signature `v(ẋ,t,x)`
* `problem`: [`ODEProblem`](@ref) whose vector field to use

See [`NormalizedHermiteExtrapolation`](@ref) for a version of this
extrapolation that is normalised to the interval ``[0,1]``.


#### Derivation

The interpolation works as follows:
Start by defining the 3rd degree polynomial and its derivative by
```math
\begin{aligned}
g(x) &= a_0 + a_1 x + a_2 x^2 + a_3 x^3 , \\
g'(x) &= a_1 + 2 a_2 x + 3 a_3 x^2 ,
\end{aligned}
```
and apply the constraints
```math
\begin{aligned}
g(0) &= f_0 & & \Rightarrow & a_0 &= f_0 , \\
g(1) &= f_1 & & \Rightarrow & a_0 + a_1 + a_2 + a_3 &= f_1 , \\
g'(0) &= f'_0 & & \Rightarrow & a_1 &= f'_0 , \\
g'(1) &= f'_1 & & \Rightarrow & a_1 + 2 a_2 + 3 a_3 &= f'_1 . \\
\end{aligned}
```
Solving for ``a_0, a_1, a_2, a_3`` leads to
```math
\begin{aligned}
a_0 &= f_0 , &
a_1 &= f'_0 , &
a_2 &= - 3 f_0 + 3 f_1 - 2 f'_0 - f'_1 , &
a_3 &= 2 f_0 - 2 f_1 + f'_0 + f'_1 ,
\end{aligned}
```
so that the polynomial ``g(x)`` reads
```math
g(x) = f_0 + f'_0 x + (- 3 f_0 + 3 f_1 - 2 f'_0 - f'_1) x^2 + (2 f_0 - 2 f_1 + f'_0 + f'_1) x^3 .
```
The function and derivative values can be factored out, so that ``g(x)`` can be rewritten as
```math
g(x) = f_0 (1 - 3 x^2 + 2 x^3) + f_1 (3 x^2 - 2 x^3) + f'_0 (x - 2 x^2 + x^3) + f'_1 (- x^2 + x^3) ,
```
or in generic form as
```math
g(x) = f_0 a_0(x) + f_1 a_1(x) + f'_0 b_0(x) + f'_1 b_1(x) ,
```
with basis functions
```math
\begin{aligned}
a_0 (x) &= 1 - 3 x^2 + 2 x^3 , &
b_0 (x) &= x - 2 x^2 + x^3 , \\
a_1 (x) &= 3 x^2 - 2 x^3 , &
b_1 (x) &= - x^2 + x^3 .
\end{aligned}
```
The derivative ``g'(x)`` accordingly reads
```math
g'(x) = f_0 a'_0(x) + f_1 a'_1(x) + f'_0 b'_0(x) + f'_1 b'_1(x) ,
```
with
```math
\begin{aligned}
a'_0 (x) &= - 6 x + 6 x^2 , &
b'_0 (x) &= 1 - 4 x + 3 x^2 , \\
a'_1 (x) &= 6 x - 6 x^2 , &
b'_1 (x) &= - 2 x + 3 x^2 .
\end{aligned}
```
The basis functions ``a_0``and ``a_1`` are associated with the function
values at ``x_0`` and ``x_1``, respectively, while the basis functions
``b_0`` and ``b_1`` are associated with the derivative values at
``x_0`` and ``x_1``.
The basis functions satisfy the following relations,
```math
\begin{aligned}
a_i (x_j) &= \delta_{ij} , &
b_i (x_j) &= 0 , &
a'_i (x_j) &= 0 , &
b'_i (x_j) &= \delta_{ij} , &
i,j &= 0, 1 ,
\end{aligned}
```
where ``\delta_{ij}`` denotes the Kronecker-delta, so that
```math
\begin{aligned}
g(0) &= f_0 , &
g(1) &= f_1 , &
g'(0) &= f'_0 , &
g'(1) &= f'_1 .
\end{aligned}
```
"""
struct HermiteExtrapolation <: Extrapolation end


# Evaluate the Hermite polynomial through the two samples (x₀,ẋ₀) and (x₁,ẋ₁) at the
# normalised time c, that is at t = t₁ + c ⋅ Δt with Δt = t₁ - t₀, so that c = -1
# corresponds to the first and c = 0 to the second sample.
# The coefficients are the basis functions of the derivation above, evaluated at
# s = 1 + c. The derivative values ẋ₀ and ẋ₁ are with respect to the time t, hence
# the scaling by Δt. Passing Δt = 1 amounts to a fully normalised interpolation.
function _extrapolate_hermite!(x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    c::TT, Δt::TT, xᵢ::AbstractArray{DT}) where {DT,TT}

    a₁ = 1 - 3c^2 - 2c^3
    a₀ = 1 - a₁
    b₁ = c * (1 + c)^2
    b₀ = c^2 * (1 + c)
    xᵢ .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* Δt .* ẋ₀ .+ b₁ .* Δt .* ẋ₁

    return xᵢ
end

function _extrapolate_hermite!(x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    c::TT, Δt::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT}) where {DT,TT}

    _extrapolate_hermite!(x₀, ẋ₀, x₁, ẋ₁, c, Δt, xᵢ)

    a₁ = -6c * (1 + c) / Δt
    a₀ = -a₁
    b₁ = (1 + c) * (1 + 3c)
    b₀ = c * (2 + 3c)
    ẋᵢ .= a₀ .* x₀ .+ a₁ .* x₁ .+ b₀ .* ẋ₀ .+ b₁ .* ẋ₁

    return (xᵢ, ẋᵢ)
end


function extrapolate!(t₀::TT, x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    t₁::TT, x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    tᵢ::TT, xᵢ::AbstractArray{DT},
    ::HermiteExtrapolation) where {DT,TT}

    t₀ == t₁ && throw(ArgumentError("t₀ and t₁ in Hermite extrapolation are identical!"))

    Δt::TT = t₁ - t₀

    return _extrapolate_hermite!(x₀, ẋ₀, x₁, ẋ₁, (tᵢ - t₁) / Δt, Δt, xᵢ)
end

function extrapolate!(t₀::TT, x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    t₁::TT, x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    tᵢ::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT},
    ::HermiteExtrapolation) where {DT,TT}

    t₀ == t₁ && throw(ArgumentError("t₀ and t₁ in Hermite extrapolation are identical!"))

    Δt::TT = t₁ - t₀

    return _extrapolate_hermite!(x₀, ẋ₀, x₁, ẋ₁, (tᵢ - t₁) / Δt, Δt, xᵢ, ẋᵢ)
end


function solutionstep!(sol, history, problem::Union{AbstractProblemODE,SODEProblem}, extrap::HermiteExtrapolation; nowarn=false)
    t₀, q₀, q̇₀ = history[2].t, history[2].q, history[2].q̇
    t₁, q₁, q̇₁ = history[1].t, history[1].q, history[1].q̇

    if q₀ == q₁
        nowarn || @warn "Hermite Extrapolation: q's history[1] and history[2] are identical!"
        sol.q .= q₁
        sol.q̇ .= q̇₁
    else
        extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, sol.t, sol.q, sol.q̇, extrap)
    end

    return sol
end

function solutionstep!(sol, history, problem::Union{AbstractProblemPODE,AbstractProblemIODE}, extrap::HermiteExtrapolation; nowarn=false)
    t₀, q₀, v₀, p₀, f₀ = history[2].t, history[2].q, history[2].q̇, history[2].p, history[2].ṗ
    t₁, q₁, v₁, p₁, f₁ = history[1].t, history[1].q, history[1].q̇, history[1].p, history[1].ṗ

    if q₀ == q₁
        nowarn || @warn "Hermite Extrapolation: q's history[1] and history[2] are identical!"
        sol.q .= q₁
        sol.q̇ .= v₁
    else
        extrapolate!(t₀, q₀, v₀, t₁, q₁, v₁, sol.t, sol.q, sol.q̇, extrap)
    end

    if p₀ == p₁
        nowarn || @warn "Hermite Extrapolation: p's history[1] and history[2] are identical!"
        sol.p .= p₁
        sol.ṗ .= f₁
    else
        extrapolate!(t₀, p₀, f₀, t₁, p₁, f₁, sol.t, sol.p, sol.ṗ, extrap)
    end

    return sol
end


@doc raw"""
# Normalized Hermite's Interpolating Polynomials

Implements the same two point Hermite inter-/extrapolation as
[`HermiteExtrapolation`](@ref), but normalised to the interval ``[0,1]``, so that
no sample times need to be passed. Instead of the time ``t`` to extrapolate to,
the normalised time ``c`` is passed, which corresponds to
```math
t = t_1 + c \, \Delta t ,
\qquad
\Delta t = t_1 - t_0 ,
```
where ``t_0`` and ``t_1`` are the times of the first and second sample.
Hence ``c = -1`` reproduces the first and ``c = 0`` the second sample.

As the interpolation is normalised, all derivative values are with respect to the
normalised time ``c``, that is they are scaled by ``\Delta t`` compared to the
vector field values of [`HermiteExtrapolation`](@ref).

Call with one of the following methods
```julia
extrapolate!(x₀, ẋ₀, x₁, ẋ₁, c, x, NormalizedHermiteExtrapolation())
extrapolate!(x₀, ẋ₀, x₁, ẋ₁, c, x, ẋ, NormalizedHermiteExtrapolation())
```

where

* `x₀`: first  solution value $x_0 = x(t_0)$
* `ẋ₀`: first  derivative value $ẋ_0 = \Delta t \, v(t_0, x(t_0))$
* `x₁`: second solution value $x_1 = x(t_1)$
* `ẋ₁`: second derivative value $ẋ_1 = \Delta t \, v(t_1, x(t_1))$
* `c`:  normalised time $c$ to extrapolate, corresponding to $t = t_1 + c \, \Delta t$
* `x`:  extrapolated solution value $x(t)$
* `ẋ`:  extrapolated derivative value $\Delta t \, ẋ(t)$


#### Basis functions

Substituting ``x \to 1 + c`` into the basis functions ``a_i(x)`` and ``b_i(x)``
derived for [`HermiteExtrapolation`](@ref) gives
```math
\begin{aligned}
a_0 (c) &= 3 c^2 + 2 c^3 , &
b_0 (c) &= c^2 (1 + c) , \\
a_1 (c) &= 1 - 3 c^2 - 2 c^3 , &
b_1 (c) &= c (1 + c)^2 ,
\end{aligned}
```
so that
```math
g(c) = x_0 \, a_0(c) + x_1 \, a_1(c) + ẋ_0 \, b_0(c) + ẋ_1 \, b_1(c) ,
```
with derivatives
```math
\begin{aligned}
a'_0 (c) &= 6 c (1 + c) , &
b'_0 (c) &= c (2 + 3 c) , \\
a'_1 (c) &= - 6 c (1 + c) , &
b'_1 (c) &= (1 + c) (1 + 3 c) .
\end{aligned}
```
The basis functions satisfy
```math
\begin{aligned}
g(-1) &= x_0 , &
g(0) &= x_1 , &
g'(-1) &= ẋ_0 , &
g'(0) &= ẋ_1 .
\end{aligned}
```
"""
struct NormalizedHermiteExtrapolation <: Extrapolation end


function extrapolate!(x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    cᵢ::TT, xᵢ::AbstractArray{DT},
    ::NormalizedHermiteExtrapolation) where {DT,TT}

    return _extrapolate_hermite!(x₀, ẋ₀, x₁, ẋ₁, cᵢ, one(TT), xᵢ)
end

function extrapolate!(x₀::AbstractArray{DT}, ẋ₀::AbstractArray{DT},
    x₁::AbstractArray{DT}, ẋ₁::AbstractArray{DT},
    cᵢ::TT, xᵢ::AbstractArray{DT}, ẋᵢ::AbstractArray{DT},
    ::NormalizedHermiteExtrapolation) where {DT,TT}

    return _extrapolate_hermite!(x₀, ẋ₀, x₁, ẋ₁, cᵢ, one(TT), xᵢ, ẋᵢ)
end

function solutionstep!(sol, history, problem::Union{AbstractProblemODE,SODEProblem}, ::NormalizedHermiteExtrapolation; nowarn=false)
    t₀, q₀, q̇₀ = history[2].t, history[2].q, history[2].q̇
    t₁, q₁, q̇₁ = history[1].t, history[1].q, history[1].q̇

    Δt = timestep(problem)
    cᵢ = (sol.t - t₁) / Δt

    if q₀ == q₁
        nowarn || @warn "Normalized Hermite Extrapolation: q's history[1] and history[2] are identical!"
        sol.q .= q₁
        sol.q̇ .= q̇₁
    else
        _extrapolate_hermite!(q₀, q̇₀, q₁, q̇₁, cᵢ, Δt, sol.q, sol.q̇)
    end

    return sol
end

function solutionstep!(sol, history, problem::Union{AbstractProblemPODE,AbstractProblemIODE}, ::NormalizedHermiteExtrapolation; nowarn=false)
    t₀, q₀, v₀, p₀, f₀ = history[2].t, history[2].q, history[2].q̇, history[2].p, history[2].ṗ
    t₁, q₁, v₁, p₁, f₁ = history[1].t, history[1].q, history[1].q̇, history[1].p, history[1].ṗ

    Δt = timestep(problem)
    cᵢ = (sol.t - t₁) / Δt

    if q₀ == q₁
        nowarn || @warn "Normalized Hermite Extrapolation: q's history[1] and history[2] are identical!"
        sol.q .= q₁
        sol.q̇ .= v₁
    else
        _extrapolate_hermite!(q₀, v₀, q₁, v₁, cᵢ, Δt, sol.q, sol.q̇)
    end

    if p₀ == p₁
        nowarn || @warn "Normalized Hermite Extrapolation: p's history[1] and history[2] are identical!"
        sol.p .= p₁
        sol.ṗ .= f₁
    else
        _extrapolate_hermite!(p₀, f₀, p₁, f₁, cᵢ, Δt, sol.p, sol.ṗ)
    end

    return sol
end
