using GeometricIntegratorsBase
using Test

using GeometricIntegratorsBase: extrapolate!, functions, initial_conditions, initialguess, timestep, value

using ..HarmonicOscillator

ode = odeproblem()
pode = podeproblem()
iode = iodeproblem()


# Compute Reference Solution for ODEs

const Δt = timestep(ode)
const t₀ = value(initial_conditions(ode).t)
const t₁ = t₀ + Δt
const t₂ = t₁ + Δt
const t₋ = t₀ - Δt
const tₚ = t₋
const tₙ = t₁
const tᵢ = tₙ

x₀ = initial_conditions(ode).q

k = parameters(ode).k
ω = parameters(ode).ω
A = sqrt(x₀[2]^2 / k + x₀[1]^2)
ϕ = asin(x₀[1] / A)

xₚ = exact_solution(t₀ - Δt, x₀, t₀, parameters(ode))
xₙ = exact_solution(t₀ + Δt, x₀, t₀, parameters(ode))


# Create ODE Solution Arrays

x₁ = zero(x₀)
x₂ = zero(x₀)
xᵢ = zero(x₀)

ẋₚ = zero(x₀)
ẋ₀ = zero(x₀)
ẋ₁ = zero(x₀)
ẋ₂ = zero(x₀)
ẋₙ = zero(x₀)
ẋᵢ = zero(x₀)

functions(ode).v(ẋₚ, tₚ, xₚ, parameters(ode))
functions(ode).v(ẋ₀, t₀, x₀, parameters(ode))
functions(ode).v(ẋₙ, tₙ, xₙ, parameters(ode))


# Hermite Extrapolation

extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, tᵢ, xᵢ, ẋᵢ, HermiteExtrapolation())

# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, ẋₙ, ẋᵢ .- ẋₙ)

@test xᵢ ≈ xₙ atol = 1E-5
@test ẋᵢ ≈ ẋₙ atol = 1E-4

@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, HermiteExtrapolation()) == xᵢ
@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, HermiteExtrapolation()) == (xᵢ, ẋᵢ)

sol = (t=t₁, q=x₁, q̇=ẋ₁)
ref = (t=tᵢ, q=xᵢ, q̇=ẋᵢ)
history = (t=[t₀, tₚ], q=[x₀, xₚ], q̇=[ẋ₀, ẋₚ])
solutionstep!(sol, history, ode, HermiteExtrapolation())
@test sol == ref


# solution and history tuples
sol = (t=tᵢ, q=xᵢ, q̇=ẋᵢ)
history = (t=[t₀], q=[x₀], q̇=[ẋ₀])

# Euler Extrapolation for ODEs

solutionstep!(sol, history, ode, EulerExtrapolation(0))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-2
@test sol.q̇ ≈ ẋₙ atol = 5E-2

solutionstep!(sol, history, ode, EulerExtrapolation(1))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-3
@test sol.q̇ ≈ ẋₙ atol = 5E-3

solutionstep!(sol, history, ode, EulerExtrapolation(2))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-5
@test sol.q̇ ≈ ẋₙ atol = 5E-5

solutionstep!(sol, history, ode, EulerExtrapolation(3))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-6
@test sol.q̇ ≈ ẋₙ atol = 1E-6

solutionstep!(sol, history, ode, EulerExtrapolation(4))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-8
@test sol.q̇ ≈ ẋₙ atol = 1E-8

solutionstep!(sol, history, ode, EulerExtrapolation(5))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-10
@test sol.q̇ ≈ ẋₙ atol = 1E-10


# Midpoint Extrapolation for ODEs

solutionstep!(sol, history, ode, MidpointExtrapolation(0))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-5
@test sol.q̇ ≈ ẋₙ atol = 5E-5

solutionstep!(sol, history, ode, MidpointExtrapolation(1))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-8
@test sol.q̇ ≈ ẋₙ atol = 1E-8

solutionstep!(sol, history, ode, MidpointExtrapolation(2))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-12
@test sol.q̇ ≈ ẋₙ atol = 1E-12

solutionstep!(sol, history, ode, MidpointExtrapolation(3))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-16
@test sol.q̇ ≈ ẋₙ atol = 1E-16

solutionstep!(sol, history, ode, MidpointExtrapolation(4))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-16
@test sol.q̇ ≈ ẋₙ atol = 1E-16

solutionstep!(sol, history, ode, MidpointExtrapolation(5))
# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, xₙ, ẋᵢ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-15
@test sol.q̇ ≈ ẋₙ atol = 1E-15


# Create PODE Solution Arrays

q₀ = initial_conditions(pode).q
p₀ = initial_conditions(pode).p

qᵢ = zero(q₀)
pᵢ = zero(p₀)

q̇ₚ = zero(q₀)
q̇₀ = zero(q₀)
q̇ₙ = zero(q₀)
q̇ᵢ = zero(q₀)

ṗₚ = zero(p₀)
ṗ₀ = zero(p₀)
ṗₙ = zero(p₀)
ṗᵢ = zero(p₀)

sol = (t=tᵢ, q=qᵢ, p=pᵢ, q̇=q̇ᵢ, ṗ=ṗᵢ)


# Compute Reference Solution for PODEs

qₚ = [xₚ[1]]
pₚ = [xₚ[2]]

qₙ = [xₙ[1]]
pₙ = [xₙ[2]]

functions(pode).v(q̇ₚ, tₚ, qₚ, pₚ, parameters(pode))
functions(pode).v(q̇₀, t₀, q₀, p₀, parameters(pode))
functions(pode).v(q̇ₙ, tₙ, qₙ, pₙ, parameters(pode))

functions(pode).f(ṗₚ, tₚ, qₚ, pₚ, parameters(pode))
functions(pode).f(ṗ₀, t₀, q₀, p₀, parameters(pode))
functions(pode).f(ṗₙ, tₙ, qₙ, pₙ, parameters(pode))


# Hermite Extrapolation

history = (t=[t₀, tₚ], q=[q₀, qₚ], p=[p₀, pₚ], q̇=[q̇₀, q̇ₚ], ṗ=[ṗ₀, ṗₚ])

solutionstep!(sol, history, pode, HermiteExtrapolation())
# println(qᵢ, qₙ, qᵢ .- qₙ)
# println(pᵢ, pₙ, pᵢ .- pₙ)
# println(q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 5E-6
@test pᵢ ≈ pₙ atol = 5E-8
@test q̇ᵢ ≈ q̇ₙ atol = 1E-4
@test ṗᵢ ≈ ṗₙ atol = 1E-6


# Midpoint Extrapolation for PODEs

history = (t=[t₀], q=[q₀], p=[p₀], q̇=[q̇₀], ṗ=[ṗ₀])

solutionstep!(sol, history, pode, MidpointExtrapolation(0))
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
# println(0, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(0, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-6
@test pᵢ ≈ pₙ atol = 1E-4
@test q̇ᵢ ≈ q̇ₙ atol = 1E-4
@test ṗᵢ ≈ ṗₙ atol = 1E-6

solutionstep!(sol, history, pode, MidpointExtrapolation(1))
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
# println(1, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(1, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-10
@test pᵢ ≈ pₙ atol = 1E-8
@test q̇ᵢ ≈ q̇ₙ atol = 1E-8
@test ṗᵢ ≈ ṗₙ atol = 1E-10

solutionstep!(sol, history, pode, MidpointExtrapolation(2))
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
# println(2, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(2, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-14
@test pᵢ ≈ pₙ atol = 1E-12
@test q̇ᵢ ≈ q̇ₙ atol = 1E-12
@test ṗᵢ ≈ ṗₙ atol = 1E-14

solutionstep!(sol, history, pode, MidpointExtrapolation(3))
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
# println(3, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(3, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-16
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-16

solutionstep!(sol, history, pode, MidpointExtrapolation(4))
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
# println(4, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(4, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-16
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-16

solutionstep!(sol, history, pode, MidpointExtrapolation(5))
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
# println(5, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(5, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-15
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-15


# Create IODE Solution Arrays

q₀ = initial_conditions(iode).q
p₀ = initial_conditions(iode).p

qᵢ = zero(q₀)
qₚ = zero(q₀)
qₙ = zero(q₀)

pᵢ = zero(p₀)
pₚ = zero(p₀)
pₙ = zero(p₀)

q̇ₚ = zero(q₀)
q̇₀ = zero(q₀)
q̇ₙ = zero(q₀)
q̇ᵢ = zero(q₀)

ṗₚ = zero(p₀)
ṗ₀ = zero(p₀)
ṗₙ = zero(p₀)
ṗᵢ = zero(p₀)

sol = (t=tᵢ, q=qᵢ, p=pᵢ, q̇=q̇ᵢ, ṗ=ṗᵢ)


# Compute Reference Solution for IODEs

qₚ .= [xₚ[1]]
pₚ .= [xₚ[2]]

qₙ .= [xₙ[1]]
pₙ .= [xₙ[2]]

initialguess(iode).v(q̇ₚ, tₚ, qₚ, pₚ, parameters(iode))
initialguess(iode).v(q̇₀, t₀, q₀, p₀, parameters(iode))
initialguess(iode).v(q̇ₙ, tₙ, qₙ, pₙ, parameters(iode))

functions(iode).ϑ(pₚ, tₚ, qₚ, q̇ₚ, parameters(iode))
functions(iode).ϑ(pₙ, tₙ, qₙ, q̇ₙ, parameters(iode))

initialguess(iode).f(ṗₚ, tₚ, qₚ, q̇ₚ, parameters(iode))
initialguess(iode).f(ṗ₀, t₀, q₀, q̇₀, parameters(iode))
initialguess(iode).f(ṗₙ, tₙ, qₙ, q̇ₙ, parameters(iode))


# Hermite Extrapolation

history = (t=[t₀, tₚ], q=[q₀, qₚ], p=[p₀, pₚ], q̇=[q̇₀, q̇ₚ], ṗ=[ṗ₀, ṗₚ])

solutionstep!(sol, history, iode, HermiteExtrapolation())
# println(qᵢ, qₙ, qᵢ .- qₙ)
# println(pᵢ, pₙ, pᵢ .- pₙ)
# println(q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 5E-6
@test pᵢ ≈ pₙ atol = 5E-8
@test q̇ᵢ ≈ q̇ₙ atol = 1E-4
@test ṗᵢ ≈ ṗₙ atol = 1E-6


# Midpoint Extrapolation for IODEs

history = (t=[t₀], q=[q₀], p=[p₀], q̇=[q̇₀], ṗ=[ṗ₀])

solutionstep!(sol, history, iode, MidpointExtrapolation(0))
# println(0, qᵢ, qₙ, qᵢ .- qₙ)
# println(0, pᵢ, pₙ, pᵢ .- pₙ)
# println(0, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(0, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-6
@test pᵢ ≈ pₙ atol = 1E-4
@test q̇ᵢ ≈ q̇ₙ atol = 1E-4
@test ṗᵢ ≈ ṗₙ atol = 1E-6

solutionstep!(sol, history, iode, MidpointExtrapolation(1))
# println(1, qᵢ, qₙ, qᵢ .- qₙ)
# println(1, pᵢ, pₙ, pᵢ .- pₙ)
# println(1, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(1, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-10
@test pᵢ ≈ pₙ atol = 1E-8
@test q̇ᵢ ≈ q̇ₙ atol = 1E-8
@test ṗᵢ ≈ ṗₙ atol = 1E-10

solutionstep!(sol, history, iode, MidpointExtrapolation(2))
# println(2, qᵢ, qₙ, qᵢ .- qₙ)
# println(2, pᵢ, pₙ, pᵢ .- pₙ)
# println(2, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(2, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-14
@test pᵢ ≈ pₙ atol = 1E-12
@test q̇ᵢ ≈ q̇ₙ atol = 1E-12
@test ṗᵢ ≈ ṗₙ atol = 1E-14

solutionstep!(sol, history, iode, MidpointExtrapolation(3))
# println(3, qᵢ, qₙ, qᵢ .- qₙ)
# println(3, pᵢ, pₙ, pᵢ .- pₙ)
# println(3, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(3, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-16
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-16

solutionstep!(sol, history, iode, MidpointExtrapolation(4))
# println(4, qᵢ, qₙ, qᵢ .- qₙ)
# println(4, pᵢ, pₙ, pᵢ .- pₙ)
# println(4, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(4, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-16
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-16

solutionstep!(sol, history, iode, MidpointExtrapolation(5))
# println(5, qᵢ, qₙ, qᵢ .- qₙ)
# println(5, pᵢ, pₙ, pᵢ .- pₙ)
# println(5, q̇ᵢ, q̇ₙ, q̇ᵢ .- q̇ₙ)
# println(5, ṗᵢ, ṗₙ, ṗᵢ .- ṗₙ)
@test qᵢ ≈ qₙ atol = 1E-15
@test pᵢ ≈ pₙ atol = 1E-16
@test q̇ᵢ ≈ q̇ₙ atol = 1E-16
@test ṗᵢ ≈ ṗₙ atol = 1E-15
