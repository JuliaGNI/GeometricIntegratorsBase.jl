using GeometricIntegratorsBase
using Test

using GeometricIntegratorsBase: functions, timestep, value
using GeometricIntegratorsBase: extrapolate!, initialguess, initialstate, initialtime
using GeometricIntegratorsBase: StateVariable, VectorfieldVariable

using ..HarmonicOscillator

ode = odeproblem()
pode = podeproblem()
iode = iodeproblem()


# Compute Reference Solution for ODEs

const Δt = timestep(ode)
const t₀ = initialtime(ode)
const t₁ = t₀ + Δt
const t₂ = t₁ + Δt
const t₋ = t₀ - Δt
const tₚ = t₋
const tₙ = t₁
const tᵢ = tₙ

x₀ = initialstate(ode).q

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

ẋ₀ = VectorfieldVariable(x₀)
ẋ₁ = VectorfieldVariable(x₁)
ẋ₂ = VectorfieldVariable(x₂)
ẋᵢ = VectorfieldVariable(xᵢ)
ẋₙ = VectorfieldVariable(xₙ)
ẋₚ = VectorfieldVariable(xₚ)

functions(ode).v(ẋₚ, tₚ, xₚ, parameters(ode))
functions(ode).v(ẋ₀, t₀, x₀, parameters(ode))
functions(ode).v(ẋₙ, tₙ, xₙ, parameters(ode))


# Create SolutionStep for ODE Tests
sol = SolutionStep(ode; nhistory=2)

copy!(sol, tₚ, (q=xₚ, q̇=ẋₚ))
reset!(sol, Δt)

copy!(sol, t₀, (q=x₀, q̇=ẋ₀))
reset!(sol, Δt)


# Hermite Extrapolation

extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, tᵢ, xᵢ, ẋᵢ, HermiteExtrapolation())

# println(xᵢ, xₙ, xᵢ .- xₙ)
# println(ẋᵢ, ẋₙ, ẋᵢ .- ẋₙ)

@test xᵢ ≈ xₙ atol = 1E-5
@test ẋᵢ ≈ ẋₙ atol = 1E-4

@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, HermiteExtrapolation()) == xᵢ
@test extrapolate!(tₚ, xₚ, ẋₚ, t₀, x₀, ẋ₀, t₁, x₁, ẋ₁, HermiteExtrapolation()) == (xᵢ, ẋᵢ)


# Hermite Extrapolation for ODE solutionstep
copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, HermiteExtrapolation())
@test sol.t == tᵢ
@test sol.q == xᵢ
@test sol.q̇ == ẋᵢ


# Euler Extrapolation for ODEs

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(0))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-2
@test sol.q̇ ≈ ẋₙ atol = 5E-2

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(1))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-3
@test sol.q̇ ≈ ẋₙ atol = 5E-3

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(2))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-5
@test sol.q̇ ≈ ẋₙ atol = 5E-5

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(3))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-6
@test sol.q̇ ≈ ẋₙ atol = 1E-6

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(4))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-8
@test sol.q̇ ≈ ẋₙ atol = 1E-8

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, EulerExtrapolation(5))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-10
@test sol.q̇ ≈ ẋₙ atol = 1E-10


# Midpoint Extrapolation for ODEs

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(0))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 5E-5
@test sol.q̇ ≈ ẋₙ atol = 5E-5

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(1))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-8
@test sol.q̇ ≈ ẋₙ atol = 1E-8

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(2))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-12
@test sol.q̇ ≈ ẋₙ atol = 1E-12

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(3))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-16
@test sol.q̇ ≈ ẋₙ atol = 1E-16

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(4))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-16
@test sol.q̇ ≈ ẋₙ atol = 1E-16

copy!(sol, State(t₁, (q=x₀, q̇=ẋ₀)))
solutionstep!(current(sol), state(sol), ode, MidpointExtrapolation(5))
# println(sol.q, xₙ, sol.q .- xₙ)
# println(sol.q̇, xₙ, sol.q̇ .- ẋₙ)
@test sol.q ≈ xₙ atol = 1E-15
@test sol.q̇ ≈ ẋₙ atol = 1E-15


# Create PODE Solution Arrays

q₀ = initialstate(pode).q
p₀ = initialstate(pode).p

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


# Create SolutionStep for PODE Tests

sol = SolutionStep(pode; nhistory=2)
copy!(sol, tₚ, (q=qₚ, p=pₚ, q̇=q̇ₚ, ṗ=ṗₚ))
reset!(sol, Δt)

copy!(sol, t₀, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
reset!(sol, Δt)


# Hermite Extrapolation

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, HermiteExtrapolation())
# println(sol.q, qₙ, sol.q .- qₙ)
# println(sol.p, pₙ, sol.p .- pₙ)
# println(sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 5E-6
@test sol.p ≈ pₙ atol = 5E-8
@test sol.q̇ ≈ q̇ₙ atol = 1E-4
@test sol.ṗ ≈ ṗₙ atol = 1E-6


# Midpoint Extrapolation for PODEs

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(0))
# println(0, sol.q, qₙ, sol.q .- qₙ)
# println(0, sol.p, pₙ, sol.p .- pₙ)
# println(0, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(0, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-6
@test sol.p ≈ pₙ atol = 1E-4
@test sol.q̇ ≈ q̇ₙ atol = 1E-4
@test sol.ṗ ≈ ṗₙ atol = 1E-6

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(1))
# println(1, sol.q, qₙ, sol.q .- qₙ)
# println(1, sol.p, pₙ, sol.p .- pₙ)
# println(1, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(1, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-10
@test sol.p ≈ pₙ atol = 1E-8
@test sol.q̇ ≈ q̇ₙ atol = 1E-8
@test sol.ṗ ≈ ṗₙ atol = 1E-10

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(2))
# println(2, sol.q, qₙ, sol.q .- qₙ)
# println(2, sol.p, pₙ, sol.p .- pₙ)
# println(2, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(2, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-14
@test sol.p ≈ pₙ atol = 1E-12
@test sol.q̇ ≈ q̇ₙ atol = 1E-12
@test sol.ṗ ≈ ṗₙ atol = 1E-14

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(3))
# println(3, sol.q, qₙ, sol.q .- qₙ)
# println(3, sol.p, pₙ, sol.p .- pₙ)
# println(3, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(3, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-16
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-16

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(4))
# println(4, sol.q, qₙ, sol.q .- qₙ)
# println(4, sol.p, pₙ, sol.p .- pₙ)
# println(4, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(4, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-16
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-16

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), pode, MidpointExtrapolation(5))
# println(5, sol.q, qₙ, sol.q .- qₙ)
# println(5, sol.p, pₙ, sol.p .- pₙ)
# println(5, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(5, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-15
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-15


# Create IODE Solution Arrays

q₀ = initialstate(iode).q
p₀ = initialstate(iode).p

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


# Create SolutionStep for IODE Tests

sol = SolutionStep(iode; nhistory=2)
copy!(sol, tₚ, (q=qₚ, p=pₚ, q̇=q̇ₚ, ṗ=ṗₚ))
reset!(sol, Δt)

copy!(sol, t₀, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
reset!(sol, Δt)


# Hermite Extrapolation

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, HermiteExtrapolation())
# println(sol.q, qₙ, sol.q .- qₙ)
# println(sol.p, pₙ, sol.p .- pₙ)
# println(sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 5E-6
@test sol.p ≈ pₙ atol = 5E-8
@test sol.q̇ ≈ q̇ₙ atol = 1E-4
@test sol.ṗ ≈ ṗₙ atol = 1E-6


# Midpoint Extrapolation for IODEs

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(0))
# println(0, sol.q, qₙ, sol.q .- qₙ)
# println(0, sol.p, pₙ, sol.p .- pₙ)
# println(0, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(0, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-6
@test sol.p ≈ pₙ atol = 1E-4
@test sol.q̇ ≈ q̇ₙ atol = 1E-4
@test sol.ṗ ≈ ṗₙ atol = 1E-6

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(1))
# println(1, sol.q, qₙ, sol.q .- qₙ)
# println(1, sol.p, pₙ, sol.p .- pₙ)
# println(1, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(1, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-10
@test sol.p ≈ pₙ atol = 1E-8
@test sol.q̇ ≈ q̇ₙ atol = 1E-8
@test sol.ṗ ≈ ṗₙ atol = 1E-10

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(2))
# println(2, sol.q, qₙ, sol.q .- qₙ)
# println(2, sol.p, pₙ, sol.p .- pₙ)
# println(2, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(2, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-14
@test sol.p ≈ pₙ atol = 1E-12
@test sol.q̇ ≈ q̇ₙ atol = 1E-12
@test sol.ṗ ≈ ṗₙ atol = 1E-14

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(3))
# println(3, sol.q, qₙ, sol.q .- qₙ)
# println(3, sol.p, pₙ, sol.p .- pₙ)
# println(3, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(3, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-16
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-16

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(4))
# println(4, sol.q, qₙ, sol.q .- qₙ)
# println(4, sol.p, pₙ, sol.p .- pₙ)
# println(4, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(4, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-16
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-16

copy!(sol, t₁, (q=q₀, p=p₀, q̇=q̇₀, ṗ=ṗ₀))
solutionstep!(current(sol), state(sol), iode, MidpointExtrapolation(5))
# println(5, sol.q, qₙ, sol.q .- qₙ)
# println(5, sol.p, pₙ, sol.p .- pₙ)
# println(5, sol.q̇, q̇ₙ, sol.q̇ .- q̇ₙ)
# println(5, sol.ṗ, ṗₙ, sol.ṗ .- ṗₙ)
@test sol.q ≈ qₙ atol = 1E-15
@test sol.p ≈ pₙ atol = 1E-16
@test sol.q̇ ≈ q̇ₙ atol = 1E-16
@test sol.ṗ ≈ ṗₙ atol = 1E-15
