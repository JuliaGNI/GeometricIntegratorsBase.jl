using Profile
using GeometricIntegratorsBase
using GeometricIntegratorsBase: Solution

include("harmonic_oscillator.jl")

using ..HarmonicOscillator

const Δt = 0.1
const nt = 10_000
const t₀ = 0.0
const t₁ = nt * Δt


function profile(method)
    ode = odeproblem(timespan=(t₀, t₁), timestep=Δt)
    int = GeometricIntegrator(ode, method)
    sol = Solution(ode)

    @time integrate!(sol, int)
    @time integrate!(sol, int)

    Profile.clear()
    Profile.clear_malloc_data()

    Profile.Allocs.@profile integrate!(sol, int)
end

profile(ExplicitEuler())
profile(ImplicitEuler())
