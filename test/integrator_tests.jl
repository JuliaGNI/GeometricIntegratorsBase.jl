using GeometricIntegratorsBase
using Test

using GeometricIntegratorsBase: ODEMethod
using GeometricIntegratorsBase: equations, timestep
using GeometricEquations: GeometricProblem

import GeometricIntegratorsBase: integrate_step!

import ..HarmonicOscillator: odeproblem


struct ExplicitEulerTest <: ODEMethod end

function integrate_step!(sol, hist, params, int::GeometricIntegrator{<:ExplicitEulerTest, <:GeometricProblem})
    # compute vector field
    equations(int).v(sol.q̇, sol.t, sol.q, params)

    # compute update
    sol.q .+= timestep(int) .* sol.q̇

    return (
        solution = sol,
    )
end

sol = integrate(odeproblem(), ExplicitEulerTest())
