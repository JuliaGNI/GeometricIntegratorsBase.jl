using GeometricIntegratorsBase
using Test

using GeometricIntegratorsBase: ODEMethod
using GeometricIntegratorsBase: equations, timestep
using GeometricEquations: GeometricProblem, ODEProblem

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


# A nonlinear solve that breaks down throws a `NonlinearSolverException` rather than returning a
# bad iterate. The time-stepping loop has to catch it, or `integrate` — which allocates the
# solution internally — would discard the whole trajectory instead of the part that is valid.
#
# The vector field goes non-finite past a fixed time rather than blowing up gradually, so the step
# at which the solve breaks down is pinned down (with a timestep of 0.1 and the stage time of the
# implicit method at the end of the step, that is n = 4) instead of depending on how an overflow
# happens to round.
function breakdown_v(v, t, q, params)
    v .= t ≥ 0.35 ? Inf : q
    nothing
end

breakdown = ODEProblem(breakdown_v, (0.0, 1.0), 0.1, [1.0, 0.0])

sol_breakdown = @test_logs (:warn, r"Nonlinear solver failed at timestep n=4") match_mode = :any integrate(breakdown, ImplicitEuler())

# the steps before the breakdown survive, and the exception did not propagate
@test all(n -> all(isfinite, sol_breakdown.q[n]), 0:3)
