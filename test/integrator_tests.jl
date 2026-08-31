using GeometricIntegratorsBase
using Test

using GeometricIntegratorsBase: ODEMethod
using GeometricIntegratorsBase: equations, method, timestep
using GeometricEquations: GeometricProblem, ODEProblem
using SimpleSolvers: NonlinearSolverException

import GeometricIntegratorsBase: integrate_step!

import ..HarmonicOscillator: odeproblem

struct ExplicitEulerTest <: ODEMethod end

function integrate_step!(sol, hist, params, int::GeometricIntegrator{
        <:ExplicitEulerTest, <:GeometricProblem})
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
# The handling is tested directly first, on a method that throws at a chosen time. That pins down
# which step fails, what the surviving steps hold, what the rest of the solution holds, and that
# any other exception is rethrown — none of which then depend on how SimpleSolvers happens to
# arrive at the non-finite direction it rejects. The end-to-end test below is the complement: it
# checks that a real solve does reach that rejection.

struct FailingTest{ET} <: ODEMethod
    tfail::Float64
    exception::ET
end

function integrate_step!(sol, hist, params, int::GeometricIntegrator{
        <:FailingTest, <:GeometricProblem})
    # `reset!` sets `sol.t` to the time at the end of the step, so with a timestep of 0.1 and
    # `tfail = 0.35` the first step to throw is n = 4
    sol.t ≥ method(int).tfail && throw(method(int).exception)

    # a doubling per step: the steps that survive a breakdown are then exactly checkable, and
    # distinguishable from the zeros `Solution` allocates for the steps that never ran
    sol.q .*= 2

    return (
        solution = sol,
    )
end

failing = odeproblem(timespan = (0.0, 1.0), timestep = 0.1)

sol_failing = @test_logs (:warn, r"Nonlinear solver failed at timestep n=4") match_mode = :any integrate(
    failing, FailingTest(0.35, NonlinearSolverException("test breakdown")))

# the steps before the breakdown hold what they computed, and the exception did not propagate
@test all(n -> sol_failing.q[n] == 2^n .* sol_failing.q[0], 0:3)

# the steps from the breakdown on were never computed, so they are still the allocated zeros —
# which is what the warning has to say, since nothing in the solution itself marks the end
@test all(n -> all(iszero, sol_failing.q[n]), 4:lastindex(sol_failing.q))

# any other exception is a bug and must not be masked
@test_throws ErrorException integrate(failing, FailingTest(0.35, ErrorException("not a solver failure")))

# The complement: a real `Newton()` solve that does reach the rejection above. The vector field
# goes non-finite past a fixed time rather than blowing up gradually, so the step at which the
# solve breaks down is pinned down (with a timestep of 0.1 and the stage time of the implicit
# method at the end of the step, that is n = 4) instead of depending on how an overflow rounds.
function breakdown_v(v, t, q, params)
    v .= t ≥ 0.35 ? Inf : q
    nothing
end

breakdown = ODEProblem(breakdown_v, (0.0, 1.0), 0.1, [1.0, 0.0])

sol_breakdown = @test_logs (:warn, r"Nonlinear solver failed at timestep n=4") match_mode = :any integrate(
    breakdown, ImplicitEuler())

# the steps before the breakdown survive, and the exception did not propagate. asserting the exact
# values rather than `isfinite`, because the zeros left behind by an early `break` are finite too:
# implicit Euler on q̇ = q gives qₙ₊₁ = qₙ / (1 - h)
@test all(n -> sol_breakdown.q[n] ≈ [1 / (1 - 0.1)^n, 0.0], 0:3)
