using GeometricIntegratorsBase
using Test

using GeometricSolutions: relative_maximum_error

using ..HarmonicOscillator

ode = odeproblem()
ref = exact_solution(ode)

@testset "$(rpad("Euler integrators", 80))" begin

    sol = integrate(ode, ExplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

    sol = integrate(ode, ImplicitEuler())
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

    # ImplicitEuler with explicit (tight but reachable) solver tolerances.
    # min_iterations is retained via merge with default_options, so it need not be
    # restated here; a converged solve is now silent, so assert no solver warnings.
    sol = @test_nowarn integrate(ode, ImplicitEuler();
        x_abstol=1e-12,
        f_abstol=1e-12,
    )
    err = relative_maximum_error(sol, ref)
    @test err.q < 5E-2

end
