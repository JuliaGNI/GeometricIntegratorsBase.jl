using GeometricIntegratorsBase
using GeometricEquations
using Test

using GeometricSolutions: relative_maximum_error

using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonautonomousProblems: nonautonomous_ode_v

ode = odeproblem()
ref = exact_solution(ode)

# reference implementation of the scheme, with the stage time stated here rather than taken
# from the code under test. the times are read off the solution so that any disagreement is
# due to the scheme alone
function reference_explicit_euler(prob, sol)
    q = copy(initial_conditions(prob).q)
    v = zero(q)
    h = timestep(prob)

    for n in 1:lastindex(sol.t)
        nonautonomous_ode_v(v, sol.t[n-1], q, parameters(prob))
        q .+= h .* v
    end

    q
end

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

    # the explicit Euler method evaluates the vector field at the beginning of the step,
    # v(tₙ, qₙ). the harmonic oscillator is autonomous and cannot tell that apart from
    # v(tₙ₊₁, qₙ), so the stage time is pinned down on a non-autonomous problem
    nonautonomous = nonautonomous_odeproblem()
    sol = integrate(nonautonomous, ExplicitEuler())
    @test sol.q[end] == reference_explicit_euler(nonautonomous, sol)

    # the implicit method evaluates at the end of the step instead, and so must differ
    @test integrate(nonautonomous, ImplicitEuler()).q[end] != sol.q[end]

    @testset "Unsupported Problem Types" begin
        # neither method enforces a constraint or knows about the additional equations of a
        # differential algebraic problem, so those must be rejected rather than integrated as if
        # the constraint were absent. `DAEProblem` in particular is a member of the
        # `AbstractProblemODE` union and used to be integrated silently, which is the failure
        # mode the unconstrained unions in `src/GeometricIntegratorsBase.jl` exist to prevent
        for method in (ExplicitEuler(), ImplicitEuler())
            for prob in (daeproblem(), pdaeproblem(), hdaeproblem())
                err = try
                    integrate(prob, method)
                    nothing
                catch e
                    e
                end

                @test err isa ArgumentError
                @test occursin(string(nameof(typeof(method))), err.msg)
                @test occursin(string(nameof(typeof(equation(prob)))), err.msg)
            end
        end
    end

end
