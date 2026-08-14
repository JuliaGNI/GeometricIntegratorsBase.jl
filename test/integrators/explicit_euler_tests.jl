using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricIntegratorsBase: ExplicitEulerCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonautonomousProblems: nonautonomous_ode_v


# Reference implementation of the scheme, written out directly so that the stage time is stated
# here rather than taken from the code under test. The times are read off the solution rather
# than recomputed, so that any disagreement is due to the scheme alone.
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


# Accuracy, convergence order, data types and the rejection of unsupported problem types are
# asserted for every method of the package in `common_tests.jl`.
@testset "$(rpad("ExplicitEuler Method Tests", 80))" begin

    @testset "Method Properties" begin
        method = ExplicitEuler()

        @test isexplicit(method)
        @test !isimplicit(method)
        @test !issymmetric(method)
        @test !issymplectic(method)

        @test isodemethod(method)

        # explicit methods need neither a solver nor an initial guess
        @test default_solver(method) isa NoSolver
        @test default_iguess(method) isa NoInitialGuess
        @test solversize(method, odeproblem()) == 0
    end

    @testset "Cache Structure" begin
        ode = odeproblem()
        method = ExplicitEuler()

        cache = Cache{Float64}(ode, method)
        @test cache isa ExplicitEulerCache{Float64}

        @test axes(cache.v) == axes(initial_conditions(ode).q)

        # explicit methods have no nonlinear solver solution
        @test ismissing(nlsolution(cache))

        @test CacheType(Float64, ode, method) == ExplicitEulerCache{Float64}
        @test CacheType(Float32, ode, method) == ExplicitEulerCache{Float32}
    end

    @testset "Non-Autonomous Stage Time" begin
        # the method evaluates the vector field at the beginning of the step, v(tₙ, qₙ). the
        # harmonic oscillator is autonomous and cannot tell that apart from v(tₙ₊₁, qₙ), so the
        # stage time is pinned down against the reference implementation above
        prob = nonautonomous_odeproblem()
        sol = integrate(prob, ExplicitEuler())

        @test sol.q[end] == reference_explicit_euler(prob, sol)

        # the implicit method evaluates at the end of the step instead, and so must differ
        @test integrate(prob, ImplicitEuler()).q[end] != sol.q[end]
    end

end
