using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitEulerCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator


# Accuracy, convergence order, data types and the rejection of unsupported problem types are
# asserted for every method of the package in `common_tests.jl`.
@testset "$(rpad("ImplicitEuler Method Tests", 80))" begin

    @testset "Method Properties" begin
        method = ImplicitEuler()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test !issymmetric(method)
        @test !issymplectic(method)

        @test isodemethod(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        ode = odeproblem()
        method = ImplicitEuler()

        cache = Cache{Float64}(ode, method)
        @test cache isa ImplicitEulerCache{Float64}

        @test axes(cache.x) == axes(initial_conditions(ode).q)
        @test axes(cache.q) == axes(initial_conditions(ode).q)
        @test axes(cache.v) == axes(initial_conditions(ode).q)
        @test axes(cache.v̄) == axes(initial_conditions(ode).q)

        @test nlsolution(cache) === cache.x

        @test solversize(method, ode) == length(initial_conditions(ode).q)

        @test CacheType(Float64, ode, method) == ImplicitEulerCache{Float64}
        @test CacheType(Float32, ode, method) == ImplicitEulerCache{Float32}
    end

    @testset "Different Timesteps" begin
        # the convergence order test in `common_tests.jl` asserts this asymptotically; here it is
        # asserted at the timesteps a caller would actually reach for
        errs = [
            begin
                ode = odeproblem([0.5, 0.0]; timestep=Δt)
                relative_maximum_error(integrate(ode, ImplicitEuler()), exact_solution(ode)).q
            end for Δt in (0.05, 0.2)
        ]

        @test errs[1] < 2E-2
        @test errs[2] < 1E-1
        @test errs[1] < errs[2]
    end

    @testset "Different Initial Conditions" begin
        # note that the reference solution is singular for q₀ = 0
        for x₀ in ([1.0, 0.5], [0.2, 0.8])
            ode = odeproblem(x₀)
            @test relative_maximum_error(integrate(ode, ImplicitEuler()), exact_solution(ode)).q < 5E-2
        end
    end

    @testset "Solver Options" begin
        ode = odeproblem()

        # loose tolerances still converge, and `max_iterations` is accepted as a solver option
        for options in ((x_abstol=1e-4, f_abstol=1e-4), (max_iterations=200,))
            sol = integrate(ode, ImplicitEuler(); options...)
            @test all(x -> all(isfinite, x), sol.q)
        end
    end

    @testset "Energy Behaviour" begin
        # the method is not symplectic, so the energy of the harmonic oscillator is not preserved,
        # but neither does it run away over the default time span
        ode = odeproblem()
        @test max_energy_error(integrate(ode, ImplicitEuler()), ode) < 0.1
    end

    @testset "Comparison with ExplicitEuler" begin
        ode = odeproblem()
        ref = exact_solution(ode)

        err_implicit = relative_maximum_error(integrate(ode, ImplicitEuler()), ref).q
        err_explicit = relative_maximum_error(integrate(ode, ExplicitEuler()), ref).q

        # both are first order methods, so their errors are of the same magnitude on this problem
        @test min(err_implicit, err_explicit) / max(err_implicit, err_explicit) > 0.2
    end

    @testset "Stiff Problem Behaviour" begin
        # the point of an implicit method: q̇ = -100 q at Δt = 0.01 is far outside the stability
        # region of the explicit method, whereas implicit Euler gives qₙ₊₁ = qₙ / (1 + 100 Δt)
        # and so decays monotonically at any timestep
        stiff_v(v, t, q, params) = (v[1] = -100.0 * q[1]; nothing)
        sol = integrate(ODEProblem(stiff_v, (0.0, 0.1), 0.01, [1.0]), ImplicitEuler())

        @test all(n -> sol.q[n][1] ≈ 1 / 2^n, axes(sol.q, 1))
        @test issorted([sol.q[n][1] for n in axes(sol.q, 1)]; rev=true)
    end

    @testset "Edge Cases" begin
        # a time span resolved far more finely than it needs to be, and one of a single step
        for (timestep, timespan) in ((1e-4, (0.0, 0.01)), (0.1, (0.0, 0.1)))
            sol = integrate(odeproblem([0.5, 0.0]; timestep, timespan), ImplicitEuler())

            @test length(sol.t) == round(Int, (timespan[2] - timespan[1]) / timestep) + 1
            @test all(x -> all(isfinite, x), sol.q)
        end
    end

    @testset "Method Interface" begin
        ode = odeproblem()
        int = GeometricIntegrator(ode, ImplicitEuler())

        @test int isa GeometricIntegrator
        @test typeof(int).parameters[1] == ImplicitEuler
        @test typeof(int).parameters[2] == typeof(ode)
    end

end
