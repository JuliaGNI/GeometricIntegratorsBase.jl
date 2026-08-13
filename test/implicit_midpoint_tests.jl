using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitMidpointCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator


# maximum relative energy error along an ODE solution
function max_energy_error(sol, prob)
    maximum(abs, compute_energy_error(sol.t, sol.q, parameters(prob))[2])
end


@testset "$(rpad("ImplicitMidpoint Method Tests", 80))" begin

    @testset "Method Properties" begin
        method = ImplicitMidpoint()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test issymmetric(method)
        @test issymplectic(method)

        @test isodemethod(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        ode = odeproblem()
        method = ImplicitMidpoint()

        cache = Cache{Float64}(ode, method)
        @test cache isa ImplicitMidpointCache{Float64}

        @test axes(cache.x) == axes(initial_conditions(ode).q)
        @test axes(cache.q) == axes(initial_conditions(ode).q)
        @test axes(cache.v) == axes(initial_conditions(ode).q)
        @test axes(cache.v̄) == axes(initial_conditions(ode).q)

        @test nlsolution(cache) === cache.x

        @test solversize(method, ode) == length(initial_conditions(ode).q)

        @test CacheType(Float64, ode, method) == ImplicitMidpointCache{Float64}
        @test CacheType(Float32, ode, method) == ImplicitMidpointCache{Float32}
    end

    @testset "Integration Accuracy" begin
        ode = odeproblem()
        ref = exact_solution(ode)

        sol = integrate(ode, ImplicitMidpoint())
        err = relative_maximum_error(sol, ref)
        @test err.q < 1E-3

        # a converged solve is silent, so tight tolerances must not produce solver warnings
        sol_tight = @test_nowarn integrate(ode, ImplicitMidpoint();
            x_abstol=1e-12,
            f_abstol=1e-12,
        )
        @test relative_maximum_error(sol_tight, ref).q < 1E-3

        # considerably more accurate than the first order Euler methods
        @test err.q < relative_maximum_error(integrate(ode, ImplicitEuler()), ref).q / 10
    end

    @testset "Convergence Order" begin
        # second order, so halving the timestep should reduce the error by a factor of four
        errs = [
            begin
                prob = odeproblem(; timestep=Δt)
                relative_maximum_error(integrate(prob, ImplicitMidpoint()), exact_solution(prob)).q
            end for Δt in (0.1, 0.05, 0.025)
        ]

        @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=2E-1))
    end

    @testset "Energy Conservation" begin
        # the implicit midpoint method conserves quadratic invariants exactly, so the energy
        # of the harmonic oscillator is preserved up to round-off, even over long times
        ode = odeproblem(; timespan=(0.0, 100.0))
        @test max_energy_error(integrate(ode, ImplicitMidpoint()), ode) < 1E-13
    end

    @testset "Different Initial Conditions" begin
        # note that the reference solution is singular for q₀ = 0
        for x₀ in ([1.0, 0.0], [0.2, 0.8], [0.5, 0.5])
            ode = odeproblem(x₀)
            sol = integrate(ode, ImplicitMidpoint())
            @test relative_maximum_error(sol, exact_solution(ode)).q < 1E-3
        end
    end

    @testset "Data Type Consistency" begin
        ode = odeproblem()
        sol = integrate(ode, ImplicitMidpoint())

        @test eltype(sol.q[0]) == Float64
        @test all(isfinite, sol.q[end])
    end

end
