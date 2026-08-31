using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitMidpointCache, ImplicitMidpointPODECache,
                                ImplicitMidpointIODECache
using GeometricIntegratorsBase: CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator
using ..NonlinearProblems

# Accuracy, convergence order, data types, the agreement of the formulations of one and the same
# equation and the rejection of unsupported problem types are asserted for every method of the
# package in `common_tests.jl`.
@testset "$(rpad("ImplicitMidpoint Method Tests", 80))" begin
    @testset "Method Properties" begin
        method = ImplicitMidpoint()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test issymmetric(method)
        @test issymplectic(method)

        @test isodemethod(method)
        @test ispodemethod(method)
        @test ishodemethod(method)
        @test isiodemethod(method)
        @test islodemethod(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        method = ImplicitMidpoint()

        @testset "ODE Problems" begin
            ode = odeproblem()

            cache = Cache{Float64}(ode, method)
            @test cache isa ImplicitMidpointCache{Float64}

            @test axes(cache.x) == axes(initial_conditions(ode).q)
            @test axes(cache.q) == axes(initial_conditions(ode).q)
            @test axes(cache.v) == axes(initial_conditions(ode).q)

            @test nlsolution(cache) === cache.x

            @test solversize(method, ode) == length(initial_conditions(ode).q)

            @test CacheType(Float64, ode, method) == ImplicitMidpointCache{Float64}
            @test CacheType(Float32, ode, method) == ImplicitMidpointCache{Float32}
        end

        @testset "PODE/HODE Problems" begin
            pode = podeproblem()

            cache = Cache{Float64}(pode, method)
            @test cache isa ImplicitMidpointPODECache{Float64}

            # both stage vector fields are solved for, so the solution vector holds two of them
            @test length(cache.x) ==
                  length(initial_conditions(pode).q) + length(initial_conditions(pode).p)
            @test solversize(method, pode) == length(cache.x)

            @test axes(cache.q) == axes(initial_conditions(pode).q)
            @test axes(cache.p) == axes(initial_conditions(pode).p)
            @test axes(cache.v) == axes(initial_conditions(pode).q)
            @test axes(cache.f) == axes(initial_conditions(pode).p)

            @test nlsolution(cache) === cache.x

            @test CacheType(Float64, pode, method) == ImplicitMidpointPODECache{Float64}
            @test CacheType(Float32, pode, method) == ImplicitMidpointPODECache{Float32}

            # the Hamiltonian formulation is integrated with the very same cache
            @test Cache{Float64}(hodeproblem(), method) isa
                  ImplicitMidpointPODECache{Float64}
        end

        @testset "IODE/LODE Problems" begin
            iode = iodeproblem()

            cache = Cache{Float64}(iode, method)
            @test cache isa ImplicitMidpointIODECache{Float64}

            @test axes(cache.x) == axes(initial_conditions(iode).q)
            @test axes(cache.q) == axes(initial_conditions(iode).q)
            @test axes(cache.v) == axes(initial_conditions(iode).q)
            @test axes(cache.θ) == axes(initial_conditions(iode).q)
            @test axes(cache.f) == axes(initial_conditions(iode).q)

            @test nlsolution(cache) === cache.x

            # the solver variable is the stage velocity alone, so the nonlinear system has the
            # same size as for an ordinary differential equation
            @test solversize(method, iode) == length(initial_conditions(iode).q)

            @test CacheType(Float64, iode, method) == ImplicitMidpointIODECache{Float64}
            @test CacheType(Float32, iode, method) == ImplicitMidpointIODECache{Float32}

            # the Lagrangian formulation is integrated with the very same cache
            @test Cache{Float64}(lodeproblem(), method) isa
                  ImplicitMidpointIODECache{Float64}
        end
    end

    @testset "Energy Conservation" begin
        # the method conserves quadratic invariants exactly, so the energy of the harmonic
        # oscillator is preserved up to round-off, even over long times and in every formulation
        for problem in (odeproblem, podeproblem, hodeproblem, iodeproblem, lodeproblem)
            prob = problem(; timespan = (0.0, 100.0))
            @test max_energy_error(integrate(prob, ImplicitMidpoint()), prob) < 1E-13
        end
    end

    @testset "Different Initial Conditions" begin
        # note that the reference solution is singular for q₀ = 0
        for x₀ in ([1.0, 0.0], [0.2, 0.8], [0.5, 0.5])
            ode = odeproblem(x₀)
            @test relative_maximum_error(integrate(ode, ImplicitMidpoint()), exact_solution(ode)).q <
                  1E-3
        end
    end

    @testset "Accuracy against the Euler Methods" begin
        # second order, so considerably more accurate than the first order methods
        ode = odeproblem()
        ref = exact_solution(ode)

        err = relative_maximum_error(integrate(ode, ImplicitMidpoint()), ref).q
        @test err < relative_maximum_error(integrate(ode, ImplicitEuler()), ref).q / 10
    end

    @testset "Nonlinear Problem" begin
        # the midpoint rule and the trapezoidal rule are the same map on a linear problem, so
        # they are only distinguishable on a nonlinear one, where both remain second order
        errs = [riccati_error(integrate(riccati_problem(; timestep = Δt), ImplicitMidpoint()))
                for Δt in (0.1, 0.05, 0.025)]

        @test all(isapprox.(errs[1:(end - 1)] ./ errs[2:end], 4; atol = 2E-1))
        @test errs[1] < 1E-3
    end
end
