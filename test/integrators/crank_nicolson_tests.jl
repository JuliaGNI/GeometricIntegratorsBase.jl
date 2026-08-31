using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: CrankNicolsonCache, CrankNicolsonPODECache,
                                CrankNicolsonIODECache
using GeometricIntegratorsBase: CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonlinearProblems

# Accuracy, convergence order, data types, the agreement of the formulations of one and the same
# equation and the rejection of unsupported problem types are asserted for every method of the
# package in `common_tests.jl`.
@testset "$(rpad("CrankNicolson Method Tests", 80))" begin
    @testset "Method Properties" begin
        method = CrankNicolson()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test issymmetric(method)
        # the trapezoidal rule is conjugate to a symplectic method, but not symplectic itself
        @test !issymplectic(method)

        @test isodemethod(method)
        @test ispodemethod(method)
        @test ishodemethod(method)
        @test isiodemethod(method)
        @test islodemethod(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        method = CrankNicolson()

        @testset "ODE Problems" begin
            ode = odeproblem()

            cache = Cache{Float64}(ode, method)
            @test cache isa CrankNicolsonCache{Float64}

            @test axes(cache.x) == axes(initial_conditions(ode).q)
            @test axes(cache.q) == axes(initial_conditions(ode).q)
            @test axes(cache.v) == axes(initial_conditions(ode).q)
            @test axes(cache.v̄) == axes(initial_conditions(ode).q)

            @test nlsolution(cache) === cache.x

            @test solversize(method, ode) == length(initial_conditions(ode).q)

            @test CacheType(Float64, ode, method) == CrankNicolsonCache{Float64}
            @test CacheType(Float32, ode, method) == CrankNicolsonCache{Float32}
        end

        @testset "PODE/HODE Problems" begin
            pode = podeproblem()

            cache = Cache{Float64}(pode, method)
            @test cache isa CrankNicolsonPODECache{Float64}

            # both stage vector fields are solved for, so the solution vector holds two of them
            @test length(cache.x) ==
                  length(initial_conditions(pode).q) + length(initial_conditions(pode).p)
            @test solversize(method, pode) == length(cache.x)

            @test axes(cache.q) == axes(initial_conditions(pode).q)
            @test axes(cache.p) == axes(initial_conditions(pode).p)
            @test axes(cache.v) == axes(initial_conditions(pode).q)
            @test axes(cache.f) == axes(initial_conditions(pode).p)

            @test axes(cache.v̄) == axes(initial_conditions(pode).q)
            @test axes(cache.f̄) == axes(initial_conditions(pode).p)

            @test nlsolution(cache) === cache.x

            @test CacheType(Float64, pode, method) == CrankNicolsonPODECache{Float64}
            @test CacheType(Float32, pode, method) == CrankNicolsonPODECache{Float32}

            # the Hamiltonian formulation is integrated with the very same cache
            @test Cache{Float64}(hodeproblem(), method) isa CrankNicolsonPODECache{Float64}
        end

        @testset "IODE/LODE Problems" begin
            iode = iodeproblem()

            cache = Cache{Float64}(iode, method)
            @test cache isa CrankNicolsonIODECache{Float64}

            # the velocity at the beginning of the time step is not given by a function
            # evaluation but implicitly by ϑ(t̄, q̄, v̄) = p̄, so it is solved for alongside the
            # velocity at the end of the time step and the solution vector is twice as long
            @test length(cache.x) == 2 * length(initial_conditions(iode).q)
            @test solversize(method, iode) == 2 * length(initial_conditions(iode).q)

            @test axes(cache.q) == axes(initial_conditions(iode).q)
            @test axes(cache.v) == axes(initial_conditions(iode).q)
            @test axes(cache.θ) == axes(initial_conditions(iode).q)
            @test axes(cache.f) == axes(initial_conditions(iode).q)
            @test axes(cache.v̄) == axes(initial_conditions(iode).q)
            @test axes(cache.θ̄) == axes(initial_conditions(iode).q)
            @test axes(cache.f̄) == axes(initial_conditions(iode).q)

            @test nlsolution(cache) === cache.x

            @test CacheType(Float64, iode, method) == CrankNicolsonIODECache{Float64}
            @test CacheType(Float32, iode, method) == CrankNicolsonIODECache{Float32}

            # the Lagrangian formulation is integrated with the very same cache
            @test Cache{Float64}(lodeproblem(), method) isa CrankNicolsonIODECache{Float64}
        end
    end

    @testset "Energy Behaviour" begin
        # the harmonic oscillator energy is quadratic and the method coincides with the midpoint
        # rule there, so it is preserved up to round-off in every formulation
        for problem in (odeproblem, podeproblem, hodeproblem, iodeproblem, lodeproblem)
            prob = problem(; timespan = (0.0, 100.0))
            @test max_energy_error(integrate(prob, CrankNicolson()), prob) < 1E-13
        end
    end

    @testset "Accuracy against the Euler Methods" begin
        # second order, so considerably more accurate than the first order methods
        ode = odeproblem()
        ref = exact_solution(ode)

        err = relative_maximum_error(integrate(ode, CrankNicolson()), ref).q
        @test err < relative_maximum_error(integrate(ode, ImplicitEuler()), ref).q / 10
    end

    @testset "Nonlinear Problem" begin
        # the harmonic oscillator is linear, where the trapezoidal rule and the midpoint rule are
        # the same map. the method is second order on a nonlinear problem too
        errs = [riccati_error(integrate(riccati_problem(; timestep = Δt), CrankNicolson()))
                for Δt in (0.1, 0.05, 0.025)]

        @test all(isapprox.(errs[1:(end - 1)] ./ errs[2:end], 4; atol = 2E-1))
        @test errs[1] < 1E-3
    end

    @testset "Comparison with ImplicitMidpoint" begin
        # for linear problems the trapezoidal rule and the midpoint rule are the same map, so the
        # two coincide on the harmonic oscillator in every formulation ...
        for problem in (odeproblem, podeproblem, hodeproblem, iodeproblem, lodeproblem)
            prob = problem()
            cn = integrate(prob, CrankNicolson())
            mp = integrate(prob, ImplicitMidpoint())
            @test maximum(maximum(abs.(cn.q[n] .- mp.q[n])) for n in axes(cn.q, 1)) < 1E-14
        end

        # ... but they are genuinely different methods on a nonlinear problem, and once the
        # coefficients of a linear one depend on time
        for prob in (riccati_problem(), nonlinear_podeproblem(),
            nonautonomous_podeproblem(), nonautonomous_iodeproblem())
            @test integrate(prob, CrankNicolson()).q !=
                  integrate(prob, ImplicitMidpoint()).q
        end
    end
end
