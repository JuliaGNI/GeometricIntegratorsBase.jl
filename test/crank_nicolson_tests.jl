using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: CrankNicolsonCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator
using ..NonautonomousProblems

# nonlinear test problem q̇ = -q² with exact solution q(t) = 1 / (1 + t)
riccati_v(v, t, q, params) = (v[1] = -q[1]^2; nothing)
riccati_problem(; timestep=0.1) = ODEProblem(riccati_v, (0.0, 1.0), timestep, [1.0])
riccati_error(sol) = maximum(abs(sol.q[n][1] - 1 / (1 + sol.t[n])) for n in axes(sol.q, 1))


@testset "$(rpad("CrankNicolson Method Tests", 80))" begin

    @testset "Method Properties" begin
        method = CrankNicolson()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test issymmetric(method)
        # the trapezoidal rule is conjugate to a symplectic method, but not symplectic itself
        @test !issymplectic(method)

        @test isodemethod(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        ode = odeproblem()
        method = CrankNicolson()

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

    @testset "Integration Accuracy" begin
        ode = odeproblem()
        ref = exact_solution(ode)

        sol = integrate(ode, CrankNicolson())
        err = relative_maximum_error(sol, ref)
        @test err.q < 1E-3

        # a converged solve is silent, so tight tolerances must not produce solver warnings
        sol_tight = @test_nowarn integrate(ode, CrankNicolson();
            x_abstol=1e-12,
            f_abstol=1e-12,
        )
        @test relative_maximum_error(sol_tight, ref).q < 1E-3

        @test err.q < relative_maximum_error(integrate(ode, ImplicitEuler()), ref).q / 10
    end

    @testset "Convergence Order" begin
        # second order on the linear problem
        errs = [
            begin
                prob = odeproblem(; timestep=Δt)
                relative_maximum_error(integrate(prob, CrankNicolson()), exact_solution(prob)).q
            end for Δt in (0.1, 0.05, 0.025)
        ]
        @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=2E-1))

        # and on a nonlinear problem
        errs = [riccati_error(integrate(riccati_problem(; timestep=Δt), CrankNicolson()))
                for Δt in (0.1, 0.05, 0.025)]
        @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=2E-1))
    end

    @testset "Non-Autonomous Convergence Order" begin
        # both the harmonic oscillator and the Riccati problem above are autonomous, so
        # neither detects a wrong stage time. taking v̄ or v at the wrong end of the interval
        # costs an order of accuracy, which this problem does detect
        errs = [
            begin
                prob = nonautonomous_odeproblem(; timestep=Δt)
                relative_maximum_error(integrate(prob, CrankNicolson()),
                    nonautonomous_ode_solution(prob)).q
            end for Δt in (0.1, 0.05, 0.025)
        ]

        @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=3E-1))
    end

    @testset "Comparison with ImplicitMidpoint" begin
        # for linear problems the trapezoidal rule and the midpoint rule are the same map,
        # so both give identical results on the harmonic oscillator ...
        ode = odeproblem()
        cn = integrate(ode, CrankNicolson())
        mp = integrate(ode, ImplicitMidpoint())
        @test maximum(maximum(abs.(cn.q[n] .- mp.q[n])) for n in axes(cn.q, 1)) < 1E-14

        # ... but they are genuinely different methods on a nonlinear problem
        riccati = riccati_problem()
        sol_cn = integrate(riccati, CrankNicolson())
        sol_mp = integrate(riccati, ImplicitMidpoint())
        @test sol_cn.q != sol_mp.q
        @test riccati_error(sol_cn) < 1E-3
        @test riccati_error(sol_mp) < 1E-3
    end

    @testset "Energy Behaviour" begin
        # the harmonic oscillator energy is quadratic and the method coincides with the
        # midpoint rule there, so it is preserved up to round-off
        ode = odeproblem(; timespan=(0.0, 100.0))
        @test max_energy_error(integrate(ode, CrankNicolson()), ode) < 1E-13
    end

    @testset "Data Type Consistency" begin
        ode = odeproblem()
        sol = integrate(ode, CrankNicolson())

        @test eltype(sol.q[0]) == Float64
        @test all(isfinite, sol.q[end])
    end

end
