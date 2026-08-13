using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitMidpointCache, ImplicitMidpointIODECache
using GeometricIntegratorsBase: CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator
using ..NonautonomousProblems


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

    @testset "Non-Autonomous Convergence Order" begin
        # the harmonic oscillator is autonomous, so it cannot detect a wrong stage time.
        # evaluating the vector field anywhere but at t̄ + Δt/2 costs an order of accuracy,
        # which this problem does detect
        errs = [
            begin
                prob = nonautonomous_odeproblem(; timestep=Δt)
                relative_maximum_error(integrate(prob, ImplicitMidpoint()),
                    nonautonomous_ode_solution(prob)).q
            end for Δt in (0.1, 0.05, 0.025)
        ]

        @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=3E-1))
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

    @testset "IODE/LODE Problems" begin

        @testset "Cache Structure" begin
            iode = iodeproblem()
            method = ImplicitMidpoint()

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
        end

        @testset "Integration Accuracy" begin
            iode = iodeproblem()
            ref = exact_solution(iode)

            sol = integrate(iode, ImplicitMidpoint())
            err = relative_maximum_error(sol, ref)
            @test err.q < 1E-3
            @test err.p < 1E-3

            # a converged solve is silent, so tight tolerances must not produce solver warnings
            sol_tight = @test_nowarn integrate(iode, ImplicitMidpoint();
                x_abstol=1e-12,
                f_abstol=1e-12,
            )
            @test relative_maximum_error(sol_tight, ref).q < 1E-3

            # the Lagrangian formulation adds ω and the Lagrangian, but describes the same
            # equation, so it integrates to the same solution
            lode = lodeproblem()
            @test integrate(lode, ImplicitMidpoint()).q == sol.q
            @test integrate(lode, ImplicitMidpoint()).p == sol.p
        end

        @testset "Equivalence with the ODE Formulation" begin
            # the momentum map of both problems is regular, so the implicit equation is
            # equivalent to a first order system, and the implicit midpoint method applied to
            # either of the two is the very same map. this pins the stage times down exactly,
            # whereas a convergence order test only detects them to leading order
            for (iode, ode) in ((iodeproblem(), odeproblem()),
                (nonautonomous_iodeproblem(), nonautonomous_iode_odeproblem()))

                si = integrate(iode, ImplicitMidpoint(); x_abstol=1e-14, f_abstol=1e-14)
                so = integrate(ode, ImplicitMidpoint(); x_abstol=1e-14, f_abstol=1e-14)

                @test maximum(abs(si.q[n][1] - so.q[n][1]) for n in axes(si.q, 1)) < 1E-12
                @test maximum(abs(si.p[n][1] - so.q[n][2]) for n in axes(si.q, 1)) < 1E-12
            end
        end

        @testset "Convergence Order" begin
            # second order, so halving the timestep should reduce the error by a factor of four
            errs = [
                begin
                    prob = iodeproblem(; timestep=Δt)
                    relative_maximum_error(integrate(prob, ImplicitMidpoint()), exact_solution(prob)).q
                end for Δt in (0.1, 0.05, 0.025)
            ]

            @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=2E-1))
        end

        @testset "Non-Autonomous Convergence Order" begin
            # the harmonic oscillator is autonomous, so it cannot detect a wrong stage time.
            # both ϑ and f of this problem depend on time explicitly, and the reference is
            # computed from the equivalent first order system at a much smaller timestep
            ref = integrate(nonautonomous_iode_odeproblem(; timestep=1E-3), ImplicitMidpoint())
            qref, pref = ref.q[end][1], ref.q[end][2]

            errs = [
                begin
                    sol = integrate(nonautonomous_iodeproblem(; timestep=Δt), ImplicitMidpoint())
                    max(abs(sol.q[end][1] - qref), abs(sol.p[end][1] - pref))
                end for Δt in (0.1, 0.05, 0.025)
            ]

            @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=3E-1))
        end

        @testset "Energy Conservation" begin
            # the implicit midpoint method conserves quadratic invariants exactly, so the energy
            # of the harmonic oscillator is preserved up to round-off, even over long times
            iode = iodeproblem(; timespan=(0.0, 100.0))
            @test max_energy_error(integrate(iode, ImplicitMidpoint()), iode) < 1E-13
        end

        @testset "Data Type Consistency" begin
            iode = iodeproblem()
            sol = integrate(iode, ImplicitMidpoint())

            @test eltype(sol.q[0]) == Float64
            @test eltype(sol.p[0]) == Float64
            @test all(isfinite, sol.q[end])
            @test all(isfinite, sol.p[end])
        end

    end

end
