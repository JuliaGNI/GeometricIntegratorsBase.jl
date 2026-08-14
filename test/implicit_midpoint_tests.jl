using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitMidpointCache, ImplicitMidpointPODECache, ImplicitMidpointIODECache
using GeometricIntegratorsBase: CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using SimpleSolvers: Newton
using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonlinearProblems


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

    @testset "Unsupported Problem Types" begin
        # the method neither enforces a constraint nor knows about the additional equations of a
        # differential algebraic problem, so those must be rejected rather than integrated as if
        # the constraint were absent. the same holds for problem types the method does not
        # implement at all, and either way the error has to say so
        for prob in (daeproblem(), pdaeproblem(), hdaeproblem(), idaeproblem(), ldaeproblem())
            err = try
                integrate(prob, ImplicitMidpoint())
                nothing
            catch e
                e
            end

            @test err isa ArgumentError
            @test occursin("ImplicitMidpoint", err.msg)
            @test occursin(string(nameof(typeof(equation(prob)))), err.msg)
        end
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

    @testset "PODE/HODE Problems" begin

        @testset "Cache Structure" begin
            pode = podeproblem()
            method = ImplicitMidpoint()

            cache = Cache{Float64}(pode, method)
            @test cache isa ImplicitMidpointPODECache{Float64}

            # both stage vector fields are solved for, so the solution vector holds two of them
            @test length(cache.x) == length(initial_conditions(pode).q) + length(initial_conditions(pode).p)
            @test solversize(method, pode) == length(cache.x)

            @test axes(cache.q) == axes(initial_conditions(pode).q)
            @test axes(cache.p) == axes(initial_conditions(pode).p)
            @test axes(cache.v) == axes(initial_conditions(pode).q)
            @test axes(cache.f) == axes(initial_conditions(pode).p)

            @test nlsolution(cache) === cache.x

            @test CacheType(Float64, pode, method) == ImplicitMidpointPODECache{Float64}
            @test CacheType(Float32, pode, method) == ImplicitMidpointPODECache{Float32}

            # the Hamiltonian formulation is integrated with the very same cache
            @test Cache{Float64}(hodeproblem(), method) isa ImplicitMidpointPODECache{Float64}
        end

        @testset "Integration Accuracy" begin
            pode = podeproblem()
            ref = exact_solution(pode)

            sol = integrate(pode, ImplicitMidpoint())
            err = relative_maximum_error(sol, ref)
            @test err.q < 1E-3
            @test err.p < 1E-3

            # a converged solve is silent, so tight tolerances must not produce solver warnings
            sol_tight = @test_nowarn integrate(pode, ImplicitMidpoint();
                x_abstol=1e-12,
                f_abstol=1e-12,
            )
            @test relative_maximum_error(sol_tight, ref).q < 1E-3

            # the Hamiltonian formulation adds the Hamiltonian, but describes the same equation,
            # so it integrates to the same solution
            hode = hodeproblem()
            @test integrate(hode, ImplicitMidpoint()).q == sol.q
            @test integrate(hode, ImplicitMidpoint()).p == sol.p
        end

        @testset "Equivalence with the ODE Formulation" begin
            # a partitioned problem whose two components are collected into a single state vector
            # is an ordinary differential equation, and the method applied to either of the two is
            # the very same map. this pins the stage times down exactly, whereas a convergence
            # order test only detects them to leading order.
            # the harmonic oscillator and the non-autonomous problem are both linear and separable,
            # so the solve converges in a single Newton step and the blocks of the residual
            # Jacobian that couple the two components vanish. the nonlinear problem has neither
            # property and pins down the coupled residual and the solver iterates as well
            for (pode, ode) in ((podeproblem(), odeproblem()),
                (nonautonomous_podeproblem(), nonautonomous_pode_odeproblem()),
                (nonlinear_podeproblem(), nonlinear_pode_odeproblem()))

                sp = integrate(pode, ImplicitMidpoint(); x_abstol=1e-14, f_abstol=1e-14)
                so = integrate(ode, ImplicitMidpoint(); x_abstol=1e-14, f_abstol=1e-14)

                @test maximum(abs(sp.q[n][1] - so.q[n][1]) for n in axes(sp.q, 1)) < 1E-12
                @test maximum(abs(sp.p[n][1] - so.q[n][2]) for n in axes(sp.q, 1)) < 1E-12
            end
        end

        @testset "Convergence Order" begin
            # second order, so halving the timestep should reduce the error by a factor of four
            errs = [
                begin
                    prob = podeproblem(; timestep=Δt)
                    relative_maximum_error(integrate(prob, ImplicitMidpoint()), exact_solution(prob)).q
                end for Δt in (0.1, 0.05, 0.025)
            ]

            @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=2E-1))
        end

        @testset "Non-Autonomous Convergence Order" begin
            # the harmonic oscillator is autonomous, so it cannot detect a wrong stage time.
            # the force of this problem depends on time explicitly, and the reference is computed
            # from the equivalent first order system at a much smaller timestep
            ref = integrate(nonautonomous_pode_odeproblem(; timestep=1E-3), ImplicitMidpoint())
            qref, pref = ref.q[end][1], ref.q[end][2]

            errs = [
                begin
                    sol = integrate(nonautonomous_podeproblem(; timestep=Δt), ImplicitMidpoint())
                    max(abs(sol.q[end][1] - qref), abs(sol.p[end][1] - pref))
                end for Δt in (0.1, 0.05, 0.025)
            ]

            @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 4; atol=3E-1))
        end

        @testset "Energy Conservation" begin
            # the implicit midpoint method conserves quadratic invariants exactly, so the energy
            # of the harmonic oscillator is preserved up to round-off, even over long times
            pode = podeproblem(; timespan=(0.0, 100.0))
            @test max_energy_error(integrate(pode, ImplicitMidpoint()), pode) < 1E-13
        end

        @testset "Data Type Consistency" begin
            pode = podeproblem()
            sol = integrate(pode, ImplicitMidpoint())

            @test eltype(sol.q[0]) == Float64
            @test eltype(sol.p[0]) == Float64
            @test all(isfinite, sol.q[end])
            @test all(isfinite, sol.p[end])
        end

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
                (nonautonomous_iodeproblem(), nonautonomous_iode_odeproblem()),
                (nonautonomous_lodeproblem(), nonautonomous_iode_odeproblem()))

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
