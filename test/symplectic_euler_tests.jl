using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: SymplecticEulerCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using ..HarmonicOscillator


# maximum relative energy error along a partitioned solution
function max_energy_error(sol, prob)
    maximum(abs, compute_energy_error(sol.t, sol.q, sol.p, parameters(prob))[2])
end


@testset "$(rpad("SymplecticEuler Method Tests", 80))" begin

    @testset "Method Properties" begin
        for method in (SymplecticEulerA(), SymplecticEulerB())
            @test method isa GeometricIntegratorsBase.SymplecticEulerMethod

            @test isexplicit(method)
            @test !isimplicit(method)
            @test !issymmetric(method)
            @test issymplectic(method)

            @test ispodemethod(method)
            @test !isodemethod(method)

            # explicit methods need neither a solver nor an initial guess
            @test default_solver(method) isa NoSolver
            @test default_iguess(method) isa NoInitialGuess
            @test solversize(method, podeproblem()) == 0
        end
    end

    @testset "Cache Structure" begin
        pode = podeproblem()

        for method in (SymplecticEulerA(), SymplecticEulerB())
            cache = Cache{Float64}(pode, method)
            @test cache isa SymplecticEulerCache{Float64}

            @test axes(cache.v) == axes(initial_conditions(pode).q)
            @test axes(cache.f) == axes(initial_conditions(pode).p)

            # explicit methods have no nonlinear solver solution
            @test ismissing(nlsolution(cache))

            @test CacheType(Float64, pode, method) == SymplecticEulerCache{Float64}
            @test CacheType(Float32, pode, method) == SymplecticEulerCache{Float32}
        end
    end

    @testset "Integration Accuracy" begin
        pode = podeproblem()
        ref = exact_solution(pode)

        for method in (SymplecticEulerA(), SymplecticEulerB())
            sol = integrate(pode, method)
            err = relative_maximum_error(sol, ref)
            @test err.q < 5E-2
            @test err.p < 5E-2
        end
    end

    @testset "HODE and PODE Agreement" begin
        pode = podeproblem()
        hode = hodeproblem()

        for method in (SymplecticEulerA(), SymplecticEulerB())
            psol = integrate(pode, method)
            hsol = integrate(hode, method)
            @test psol.q == hsol.q
            @test psol.p == hsol.p
        end
    end

    @testset "A and B are Different Methods" begin
        # note that for p₀ = 0 the position sequence of B is that of A shifted by one step,
        # which makes the momentum sequences of the two methods coincide, so a generic initial
        # condition is needed here
        pode = podeproblem([0.5], [0.5])

        solA = integrate(pode, SymplecticEulerA())
        solB = integrate(pode, SymplecticEulerB())

        @test solA.q != solB.q
        @test solA.p != solB.p
    end

    @testset "Convergence Order" begin
        # both methods are first order, so halving the timestep should halve the error
        for method in (SymplecticEulerA(), SymplecticEulerB())
            errs = [
                begin
                    prob = podeproblem(; timestep=Δt)
                    relative_maximum_error(integrate(prob, method), exact_solution(prob)).q
                end for Δt in (0.1, 0.05, 0.025)
            ]

            @test all(isapprox.(errs[1:end-1] ./ errs[2:end], 2; atol=2E-1))
        end
    end

    @testset "Energy Conservation" begin
        # symplectic methods show a bounded, oscillatory energy error, i.e., the error over a
        # ten times longer integration is essentially the same, whereas the energy error of a
        # non-symplectic method of the same order grows without bound
        pode_short = podeproblem(; timespan=(0.0, 100.0))
        pode_long = podeproblem(; timespan=(0.0, 1000.0))

        for method in (SymplecticEulerA(), SymplecticEulerB())
            err_short = max_energy_error(integrate(pode_short, method), pode_short)
            err_long = max_energy_error(integrate(pode_long, method), pode_long)

            @test err_short < 5E-2
            @test err_long ≈ err_short atol = 1E-6
        end

        # for reference: the explicit Euler method on the equivalent ODE blows up
        ode_long = odeproblem(; timespan=(0.0, 1000.0))
        sol = integrate(ode_long, ExplicitEuler())
        @test maximum(abs, compute_energy_error(sol.t, sol.q, parameters(ode_long))[2]) > 1E3
    end

    @testset "Data Type Consistency" begin
        pode = podeproblem()

        for method in (SymplecticEulerA(), SymplecticEulerB())
            sol = integrate(pode, method)
            @test eltype(sol.q[0]) == Float64
            @test eltype(sol.p[0]) == Float64
            @test all(isfinite, sol.q[end])
            @test all(isfinite, sol.p[end])
        end
    end

end
