using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: SymplecticEulerCache, CacheType, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic
using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonautonomousProblems: nonautonomous_pode_v, nonautonomous_pode_f


# Reference implementations of the two schemes, written out directly so that the stage times
# are stated here rather than taken from the code under test. The times are read off the
# solution rather than recomputed, so that any disagreement is due to the scheme alone.
function reference_symplectic_euler_A(prob, sol)
    q, p = copy(initial_conditions(prob).q), copy(initial_conditions(prob).p)
    v, f = zero(q), zero(p)
    h = timestep(prob)

    for n in 1:lastindex(sol.t)
        nonautonomous_pode_f(f, sol.t[n-1], q, p, parameters(prob))
        p .+= h .* f
        nonautonomous_pode_v(v, sol.t[n], q, p, parameters(prob))
        q .+= h .* v
    end

    (q=q, p=p)
end

function reference_symplectic_euler_B(prob, sol)
    q, p = copy(initial_conditions(prob).q), copy(initial_conditions(prob).p)
    v, f = zero(q), zero(p)
    h = timestep(prob)

    for n in 1:lastindex(sol.t)
        nonautonomous_pode_v(v, sol.t[n-1], q, p, parameters(prob))
        q .+= h .* v
        nonautonomous_pode_f(f, sol.t[n], q, p, parameters(prob))
        p .+= h .* f
    end

    (q=q, p=p)
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

    @testset "Unsupported Problem Types" begin
        # the methods neither enforce a constraint nor know about the additional equations of a
        # differential algebraic problem, so those must be rejected rather than integrated as if
        # the constraint were absent. in contrast to the implicit methods the error is a
        # MethodError from `integrate_step!` rather than the ArgumentError of `initsolver`, which
        # an explicit method never reaches: its `default_solver` is `NoSolver`
        for method in (SymplecticEulerA(), SymplecticEulerB())
            for prob in (pdaeproblem(), hdaeproblem())
                @test_throws MethodError integrate(prob, method)
            end
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

    @testset "Non-Autonomous Stage Times" begin
        # the harmonic oscillator is autonomous and cannot tell t̄ from t̄ + Δt apart, so the
        # stage times are pinned down here against a reference implementation of each scheme
        pode = nonautonomous_podeproblem()

        solA = integrate(pode, SymplecticEulerA())
        solB = integrate(pode, SymplecticEulerB())

        refA = reference_symplectic_euler_A(pode, solA)
        refB = reference_symplectic_euler_B(pode, solB)

        @test solA.q[end] == refA.q
        @test solA.p[end] == refA.p
        @test solB.q[end] == refB.q
        @test solB.p[end] == refB.p

        # the reference implementations must themselves be distinguishable, otherwise the
        # comparison above would hold for either choice of stage times
        @test refA.q != refB.q
        @test refA.p != refB.p

        # HODE and PODE agree on a non-autonomous problem too
        hode = nonautonomous_hodeproblem()
        for method in (SymplecticEulerA(), SymplecticEulerB())
            @test integrate(pode, method).q == integrate(hode, method).q
            @test integrate(pode, method).p == integrate(hode, method).p
        end
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

            # the longer run samples the same oscillation more densely and so comes closer to
            # the peak of the envelope, which leaves the two maxima differing by ~3E-7 here.
            # a relative tolerance keeps the margin independent of the problem parameters;
            # any genuine drift would show up as a ratio of about ten, not as 1E-5
            @test err_long ≈ err_short rtol = 1E-3
        end

        # for reference: the explicit Euler method on the equivalent ODE blows up
        ode_long = odeproblem(; timespan=(0.0, 1000.0))
        @test max_energy_error(integrate(ode_long, ExplicitEuler()), ode_long) > 1E3
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
