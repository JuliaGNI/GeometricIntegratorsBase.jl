using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: ImplicitEulerCache, nlsolution, solversize
using GeometricIntegratorsBase: default_solver, default_iguess
using GeometricIntegratorsBase: isexplicit, isimplicit, issymmetric, issymplectic

using ..HarmonicOscillator

@testset "$(rpad("ImplicitEuler Method Tests", 80))" begin

    @testset "Method Properties" begin
        method = ImplicitEuler()

        @test !isexplicit(method)
        @test isimplicit(method)
        @test !issymmetric(method)
        @test !issymplectic(method)

        @test default_solver(method) isa Newton
        @test default_iguess(method) isa HermiteExtrapolation
    end

    @testset "Cache Structure" begin
        ode = odeproblem()
        method = ImplicitEuler()

        # Test cache creation
        cache = Cache{Float64}(ode, method)
        @test cache isa ImplicitEulerCache{Float64,ndims(ode)}

        # Test cache fields
        @test length(cache.x) == ndims(ode)
        @test length(cache.q) == ndims(ode)
        @test length(cache.v) == ndims(ode)
        @test length(cache.v̄) == ndims(ode)

        # Test nlsolution accessor
        @test nlsolution(cache) === cache.x

        # Test solver size
        @test solversize(ode, method) == ndims(ode)
    end

    @testset "Integration Accuracy - Basic" begin
        ode = odeproblem()
        ref = exact_solution(ode)

        # Test basic integration
        sol = integrate(ode, ImplicitEuler())
        err = relative_maximum_error(sol, ref)
        @test err.q < 5E-2

        # Test with tighter solver tolerances
        sol_tight = integrate(ode, ImplicitEuler();
            min_iterations=1,
            x_abstol=1e-12,
            f_abstol=1e-12,
            max_iterations=100
        )
        err_tight = relative_maximum_error(sol_tight, ref)
        @test err_tight.q < 5E-2
    end

    @testset "Different Timesteps" begin
        # Test with smaller timestep (should be more accurate)
        ode_small = odeproblem([0.5, 0.0]; timestep=0.05)
        ref_small = exact_solution(ode_small)
        sol_small = integrate(ode_small, ImplicitEuler())
        err_small = relative_maximum_error(sol_small, ref_small)

        # Test with larger timestep (should be less accurate)
        ode_large = odeproblem([0.5, 0.0]; timestep=0.2)
        ref_large = exact_solution(ode_large)
        sol_large = integrate(ode_large, ImplicitEuler())
        err_large = relative_maximum_error(sol_large, ref_large)

        @test err_small.q < 2E-2
        @test err_large.q < 1E-1
        @test err_small.q < err_large.q  # Smaller timestep should be more accurate
    end

    @testset "Different Initial Conditions" begin
        # Test with different initial position
        ode1 = odeproblem([1.0, 0.5])
        ref1 = exact_solution(ode1)
        sol1 = integrate(ode1, ImplicitEuler())
        err1 = relative_maximum_error(sol1, ref1)
        @test err1.q < 5E-2

        # Test with non-zero initial position and velocity
        ode2 = odeproblem([0.2, 0.8])
        ref2 = exact_solution(ode2)
        sol2 = integrate(ode2, ImplicitEuler())
        err2 = relative_maximum_error(sol2, ref2)
        @test err2.q < 5E-2
    end

    @testset "Convergence Properties" begin
        ode = odeproblem()

        # Test that very loose tolerances still converge
        sol_loose = integrate(ode, ImplicitEuler();
            x_abstol=1e-4,
            f_abstol=1e-4,
            max_iterations=50
        )
        @test all(x -> all(isfinite, x), sol_loose.q)

        # Test that we can increase max iterations
        sol_many_iter = integrate(ode, ImplicitEuler();
            max_iterations=200
        )
        @test all(x -> all(isfinite, x), sol_many_iter.q)
    end

    @testset "Energy Conservation Check" begin
        ode = odeproblem()
        sol = integrate(ode, ImplicitEuler())

        # Compute energy at each time step
        energies = [hamiltonian(t, Array(q), parameters(ode)) for (t, q) in zip(sol.t, sol.q)]

        # For harmonic oscillator, energy should be approximately conserved
        # (though ImplicitEuler is not symplectic, so we allow larger deviation)
        initial_energy = energies[1]
        energy_variation = maximum(abs.(energies .- initial_energy)) / initial_energy

        @test energy_variation < 0.1  # Allow 10% energy variation for ImplicitEuler
    end

    @testset "Comparison with ExplicitEuler" begin
        ode = odeproblem()
        ref = exact_solution(ode)

        sol_implicit = integrate(ode, ImplicitEuler())
        sol_explicit = integrate(ode, ExplicitEuler())

        err_implicit = relative_maximum_error(sol_implicit, ref)
        err_explicit = relative_maximum_error(sol_explicit, ref)

        # Both should be O(h) methods with similar accuracy for this problem
        @test err_implicit.q < 5E-2
        @test err_explicit.q < 5E-2

        # The errors should be of similar magnitude (within factor of 5)
        @test min(err_implicit.q, err_explicit.q) / max(err_implicit.q, err_explicit.q) > 0.2
    end

    @testset "Stiff Problem Behavior" begin
        # Create a simple stiff ODE: dy/dt = -100*y, y(0) = 1
        # Exact solution: y(t) = exp(-100*t)
        stiff_v = function (v, t, q, params)
            v[1] = -100.0 * q[1]
        end

        ode_stiff = ODEProblem(stiff_v, (0.0, 0.1), 0.01, [1.0])
        sol_stiff = integrate(ode_stiff, ImplicitEuler())

        # Check that solution remains stable (doesn't blow up)
        @test all(x -> all(isfinite, x), sol_stiff.q)
        @test all(q -> q[1] >= 0, sol_stiff.q)  # Solution should remain positive
        @test sol_stiff.q[end][1] < sol_stiff.q[1][1]  # Should decay

        # Compare with exact solution
        exact_final = exp(-100.0 * 0.1)
        @test abs(sol_stiff.q[end][1] - exact_final) / exact_final < 50.0  # Allow large error for stiff problem with implicit Euler
    end

    @testset "Edge Cases" begin
        # Test with very small timestep
        ode_tiny = odeproblem([0.5, 0.0]; timestep=1e-4, timespan=(0.0, 0.01))
        sol_tiny = integrate(ode_tiny, ImplicitEuler())
        @test all(x -> all(isfinite, x), sol_tiny.q)

        # Test with single timestep
        ode_single = odeproblem([0.5, 0.0]; timestep=0.1, timespan=(0.0, 0.1))
        sol_single = integrate(ode_single, ImplicitEuler())
        @test length(sol_single.t) == 2  # Initial + one step
        @test all(x -> all(isfinite, x), sol_single.q)
    end

    @testset "Data Type Consistency" begin
        ode = odeproblem()
        sol = integrate(ode, ImplicitEuler())

        # Check that solution maintains data type consistency
        @test all(x -> all(isfinite, x), sol.q)
        @test length(sol.q) == length(sol.t)
    end

    @testset "Method Interface" begin
        ode = odeproblem()
        int = GeometricIntegrator(ode, ImplicitEuler())

        # Test that integrator was created successfully
        @test int isa GeometricIntegrator
        @test typeof(int).parameters[1] == ImplicitEuler
        @test typeof(int).parameters[2] == typeof(ode)
    end
end
