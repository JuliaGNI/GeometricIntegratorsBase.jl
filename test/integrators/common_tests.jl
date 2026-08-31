using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions
using Test

using GeometricSolutions: relative_maximum_error
using GeometricIntegratorsBase: isimplicit
using ..HarmonicOscillator
using ..NonautonomousProblems
using ..NonlinearProblems

# Properties every integrator in this package has to have, asserted for all of them from one
# place rather than restated in each method's test file. What belongs here is what is stated in
# the same words for every method — that unsupported problems are rejected, that a supported one
# integrates to a finite solution of the right type and at the method's order, that formulations
# of one and the same equation give one and the same solution. What does not is anything whose
# rationale is particular to a method: its traits, its cache layout, its energy behaviour, and
# the schemes it is or is not supposed to coincide with. Those stay in the per-method files.

# All problem types the test problems cover, by the name this file refers to them by.
const PROBLEMS = (
    ode = odeproblem,
    pode = podeproblem,
    hode = hodeproblem,
    iode = iodeproblem,
    lode = lodeproblem,
    sode = sodeproblem,
    dae = daeproblem,
    pdae = pdaeproblem,
    hdae = hdaeproblem,
    idae = idaeproblem,
    ldae = ldaeproblem
)

# Every integrator, with the problem types it implements, the order it converges at, and the
# relative error it integrates the harmonic oscillator to at the default timestep. Adding a
# method here is what subjects it to everything below, so a new integrator that forgets to
# reject a constrained problem, or that is a half order off, is caught without writing a test.
const INTEGRATORS = (
    (method = ExplicitEuler(), order = 1, accuracy = 5E-2, problems = (:ode,)),
    (method = ImplicitEuler(), order = 1, accuracy = 5E-2, problems = (:ode,)),
    (method = SymplecticEulerA(), order = 1, accuracy = 5E-2, problems = (:pode, :hode)),
    (method = SymplecticEulerB(), order = 1, accuracy = 5E-2, problems = (:pode, :hode)),
    (method = ImplicitMidpoint(), order = 2, accuracy = 1E-3,
        problems = (:ode, :pode, :hode, :iode, :lode)),
    (method = CrankNicolson(), order = 2, accuracy = 1E-3,
        problems = (:ode, :pode, :hode, :iode, :lode))
)

# Pairs of formulations of one and the same equation: a partitioned or implicit problem and the
# first order system it is equivalent to. A method implemented for both has to produce the very
# same map on the two, which pins every stage time down exactly, whereas a convergence order test
# only detects them to leading order. The harmonic oscillator and the non-autonomous problems are
# linear and separable, so their solves converge in a single Newton step and the blocks of the
# residual Jacobian that couple the two components vanish; the nonlinear problem has neither
# property and pins down the coupled residual and the solver iterates as well.
const EQUIVALENT_FORMULATIONS = (
    (kind = :pode, problem = podeproblem, ode = odeproblem, autonomous = true),
    (kind = :pode, problem = nonautonomous_podeproblem,
        ode = nonautonomous_pode_odeproblem, autonomous = false),
    (kind = :pode, problem = nonlinear_podeproblem,
        ode = nonlinear_pode_odeproblem, autonomous = true),
    (kind = :hode, problem = nonlinear_hodeproblem,
        ode = nonlinear_pode_odeproblem, autonomous = true),
    (kind = :iode, problem = iodeproblem, ode = odeproblem, autonomous = true),
    (kind = :iode, problem = nonautonomous_iodeproblem,
        ode = nonautonomous_iode_odeproblem, autonomous = false),
    (kind = :lode, problem = nonautonomous_lodeproblem,
        ode = nonautonomous_iode_odeproblem, autonomous = false)
)

# Pairs of a problem and the same problem enriched by a Hamiltonian, or by ω and a Lagrangian.
# The additional data does not enter the integration, so both have to give the same solution.
const ENRICHED_FORMULATIONS = (
    (plain = :pode, problem = podeproblem, enriched = :hode, enrichment = hodeproblem),
    (plain = :pode, problem = nonautonomous_podeproblem,
        enriched = :hode, enrichment = nonautonomous_hodeproblem),
    (plain = :pode, problem = nonlinear_podeproblem,
        enriched = :hode, enrichment = nonlinear_hodeproblem),
    (plain = :iode, problem = iodeproblem, enriched = :lode, enrichment = lodeproblem),
    (plain = :iode, problem = nonautonomous_iodeproblem,
        enriched = :lode, enrichment = nonautonomous_lodeproblem)
)

# The timesteps the convergence order is measured over, and the factor the error is expected to
# shrink by from one to the next, which is 2^order for a method of that order.
const TIMESTEPS = (0.1, 0.05, 0.025)

methodname(method) = string(nameof(typeof(method)))

# Maximum relative error against a reference solution, over the momentum as well where the
# problem has one.
function maxerror(sol, ref)
    err = relative_maximum_error(sol, ref)
    hasproperty(err, :p) ? max(err.q, err.p) : err.q
end

# The order of a method is the factor its error shrinks by when the timestep is halved.
function convergence_ratios(method, problem)
    errs = [begin
                prob = problem(; timestep = Δt)
                maxerror(integrate(prob, method), exact_solution(prob))
            end
            for Δt in TIMESTEPS]

    errs[1:(end - 1)] ./ errs[2:end]
end

# The final state of a reference solution of a first order system, accurate enough not to enter
# the error at any of the timesteps above: the implicit midpoint method at Δt = 1E-3 leaves an
# error of order 1E-6, some three orders below the error of a second order method at the smallest
# timestep measured. It is a *fixed* method rather than the one under test, so that the reference
# is available for the methods that do not implement the first order formulation themselves.
function reference_final_state(problem)
    integrate(problem(; timestep = 1E-3), ImplicitMidpoint()).q[end]
end

@testset "$(rpad("Common Integrator Tests", 80))" begin
    @testset "Unsupported Problem Types" begin
        # no method in this package enforces a constraint or knows about the additional equations
        # of a differential algebraic problem, so those must be rejected rather than integrated as
        # if the constraint were absent. the same holds for problem types a method does not
        # implement at all, and either way the error has to name both the method and the problem:
        # `initsolver` raises it for a method that solves, the `integrate_step!` fallback for one
        # that does not
        for integrator in INTEGRATORS, (kind, problem) in pairs(PROBLEMS)

            kind in integrator.problems && continue

            prob = problem()
            err = try
                integrate(prob, integrator.method)
                nothing
            catch e
                e
            end

            @test err isa ArgumentError
            @test occursin(methodname(integrator.method), err.msg)
            @test occursin(string(nameof(typeof(equation(prob)))), err.msg)
        end
    end

    @testset "Supported Problem Types" begin
        # the solution covers the whole time span, is finite throughout, and is of the data type
        # of the problem rather than of whatever the solver happened to compute in
        for integrator in INTEGRATORS, kind in integrator.problems

            prob = PROBLEMS[kind]()
            sol = integrate(prob, integrator.method)

            @test length(sol.q) == length(sol.t)
            @test eltype(sol.q[0]) == datatype(prob)
            @test all(x -> all(isfinite, x), sol.q)

            if hasproperty(sol, :p)
                @test eltype(sol.p[0]) == datatype(prob)
                @test all(x -> all(isfinite, x), sol.p)
            end
        end
    end

    @testset "Integration Accuracy" begin
        for integrator in INTEGRATORS, kind in integrator.problems

            prob = PROBLEMS[kind]()
            @test maxerror(integrate(prob, integrator.method), exact_solution(prob)) <
                  integrator.accuracy

            # a converged solve is silent, so tight tolerances must not produce solver warnings
            isimplicit(integrator.method) || continue

            sol = @test_nowarn integrate(prob, integrator.method;
                x_abstol = 1e-12,
                f_abstol = 1e-12
            )
            @test maxerror(sol, exact_solution(prob)) < integrator.accuracy
        end
    end

    @testset "Convergence Order" begin
        # halving the timestep reduces the error by 2^order, so the ratio identifies the order
        for integrator in INTEGRATORS, kind in integrator.problems

            ratios = convergence_ratios(integrator.method, PROBLEMS[kind])
            @test all(isapprox.(ratios, 2^integrator.order; atol = 2E-1))
        end
    end

    @testset "Non-Autonomous Convergence Order" begin
        # every problem of the harmonic oscillator is autonomous, so none of them can tell the
        # stage times of a method apart: evaluating a vector field at tₙ or at tₙ₊₁ gives the same
        # result there, and a method that gets them wrong still converges at its nominal order.
        # on an explicitly time dependent problem it loses an order, which is what is measured
        for integrator in INTEGRATORS
            :ode in integrator.problems || continue

            errs = [begin
                        prob = nonautonomous_odeproblem(; timestep = Δt)
                        maxerror(integrate(prob, integrator.method), nonautonomous_ode_solution(prob))
                    end
                    for Δt in TIMESTEPS]

            @test all(isapprox.(errs[1:(end - 1)] ./ errs[2:end], 2^integrator.order; atol = 3E-1))
        end

        # the partitioned and implicit formulations have no closed form solution here, so the
        # reference is the equivalent first order system integrated to far higher accuracy
        for formulation in EQUIVALENT_FORMULATIONS
            formulation.autonomous && continue

            ref = reference_final_state(formulation.ode)

            for integrator in INTEGRATORS
                formulation.kind in integrator.problems || continue

                errs = [begin
                            sol = integrate(formulation.problem(; timestep = Δt), integrator.method)
                            max(abs(sol.q[end][1] - ref[1]), abs(sol.p[end][1] - ref[2]))
                        end
                        for Δt in TIMESTEPS]

                @test all(isapprox.(errs[1:(end - 1)] ./ errs[2:end], 2^integrator.order; atol = 3E-1))
            end
        end
    end

    @testset "Equivalence with the ODE Formulation" begin
        for integrator in INTEGRATORS, formulation in EQUIVALENT_FORMULATIONS
            # only a method implemented for both formulations can be compared on them
            (:ode in integrator.problems && formulation.kind in integrator.problems) ||
                continue

            sp = integrate(formulation.problem(), integrator.method; x_abstol = 1e-14, f_abstol = 1e-14)
            so = integrate(formulation.ode(), integrator.method; x_abstol = 1e-14, f_abstol = 1e-14)

            @test maximum(abs(sp.q[n][1] - so.q[n][1]) for n in axes(sp.q, 1)) < 1E-12
            @test maximum(abs(sp.p[n][1] - so.q[n][2]) for n in axes(sp.q, 1)) < 1E-12
        end
    end

    @testset "Hamiltonian and Lagrangian Formulations" begin
        for integrator in INTEGRATORS, formulation in ENRICHED_FORMULATIONS

            (formulation.plain in integrator.problems &&
             formulation.enriched in integrator.problems) || continue

            plain = integrate(formulation.problem(), integrator.method)
            enriched = integrate(formulation.enrichment(), integrator.method)

            @test plain.q == enriched.q
            @test plain.p == enriched.p
        end
    end
end
