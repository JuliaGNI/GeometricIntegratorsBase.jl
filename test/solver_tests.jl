using GeometricIntegratorsBase
using Test

import GeometricIntegratorsBase: DEFAULT_F_STALL_WINDOW
import GeometricIntegratorsBase: default_linesearch, default_options, initsolver

using SimpleSolvers: Backtracking, F_STALL_REPORT_MINIMUM
using SimpleSolvers: NewtonSolver, NonlinearSolverStatus, SolverState, isconverged, solve_with_status!
using SimpleSolvers: dominant_linesearch_outcome

using ..HarmonicOscillator


struct TestMethod <: GeometricMethod end

prob = odeproblem()

caches = CacheDict(prob, TestMethod())


@test initsolver(NoSolver(), TestMethod(), caches) == NoSolver()


# the untyped `default_linesearch` is the downstream hook and must keep working unchanged,
# while the typed one builds the method at the working precision directly
@test default_linesearch() isa Backtracking{Float64}
@test default_linesearch(TestMethod()) isa Backtracking{Float64}
@test default_linesearch(Float32) isa Backtracking{Float32}
@test default_linesearch(Float32, TestMethod()) isa Backtracking{Float32}

@test default_linesearch(BigFloat) isa Backtracking{BigFloat}

# the typed method must not intercept a call meant for the untyped hook: passing a method *type*
# instead of an instance still resolves to the untyped one, as it did before the typed one existed
@test default_linesearch(TestMethod) isa Backtracking{Float64}
@test which(default_linesearch, Tuple{Type{TestMethod}}) !== which(default_linesearch, Tuple{Type{Float32}})

# ... while the constraint stays wide enough for a real type that is not an `AbstractFloat`: the
# `ForwardDiff.Dual` a caller differentiating through an integration arrives with has to reach the
# typed method, which is what `T<:Real` rather than `T<:AbstractFloat` buys
@test which(default_linesearch, Tuple{Type{Rational{Int}}}) === which(default_linesearch, Tuple{Type{Float32}})

# the expansion phase of SimpleSolvers 0.11 is deliberately left off, since every method here
# solves with `Newton()`, whose direction is already scaled like a Newton step
@test default_linesearch().expand == false
@test default_linesearch(Float32).expand == false


# a solve whose residual goes nowhere is bounded per time step rather than spending
# `max_iterations` on every one of them
@test default_options(TestMethod(), prob).f_stall_window == DEFAULT_F_STALL_WINDOW
@test default_options(TestMethod(), prob).min_iterations == 1

# the window has to stay above the count below which SimpleSolvers holds the no-progress
# proportion to be no evidence of anything: giving up on less than it will report on is incoherent
@test DEFAULT_F_STALL_WINDOW > F_STALL_REPORT_MINIMUM


# `check_solver_status` is the one place the outcome of a step's solve is acted on. Its default is
# silent — SimpleSolvers reports a failed solve itself, and saying so again would double every
# message in a time-stepping loop — so what is asserted here is exactly that: that it says nothing
# for a converged *and* for a non-converged status, and that it hands the status straight back, so
# that a call site can chain it and an override can decide on the same value the caller sees.
# The second argument is passed as `nothing` because the fallback is untyped in it and must stay
# so; that it can be narrowed to dispatch per method family is a property of Julia, not of this
# function, and is not what these assertions are for.
let
    F!(y, x, params) = y .= x .^ 2 .- 2
    x = [1.0]
    s = NewtonSolver(x, similar(x); F=F!, verbosity=0)
    state = SolverState(s)
    st = solve_with_status!(x, s, state)

    @test st isa NonlinearSolverStatus
    @test isconverged(st)

    # returns its first argument unchanged, and says nothing while doing it
    @test @test_nowarn(check_solver_status(st, nothing)) === st

    # a non-converged status is *also* silent by default: the report is SimpleSolvers' to make
    y = [1.0]
    s2 = NewtonSolver(y, similar(y); F=F!, verbosity=0, max_iterations=1, min_iterations=1)
    state2 = SolverState(s2)
    st2 = solve_with_status!(y, s2, state2)
    @test !isconverged(st2)
    @test @test_nowarn(check_solver_status(st2, nothing)) === st2
end


# The status is the *only* channel a rejected line search has. Up to SimpleSolvers 0.11 a
# `NonlinearSolver` called `linesearch_warnings` from inside its own iteration, so the rejection
# was logged whether or not the caller looked; 0.12 removed that call, and the outcome now reaches
# the caller through the tally the status carries and nowhere else. That is the reason every
# `integrate_step!` here moved from `solve!` to `solve_with_status!`, so it is what has to be
# asserted: that what arrives at `check_solver_status` still says the line search was rejected.
#
# The question is asked with `count_floor = false`, because `linesearch_failures` alone would be
# vacuous here — a *converged* solve reaches the merit's round-off floor on its last step, which
# counts as a failure and would make a bare `> 0` pass on a solve where nothing went wrong. The
# contrast against the converged solve below is what gives the assertion its content.
let
    F!(y, x, params) = y .= x .^ 2 .- 2

    # the true Jacobian with its sign flipped, so the Newton direction points uphill and the line
    # search has nothing to accept; `max_iterations` is capped only to keep the test cheap
    badJ!(J, x, params) = J[1, 1] = -2 * x[1]

    x = [1.0]
    s = NewtonSolver(x, similar(x); F=F!, DF! = badJ!, verbosity=0, max_iterations=8)
    st = solve_with_status!(x, s, SolverState(s))

    @test !isconverged(st)
    @test dominant_linesearch_outcome(st, false) !== nothing

    # and it survives the hook the call sites hand it to, still silently
    @test @test_nowarn(check_solver_status(st, nothing)) === st

    # the same question of a solve whose line search was never rejected
    y = [1.0]
    s2 = NewtonSolver(y, similar(y); F=F!, verbosity=0)
    st2 = solve_with_status!(y, s2, SolverState(s2))

    @test isconverged(st2)
    @test dominant_linesearch_outcome(st2, false) === nothing
end
