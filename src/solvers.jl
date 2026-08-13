
"""
    default_linesearch([T,] method=nothing)

The line search a `method` uses by default, built at the working type `T` where one is given.

`Backtracking()` is always `Backtracking{Float64}`, so a `Float32` integration makes the
`Linesearch` constructor do a `change_precision` the caller could have avoided by building the
method at the working type in the first place. Both forms exist: the untyped one is the downstream
hook (`GeometricIntegrators` calls it for DIRK) and its signature must not change.

Neither form is called inside this package. They are the hook downstream methods reach for, which
is why the tests assert both signatures and the choice below.

`expand = true` is deliberately not set. The expansion phase of `Backtracking` is opt-in upstream
because it costs a merit evaluation — a full residual evaluation here — on a direction whose scale
is wrong, and gains nothing on one whose scale is right: a Newton step already sits at its model
minimum. Every method in this package solves with `Newton()`.
"""
default_linesearch(method=nothing) = Backtracking()
default_linesearch(::Type{T}, method=nothing) where {T<:Number} = Backtracking(T)

"""
The window, in iterations, after which a nonlinear solve that is not descending towards `f_abstol`
is given up on; the `f_stall_window` of [`default_options`](@ref).

Without it, a solve whose residual sits on a floor above `f_abstol` while its iterate keeps moving
normally is invisible: no convergence branch fires, `max_stalls` cannot see it (the step is not
small), and the solve spends `max_iterations` — 1000 by default — on *every* time step.

50 is the conservative point: it gives up only on a residual that is essentially flat, which no
`Newton()` solve on these residuals is, and it stays above `SimpleSolvers.F_STALL_REPORT_MINIMUM`,
below which SimpleSolvers itself holds the no-progress proportion to be no evidence of anything.
See `SimpleSolvers.F_STALL_WINDOW` for the linear rate a given window corresponds to at a given
`f_stall_factor`.

The one case it is *not* conservative for is a solve that converges linearly at a rate close to
that threshold — a `Picard()` fixed-point iteration on a stiff problem near its convergence limit,
say. Every method in this package uses `Newton()`, but `default_options` is reached with whatever
`solvermethod` the caller hands to `GeometricIntegrator`, so such a caller has to raise this.
"""
const DEFAULT_F_STALL_WINDOW = 50

default_options(method::GeometricMethod, problem::GeometricProblem) = (
    min_iterations=1,
    f_abstol=max(8, solversize(method, problem)) * eps(datatype(problem)),
    # Giving up and converging remain mutually exclusive (`no_progress` is gated on
    # `!residual_small`), the solve still reports the residual it achieved against the tolerance
    # it was asked for, and `default_options` is merged rather than replaced, so a caller who
    # disagrees overrides this with one keyword.
    f_stall_window=DEFAULT_F_STALL_WINDOW,
)

initsolver(::SolverMethod, ::GeometricMethod, ::CacheDict; kwargs...) = NoSolver()

# create nonlinear solver
function initsolver(solvermethod::NonlinearSolverMethod, method::GeometricMethod, caches::CacheDict; kwargs...)
    # a method that has no cache for the problem at hand does not implement it. without this
    # check the generic `NoCache` fallback of `Cache` propagates a `missing` solution vector
    # into the solver constructor, where it fails with an error that gives no hint as to why.
    x = nlsolution(caches)
    ismissing(x) && throw(ArgumentError(string(
        nameof(typeof(method)), " does not support problems of type ",
        nameof(typeof(equation(caches.problem))), ".")))

    NonlinearSolver(solvermethod, zero(x), residual!, zero(x); kwargs...)
end

# This accounts for the SimpleSolvers interface, expecting a single parameter argument,
# whereas the typical `residual!` methods expect a number of additional arguments.
residual!(y, x, parameters::Union{Tuple,NamedTuple}) = residual!(y, x, parameters...)
