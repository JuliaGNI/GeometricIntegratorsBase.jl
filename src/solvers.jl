
default_linesearch(method=nothing) = Backtracking()

# `Backtracking()` is always `Backtracking{Float64}`, so a `Float32` integration makes the
# `Linesearch` constructor do a `change_precision` the caller could have avoided by building the
# method at the working type in the first place. Both methods are kept: the untyped one is the
# downstream hook (`GeometricIntegrators` calls it for DIRK) and its signature must not change.
#
# `expand = true` is deliberately not set. The expansion phase of `Backtracking` is opt-in
# upstream because it costs a merit evaluation — a full residual evaluation here — on a direction
# whose scale is wrong, and gains nothing on one whose scale is right: a Newton step already sits
# at its model minimum. Every method in this package solves with `Newton()`.
default_linesearch(::Type{T}, method=nothing) where {T} = Backtracking(T)

default_options(method::GeometricMethod, problem::GeometricProblem) = (
    min_iterations=1,
    f_abstol=max(8, solversize(method, problem)) * eps(datatype(problem)),
    # Bound a solve whose residual sits on a floor above `f_abstol` while its iterate keeps
    # moving normally: no convergence branch fires, `max_stalls` cannot see it (the step is not
    # small), and the solve spends `max_iterations` — 1000 by default — on *every* time step.
    #
    # 50 is the conservative point. An iteration converging linearly with rate ρ improves by ρ^W
    # over a window W, so at SimpleSolvers' `f_stall_factor = 0.5` a window of 50 gives up only on
    # ρ > 2^(-1/50) ≈ 0.986, i.e. on a residual that is essentially flat — which no Newton solve
    # on these residuals is. It also stays above `SimpleSolvers.F_STALL_REPORT_MINIMUM`, below
    # which SimpleSolvers itself holds the no-progress proportion to be no evidence of anything;
    # giving up on less than the package will report on would be incoherent.
    #
    # Giving up and converging remain mutually exclusive (`no_progress` is gated on
    # `!residual_small`), the solve still reports the residual it achieved against the tolerance
    # it was asked for, and `default_options` is merged rather than replaced, so a caller who
    # disagrees overrides it with one keyword.
    f_stall_window=50,
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
