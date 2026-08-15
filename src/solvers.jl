
"""
    default_linesearch([T,] method=nothing)

The line search a `method` uses by default, built at the working type `T` where one is given.

`Backtracking()` is always `Backtracking{Float64}`, so a `Float32` integration makes the
`Linesearch` constructor do a `change_precision` the caller could have avoided by building the
method at the working type in the first place. Both forms exist: the untyped one is the downstream
hook (`GeometricIntegrators` calls it for DIRK) and its signature must not change.

Neither form is called inside this package. They are the hook downstream methods reach for, which
is why the tests assert both signatures and the choice below.

`T` is constrained so that the typed method cannot quietly intercept a call meant for the untyped
one: a caller passing a *method type* rather than an instance — `default_linesearch(ImplicitEuler)`
— falls through to the untyped hook and gets a `Backtracking{Float64}`, as it did before the typed
method existed, rather than an error from inside `Backtracking`. `Real` rather than `AbstractFloat`,
because `Backtracking` is buildable at any real type that `T(::Float64)` accepts, which includes
the `ForwardDiff.Dual` a caller differentiating through an integration would arrive with.

`expand = true` is deliberately not set. The expansion phase of `Backtracking` is opt-in upstream
because it costs a merit evaluation — a full residual evaluation here — on a direction whose scale
is wrong, and gains nothing on one whose scale is right: a Newton step already sits at its model
minimum. Every method in this package solves with `Newton()`.
"""
default_linesearch(method=nothing) = Backtracking()
default_linesearch(::Type{T}, method=nothing) where {T<:Real} = Backtracking(T)

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

"""
    default_options(method, problem)

The options a `method` is solved with unless the caller says otherwise. They reach the nonlinear
solver as the keyword arguments of `initsolver`.

`GeometricIntegrator` *merges* them with the options it is handed —
`merge(default_options(method, problem), options)` in `src/integrator.jl` — so a caller who sets
one of them keeps the rest. Replacing the whole bundle instead, as this package did up to 0.4.x,
meant that passing any option at all silently dropped every default not restated alongside it.

The bundle is queried on the `initmethod`-specialised method, so `solversize` below sees the
concrete integrator rather than the method the caller named: `LobattoIIIA(2)` on a 2-dof ODE
arrives here as a `DIRK` with `solversize == 2`, `RadauIIA(2)` as an `IRK` with `solversize == 4`.

* `min_iterations = 1` — SimpleSolvers tests its stopping criteria *before* the first step, so
  without this a solve whose initial guess already satisfies them takes no step at all, and the
  extrapolated initial guesses used throughout this package are regularly that good.
* `f_abstol` scales with the size of the solver system: `max(8, solversize(method, problem))`
  residual components at `eps(datatype(problem))` apiece, because the norm of a residual sitting
  on its round-off floor grows with the number of components it sums. `f_reltol` is no substitute
  — it is anchored at the *initial* residual, so a good initial guess makes the relative gate
  tighter rather than looser. The floor is 8 and not 2 because the smallest solver systems are
  genuinely tiny: at `2eps ≈ 4.4e-16` the `DIRK` case above sits below its own round-off floor
  and stagnates, whereas `8eps ≈ 1.78e-15` clears every floor measured downstream and is what
  this package defaulted to before the size factor existed.
* `f_stall_window` — see [`DEFAULT_F_STALL_WINDOW`](@ref).

Overriding this is how a downstream family of methods changes solver options in one place rather
than at every call site. The two overrides that existed for `f_abstol` — one for the implicit
Runge-Kutta families, one for SPARK — were both removed once the size scaling landed, since the
sized value lands above the residual floor of each of them.
"""
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

"""
    check_solver_status(status, int)

Act on the outcome of the nonlinear solve of one time step. Every `integrate_step!` that solves
calls this with the `SimpleSolvers.NonlinearSolverStatus` that `solve_with_status!` returned, and
returns whatever it returns.

**It is silent by default, and that is the point of it rather than an omission.** SimpleSolvers
reports a solve that did not converge itself, from the end of its own `solve!`, and since 0.12 the
status is the *programmatic* counterpart of that report rather than a replacement for it — a caller
that wants to act on a rejected line search reads `SimpleSolvers.dominant_linesearch_outcome`
instead of scraping the log. Warning here as well would say the same thing twice per time step, and
SimpleSolvers' own report already backs off (occurrences 1, 2, 4, 8, …) precisely because a
time-stepping loop is the case that floods.

What it buys is that the status is never silently dropped: it arrives somewhere named, on every one
of the call sites, and there is one place — this one — to change the policy for all of them.
Overriding it is how a downstream family of methods becomes strict about non-convergence, in the
same way that overriding [`default_options`](@ref) is how it changes solver options in one place:

```julia
function GeometricIntegratorsBase.check_solver_status(status, int::GeometricIntegrator{<:MyMethod})
    SimpleSolvers.isconverged(status) || throw(NonlinearSolverException("step did not converge"))
    status
end
```

`NonlinearSolverException` is the type to throw, because `integrate!` already catches exactly that
one, warns with the time step it failed at, and returns the trajectory computed so far instead of
discarding it — see `src/integrate.jl`. Anything else propagates and loses the solution.

Both arguments of the fallback are untyped on purpose. `status` is, so that this package does not
have to enumerate SimpleSolvers' status types to stay reachable — every solver method produces a
`NonlinearSolverStatus` today, but the hook has no reason to care. `int` is, so that *every*
integrator reaches the default; an override narrows one or both, which is what makes it a per-method
policy rather than a global one.
"""
check_solver_status(status, int) = status
