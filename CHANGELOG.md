# Release Notes

All notable changes to GeometricIntegratorsBase.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually changed,
so that a compat-only bump can be told apart from an interface change.


## 0.6.3

### Breaking Changes

* Requires SimpleSolvers 0.12.1. 0.12 is breaking for a caller of this package in one way that
  matters: a `NonlinearSolver` no longer emits line-search warnings from inside its iteration, so a
  rejected line search reaches its caller through the returned status and nobody at all if the
  caller drops it. That is what the change below is for. The other two breaks of 0.12 —
  `NonlinearSolverStatus` gaining a field, and a `LinesearchMethod` implementing
  `solve_with_status` rather than `solve` — do not reach this package, which constructs neither.

  Nothing measurable moves. 0.12's other behavioural change, `Bisection` bisecting toward a minimum
  rather than toward whichever root the value bracket held, does not apply here either:
  `default_linesearch` is `Backtracking` and every method in this package solves with `Newton()`.
  The full suite passes unedited, including the two assertions that would have caught a change in
  what a solve says — `test/integrators/common_tests.jl`, which requires a converged solve to be
  *silent*, and `test/integrator_tests.jl`, which pins the exact wording of the failure warning.

### New Features

* Every nonlinear solve now goes through `solve_with_status!` rather than `solve!`, and the
  `SimpleSolvers.NonlinearSolverStatus` it returns is handed to a new `check_solver_status`. This
  is the seven `integrate_step!` methods of `ImplicitEuler`, `ImplicitMidpoint` and
  `CrankNicolson`, and it replaces the pair of commented-out stubs (`# println(status(solver))`,
  `# println(meets_stopping_criteria(status(solver)))`) that had sat under the solve since before
  the status existed in its current form.

  `solve_with_status!` is used in its state-taking form, added in SimpleSolvers 0.12.1 at this
  package's request. The other two forms build a `NonlinearSolverState` per call, which in a
  time-stepping loop is one allocation per step and would also have detached the status from
  `solverstate(int)` — the object a caller reads a solve through afterwards, as
  GeometricIntegrators' SPARK tests do.

  **The migration costs nothing.** That is worth stating because it was not true of the first cut:
  `solve_with_status!` was `solve!` followed by `status(s, state)`, and `solve!` had *already* built
  a status internally to hand to its own warning path and then discarded it — so reading the outcome
  cost a second `NonlinearSolverStatus`, i.e. five `l2norm` passes per solve, per step. The review of
  the upstream pull request caught it, and 0.12.1 ships both entry points sharing one body that
  returns the status the loop built. So the status is now genuinely free: it is the value
  SimpleSolvers computes for its own report either way, handed back instead of thrown away.

* `check_solver_status(status, int)` is **silent by default, deliberately**, and the docstring says
  why at length: SimpleSolvers reports a solve that did not converge itself, once per solve and
  with a back-off (occurrences 1, 2, 4, 8, …) that exists precisely because a time-stepping loop is
  the case that floods. Warning here as well would say the same thing twice per step. What the
  function buys is that the status is never silently dropped and that there is *one* place to
  change the policy for all of them — overriding it is how a downstream family of methods becomes
  strict about non-convergence, in the same way that overriding `default_options` is how it changes
  solver options in one place. It is exported for that reason.

* The line-search tally is covered by `test/solver_tests.jl`. It is the one thing the migration
  exists for and the one thing nothing asserted: the suite pinned the wording of the
  non-convergence warning and required a converged solve to be silent, but a *rejected line search*
  — the diagnostic 0.12 stopped logging from inside the iteration — was checked by nothing, which
  is why widening the compat entry alone passed CI green. The new assertion drives a solve onto an
  ascent direction with a sign-flipped Jacobian and requires the rejection to be readable off the
  status that reaches `check_solver_status`.

  It asks `dominant_linesearch_outcome(status, false)` rather than `linesearch_failures(status) > 0`
  because the latter is vacuous here: a converged solve reaches the merit's round-off floor on its
  last step, `LINESEARCH_FLOOR` counts as a failure, and the bare count is therefore positive for a
  solve where nothing went wrong. The converged solve is asserted the same question alongside it,
  so the two readings differ for the reason claimed.


## 0.6.2

### Breaking Changes

* `ExplicitEuler` and `ImplicitEuler` dispatch on `ProblemODE` rather than `AbstractProblemODE`.
  The latter includes `DAEProblem`, so `integrate(daeproblem(), ExplicitEuler())` used to return
  a solution that silently ignored the constraint. `DAEProblem` is the only type removed.
* A fallback `integrate_step!` throws the `ArgumentError` of `initsolver` verbatim, so every
  unsupported method/problem combination is now rejected in the same words. Previously the
  rejection surfaced as whatever exception the first unimplemented piece of the machinery
  happened to raise — a `MethodError` for the symplectic Euler methods, a `FieldError` for
  `ImplicitEuler` on a `PDAEProblem`, an `ArgumentError` for the implicit ODE methods.

  Shipped as a patch rather than as 0.7.0; this was a maintainer's call recorded in the commit.

### New Features

* `ImplicitMidpoint` and `CrankNicolson` for `PODEProblem` and `HODEProblem`, through the same
  structs and the same traits mechanism that already carried their IODE/LODE support. For
  q̇ = v(t,q,p), ṗ = f(t,q,p) the solver solution vector holds both stage vector fields, so the
  nonlinear system is twice the size of the one for an ordinary differential equation.
* A new `ProblemPODE` union alongside `ProblemODE` and `ProblemIODE`, so that `PDAEProblem` and
  `HDAEProblem` keep being rejected rather than integrated as if their constraint were absent.
  `symplectic_euler.jl` is switched to it as well.

### Tests

* The test suite is reorganised: test problems to `test/examples/`, per-method test files to
  `test/integrators/`, mirroring `src/integrators/`. `euler_tests.jl` becomes
  `explicit_euler_tests.jl`.
* The five method files repeated the same tests almost verbatim. They now live in
  `test/integrators/common_tests.jl`, driven by a table of every method with the problem types
  it implements, the order it converges at and the accuracy it reaches. Coverage grows from 373
  to 558 assertions: unsupported-problem rejection now covers 51 method/problem combinations
  instead of 19.
* Both partitioned test problems were linear *and* separable, so the off-diagonal blocks of the
  residual Jacobian vanished and the two halves of the solver solution vector were never
  actually coupled. A nonlinear, non-separable partitioned problem and its equivalent first
  order system are added, so the equivalence tests pin the coupled residual down.


## 0.6.1

### Documentation

* `default_options` gains a docstring covering the merge semantics and the reasoning behind each
  of its defaults, and `src/solvers.jl` is rendered in the manual on a new Solvers page. Every
  `@autodocs` block here is `Pages`-filtered and `solvers.jl` was in none of them, so a dangling
  reference in these docstrings used to fail GeometricIntegrators' documentation build rather
  than this one.


## 0.6.0

### Breaking Changes

* Requires SimpleSolvers 0.11. Nothing this package used was removed — the one breaking change
  upstream is `Backtracking` losing its `α₀` key, the `DEFAULT_ARMIJO_α₀` constant and the
  positional constructor, none of which appear here.
* SimpleSolvers 0.11 tests `all(isfinite, direction)` where 0.10 tested only `isnan`, so an
  *overflowed* direction — which used to pass every guard and stall silently, `Inf * nan_factor`
  being `Inf` — now raises `NonlinearSolverException`. The time-stepping loop catches it, warns
  naming the timestep, and breaks, as the adjacent `isnan` check does. The `break` precedes the
  `copy!`, so the solution holds valid data up to `n-1`; the warning says where the valid data
  ends, since the entries past it are zeros rather than anything a caller can recognise. Only
  that one exception type is caught; anything else is rethrown.

### New Features

* `f_stall_window = 50` in `default_options`, SimpleSolvers 0.11's criterion for the case
  `max_stalls` cannot see: the residual sitting on a floor above `f_abstol` while the iterate
  keeps moving, so that no criterion fires and the solve spends `max_iterations` — 1000 by
  default — on every single time step. 50 is conservative: at `f_stall_factor = 0.5` it gives up
  only on a linear rate ρ > 2^(-1/50) ≈ 0.986, and it stays above
  `SimpleSolvers.F_STALL_REPORT_MINIMUM`. It is justified for `Newton()`; a fixed-point
  iteration converging at a rate that close to the limit would now be cut off, which is recorded
  as a limitation. Options are merged, so a caller overrides it with one keyword.
* `default_linesearch(::Type{T}, method)` for `T <: Real`, alongside the untyped hook that
  GeometricIntegrators calls for DIRK. `Backtracking()` is always `Backtracking{Float64}`, so a
  `Float32` integration made the `Linesearch` constructor do a `change_precision` the caller
  could have avoided. `Real` rather than `AbstractFloat`, so that `ForwardDiff.Dual` does not
  silently fall through to the untyped hook.

The `Backtracking` expansion phase added in SimpleSolvers 0.11 is deliberately *not* enabled: it
costs a full residual evaluation on a direction whose scale is wrong and gains nothing on one
whose scale is right, a Newton step already sitting at its model minimum.


## 0.5.3

### Breaking Changes

* `ImplicitMidpoint` and `CrankNicolson` dispatch on unions of the unconstrained problem types.
  The problem unions of GeometricEquations also cover the constrained variants, so dispatching
  on `AbstractProblemIODE` let `IDAEProblem`s and `LDAEProblem`s through, and they were
  integrated as if the constraint were absent — multipliers left at zero, ϕ, u and g never
  evaluated. The same held for `DAEProblem`s on the ODE side.

### New Features

* IODE/LODE versions of `ImplicitMidpoint` and `CrankNicolson`. The `isiodemethod` and
  `islodemethod` traits are set accordingly, which cannot follow from the supertype.

### Fixes

* Problem types a method does not implement at all used to run into the generic `NoCache`
  fallback and propagate a `missing` solution vector into the solver constructor, failing with
  an error that gave no hint as to why. `initsolver` now checks for this and names both the
  method and the problem type.


## 0.5.2

### New Features

* `SymplecticEulerA` and `SymplecticEulerB`, partitioned methods for separable Hamiltonians
  H(t,q,p) = T(t,p) + V(t,q). Under that assumption the otherwise implicit scheme decouples into
  two explicit substeps, so neither a nonlinear solver nor an initial guess is required. Both
  cover PODE and HODE problems.
* `ImplicitMidpoint` and `CrankNicolson` for `AbstractProblemODE`, following `ImplicitEuler` in
  structure. The nonlinear unknown is the stage vector field in both cases. The Crank-Nicolson
  cache holds v(t̄, q̄) at working precision, so that it survives the automatic differentiation
  of the residual, where it is a constant.

  All three reproduce the corresponding tableau-driven methods of GeometricIntegrators:
  symplectic Euler A/B bit for bit on PODE, HODE and non-autonomous problems, implicit midpoint
  and Crank-Nicolson to round-off.

* All method implementations are collected in a flat `src/integrators/` directory.

### Fixes

* `ExplicitEuler` evaluated the vector field at the *end* of the step. On entry to
  `integrate_step!`, `sol.t` already holds tₙ₊₁ while `sol.q` still holds qₙ, so
  `v(sol.t, sol.q)` gave v(tₙ₊₁, qₙ). The scheme stays first order consistent, which is why no
  test caught it, and the harmonic oscillator is autonomous, so those results are unchanged — but
  on a non-autonomous problem the method was simply not explicit Euler.
* `compute_energy_error` in the harmonic oscillator test problem referenced `DataSeries.nt`, a
  field that no longer exists, so any call threw. Two further bugs were hidden behind it:
  `axes(q, 2)` returns `OneTo(1)` for a one-dimensional `DataSeries`, so the loop covered a
  single index rather than all time steps, and `q[:, i]` extracts one component across all times
  instead of the state at time `i`. It now delegates to the `compute_invariant` and
  `compute_relative_error` utilities of GeometricSolutions.

### Tests

* None of the committed test problems was non-autonomous, so nothing exercised the time
  arguments of the integrators: swapping t̄ for `sol.t` in any of them left the entire suite
  green. `test/nonautonomous.jl` adds q̇ = cos(t) q, whose exact solution is known, and the
  separable Hamiltonian H(t,q,p) = p²/2 + (1+t) q²/2 as PODE and HODE.


## 0.5.1

### Fixes

* Restore a minimum `f_abstol` of `8eps(T)`. With `max(8, solversize(method, problem))` the
  tolerance now resolves to 1.8e-15 … 3.6e-15 for the methods here, above their residual floors.


## 0.5.0

### Breaking Changes

* `solversize` and `default_options` take their arguments in a uniform `(method, problem)` order,
  so that all methods dispatch the same way and callers no longer have to remember a per-function
  convention. The `ImplicitEuler` specialization of `solversize` and the generic fallback in
  `method.jl` are flipped accordingly, and `initmethod`'s two-argument fallback is widened to
  `::GeometricProblem` for consistency.
* `default_options` sizes `f_abstol` from `solversize(method, problem)`, which is what the new
  argument order is for.

### Fixes

* `maximum(2, ...)` threw a `MethodError`, treating `2` as a function passed to `maximum(f, itr)`;
  corrected to `max(2, ...)`.


## 0.4.3

### Breaking Changes

* Requires SimpleSolvers 0.10, which removes the root cause of the solver and line-search warning
  floods downstream packages had built log-suppression scaffolding around: a line search now
  shares its solver's `Options` (so `verbosity = 0` reaches it), converged-at-floor solves are
  silent, and a stalled solve stops after `max_stalls` instead of spinning to `max_iterations`.
* **`default_options` is merged with the caller's options rather than replaced wholesale.**
  Passing *any* option used to discard the entire default set — a documented trap downstream, and
  one that changes results rather than just performance, since silently losing `min_iterations = 1`
  controls whether a step is taken at all now that SimpleSolvers tests its stopping criteria
  before the first step.
* `f_abstol` is dropped from `default_options`. It is an *absolute* bound on ‖F(x)‖, so it has to
  sit above the round-off floor of the residual — which is the cancellation level of the terms
  `residual!` sums, not a property of the solver. `8eps() ≈ 1.78e-15` is below that floor for
  essentially every non-trivial problem, so the criterion was unreachable and Newton ran to its
  iteration limit. `f_reltol` does not rescue this: it is anchored at the *initial* residual, so
  the Hermite-extrapolation initial guess used throughout these integrators — good enough that
  Newton converges in one iteration — makes the relative gate tighter, not looser.

  (Partly reverted in 0.5.0/0.5.1, which reintroduce a *size-scaled* `f_abstol` with an `8eps(T)`
  floor rather than a flat one.)

### Tests

* The tolerances this made unnecessary are cleaned up: the below-floor and no-longer-load-bearing
  `f_abstol`/`max_iterations` values are dropped, the now-redundant `min_iterations = 1`
  restatements removed, and converged solves assert silence via `@test_nowarn`, which is usable
  as a regression tripwire again now that a converging solve is silent.


## 0.4.2

### New Features

* `NormalizedHermiteExtrapolation`, a version of the Hermite extrapolation normalised to the
  interval [0,1]. Instead of the sample times t₀ and t₁ and the time tᵢ to extrapolate to, it
  takes a normalised time cᵢ with tᵢ = t₁ + cᵢ Δt, so that cᵢ = -1 reproduces the first and
  cᵢ = 0 the second sample; all derivative values are with respect to the normalised time. Both
  versions share a single kernel; `HermiteExtrapolation` is unchanged in signature, behaviour and
  allocations.

### Fixes

* `NormalizedHermiteExtrapolation`'s `solutionstep!` methods took Δt from `timestep(problem)`
  while taking t₁ from the history. Whenever the spacing of the history differs from the nominal
  timestep the resulting polynomial was silently wrong, and the `t₀ == t₁` guard of
  `HermiteExtrapolation` was lost. Δt = t₁ - t₀ is now taken from the history, so that both
  versions agree unconditionally.


## 0.4.1

* Removes PProf and ProfileView from the explicit test dependencies.


## 0.4.0

### Breaking Changes

* `SolutionStep`'s `reset!` method takes a time rather than a timestep.


## 0.3.3

### Fixes

* Fix in midpoint extrapolation.


## 0.3.2

* Use the default linesearch as prescribed by SimpleSolvers, rather than naming one here.


## 0.3.1

### Breaking Changes

* Requires SimpleSolvers 0.9.


## 0.3.0

### Breaking Changes

* The dimension type parameter is removed from `Cache`, and the cache constructors are
  generalised.
* The type parameter order of `IntegratorCache` is changed to comply with the order of its
  subtypes.
* `NoSolver` is moved to GeometricBase.
* The `Newton` alias is removed; it is not available in SimpleSolvers.
* The `Parameters` dependency is dropped.
* The `GeometricIntegrator` fallback constructor returns `missing` rather than `nothing`.
* Two cases of type piracy are removed: `ndims` in the main module, and `prince_reference` for
  `GeometricMethod`, which was wrapped in a `try`/`catch` that swallowed unrelated failures.
* Assertions and generic `@error`s are replaced by specific `ArgumentError`s throughout, so a
  misuse now fails with a message that names what was wrong.

### Fixes

* Fixes in Euler and midpoint extrapolation, and a simplification of the Hermite interpolation.
* `update!(solstep, ...)` no longer leaves a partially mutated state behind when it fails.
* Bugfixes in `SolutionStep`.
* A fallback for `reference(::GeometricMethod)` and a missing `isdelemethod`.
* Dead code and redundant `local` keywords removed; `residual!` cleaned up.


## 0.2.2

### Fixes

* Fix solver initialisation.


## 0.2.1

* Cleanup in tests.


## 0.2.0

### Breaking Changes

* `SolutionStep` is revised, both for correctness and to fix performance problems.
* Requires GeometricBase 0.14, GeometricEquations 0.21 and GeometricSolutions 0.6.1.

### Tests

* Tests for extrapolation of velocity and force, and for Hermite extrapolation of PODE and IODE
  problems.


## 0.1.11

### New Features

* `initsolver` is generalised to arbitrary nonlinear solvers rather than just Newton.


## 0.1.10

### Breaking Changes

* Requires Julia 1.10 and SimpleSolvers 0.7.5.


## 0.1.9

### Breaking Changes

* Requires SimpleSolvers 0.7.


## 0.1.8

### New Features

* NaN detection in the main `integrate` method: the time-stepping loop stops and warns rather
  than filling the rest of the solution with garbage.


## 0.1.7

### Breaking Changes

* Requires SimpleSolvers 0.6 and Julia 1.8.

### New Features

* A custom time step for the midpoint extrapolation of PODE problems.


## 0.1.6

### Fixes

* Fix in Hermite extrapolation.


## 0.1.5

### Fixes

* Fix midpoint extrapolation.


## 0.1.4

### New Features

* `NoInitialGuess` and `NoExtrapolation` types.
* A time step limiter for the midpoint extrapolation.

### Breaking Changes

* The problem-specific solution type unions move to GeometricSolutions.


## 0.1.3

### Fixes

* Minor bugfixes and extensions around the extrapolation methods and the initialisation of
  `SolutionStep`.


## 0.1.2

### Fixes

* Import `integrate!` from GeometricBase, so that it is extended rather than shadowed.


## 0.1.1

### Fixes

* Cleanup and fixes in the `enforce_periodicity!` methods.
* Fix the `copy!` calls in `set_property!` for `SolutionStep`.
* Several further fixes in `SolutionStep`.


## 0.1.0

Initial release. The integrator framework extracted from
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl), so that packages
implementing their own methods can depend on the framework without pulling in the full method
collection:

* `SolutionStep` and the solution machinery.
* `GeometricMethod`, `GeometricIntegrator` and the trait interface.
* `Cache` and `CacheDict`.
* The extrapolation methods: `EulerExtrapolation`, `MidpointExtrapolation` and
  `HermiteExtrapolation`, with `NoExtrapolation` and `NoInitialGuess`.
* The explicit and implicit Euler integrators.


## Open Issues

The items below are the unresolved entries of `todo.md`, kept here so that they are visible
alongside the releases that created them.

### `check_solver_status` acts on nothing — open, by decision

0.6.3 routes every step's solve through `solve_with_status!` and hands the status to
`check_solver_status`, whose default body returns it and does nothing else. That was the deliberate
choice for the release — SimpleSolvers is left as the single reporting voice, so no existing run
changes what it prints — but it means the status is *available* rather than *used*, and the
interesting thing to do with it is still undone.

The obvious next step is one level up, in `integrate!`: it already recognises two ways a step can
go wrong (a `NonlinearSolverException`, and NaNs in the iterate), warns naming the time step, and
returns the trajectory computed so far rather than discarding it. A step whose solve merely failed
to converge is a third, is now detectable for the first time, and currently produces a trajectory
that continues past the point where it stopped being trustworthy with nothing in `sol` to mark it.
Doing this is a behaviour change for any run that currently limps along, which is why it was not
folded into a compat bump.

### A repeating non-convergence is reported on 1, 2, 4, 8, … — open

SimpleSolvers 0.12 replaced the `maxlog` caps on its solver report with a back-off, so a diagnosis
that repeats is reported on its 1st, 2nd, 4th, 8th … occurrence. In a time-stepping loop that is
the right shape — the alternative is one message per step — but it means occurrence 10 of a
failing solve is silent, and nothing in this package compensates or counts. A caller who wants
"how many of my 10 000 steps did not converge" cannot get it from the log, and `check_solver_status`
is where a tally would go.

### `default_options` restatements downstream — open

Several callers compensate for the pre-0.4.3 replace-not-merge behaviour by restating defaults
they did not want to change, `min_iterations = 1` in particular. Since 0.4.3 merges, those
restatements are redundant rather than wrong, so removing them is safe cleanup and not a break.
They have not been audited.

### `max_iterations` bounds downstream — open

Repositories that bounded a non-converging solve by lowering `max_iterations` no longer need to.
Since 0.6.0, `f_stall_window` bounds such a solve without also bounding one that is making
progress, which is what a low `max_iterations` could not distinguish. Worth revisiting in
GeometricIntegrators, GeometricProblems and ChargedParticleDynamics.

### Whether to reinstate `f_abstol` — superseded, no action

Recorded because the question is a natural one to re-open. 0.4.3 dropped `f_abstol` on the
grounds that an absolute residual bound below the round-off floor is unsatisfiable; the follow-up
question was whether SimpleSolvers 0.10's stall reporting made a value choosable by reading one
warning rather than by bisecting tolerances. It does not settle the case: `max_stalls` only
covers a solve whose *iterate* has frozen, while the case at issue — a residual floor above the
tolerance, with the iterate still moving — is the one `f_stall_window` was added for in 0.11.0
and adopted here in 0.6.0. 0.5.0/0.5.1 settled the tolerance itself, by sizing it from
`solversize` with an `8eps(T)` floor. Do not re-derive this from the 0.4.3 entry alone.
