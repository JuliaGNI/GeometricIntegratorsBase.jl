# Known issues

What is known to be wrong in GeometricIntegratorsBase and is not fixed. An entry leaves this file
when its fix merges, and the CHANGELOG entry of the fix names its ID.

### K1 · The manual re-renders its dependencies' docstrings — open

- **location:** `docs/src/deps/equations.md`
- **evidence:**

  `docs/src/deps/equations.md` and `deps/problems.md` render 34 GeometricEquations docstrings each,
  through explicit `@docs` blocks. GeometricEquations links within itself using plain `@ref`, which is
  correct in its own manual and dangling in this one — so **any upstream release that documents a
  newly referenced symbol breaks this build with no commit here**, and it breaks on `main` rather
  than on a pull request, because the pull request resolves the old version and the merge the new one.
  That is what cost the 0.6.5 manual.

  The fix is not a change of link syntax. `@extref` resolves a link written in *this* manual's own
  prose; it cannot reach inside an upstream docstring that an `@docs` block has rendered here. So
  removing the failure class means dropping those `@docs` blocks and linking out to the dependency's
  published manual instead — which is what the `InterLinks` plugin in `docs/make.jl` is already
  configured for, and what GeometricSolutions' own docstrings already do from the other side. That is
  a decision about what this manual should contain rather than a build fix, which is why it was not
  folded in alongside the fix above.
- **kind:** docs
- **found:** 2026-09-02

### K2 · No method declares its full property set, so `isAbstractMethod` is `false` for all of them — open

- **location:** —
- **evidence:**

  `GeometricBase.isAbstractMethod` requires all ten method properties to be non-`missing`. Since the
  0.6.7 fix the six `is*` predicates are at least *askable* through `GeometricBase`, but
  `isenergypreserving`, `isstifflyaccurate`, `order`, `name`, `description` and `reference` are
  undefined for every method in this package, so the predicate still answers `false` for all of them
  and cannot be used as the interface conformance check it is meant to be. Filling those in is per
  method, not a single change, and each one is a claim about the scheme that has to be right.
- **kind:** defect
- **found:** 2026-09-02

### K3 · `print_reference` is two generics, and neither prints for a Runge-Kutta method — open

- **location:** —
- **evidence:**

  This package defines `print_reference` and `GeometricIntegrators` defines its own without importing
  it, so they are separate functions; the downstream one is what actually prints, via
  `reference(tableau(method))`. Meanwhile `RungeKutta` attaches its reference strings to tableau types
  and nothing attaches one to a method wrapper, so this package's `print_reference` prints nothing for
  `Gauss`, `VPRK` and the rest. Unifying the two and giving the method families a `reference` belongs
  in `GeometricIntegrators`, which is why the 0.6.7 fix stopped at making this package's
  version read the shared generic.
- **kind:** defect
- **found:** 2026-09-02

### K4 · `issymmetric` is exported here and means something else in `LinearAlgebra` — open

- **location:** —
- **evidence:**

  The one collision of this kind that importing cannot resolve, because the two upstream functions are
  genuinely different. `LinearAlgebra` is a direct dependency and exports its own `issymmetric`, a
  predicate on matrices; `GeometricBase` declares an unexported `issymmetric`, a property of an
  integration method. This package imports the latter and exports it, so `using LinearAlgebra` and
  `using GeometricIntegratorsBase` together still leave a caller with two bindings for one name and no
  way to pick between them but to qualify. Pre-existing — the shadowing fix neither caused it nor made
  it worse — but it is the same failure class, and it is what the new interface guard surfaces, so it
  is recorded rather than left implicit. Resolving it means renaming the method property upstream in
  `GeometricBase`, which is not this package's call.
- **kind:** upstream
- **found:** 2026-09-02

### K5 · `check_solver_status` acts on nothing — open, by decision

- **location:** —
- **evidence:**

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
- **kind:** defect
- **found:** 2026-08-15; the same issue as GeometricIntegrators' entry *The solver status is
  available but not acted on*, which names this package as the place to act

### K6 · A repeating non-convergence is reported on 1, 2, 4, 8, … — open

- **location:** —
- **evidence:**

  SimpleSolvers 0.12 replaced the `maxlog` caps on its solver report with a back-off, so a diagnosis
  that repeats is reported on its 1st, 2nd, 4th, 8th … occurrence. In a time-stepping loop that is
  the right shape — the alternative is one message per step — but it means occurrence 10 of a
  failing solve is silent, and nothing in this package compensates or counts. A caller who wants
  "how many of my 10 000 steps did not converge" cannot get it from the log, and `check_solver_status`
  is where a tally would go.
- **kind:** upstream
- **found:** 2026-08-15

### K7 · `default_options` restatements downstream — open

- **location:** —
- **evidence:**

  Several callers compensate for the pre-0.4.3 replace-not-merge behaviour by restating defaults
  they did not want to change, `min_iterations = 1` in particular. Since 0.4.3 merges, those
  restatements are redundant rather than wrong, so removing them is safe cleanup and not a break.
  They have not been audited.
- **kind:** not verified
- **found:** 2026-08-15

### K8 · `max_iterations` bounds downstream — open

- **location:** —
- **evidence:**

  Repositories that bounded a non-converging solve by lowering `max_iterations` no longer need to.
  Since 0.6.0, `f_stall_window` bounds such a solve without also bounding one that is making
  progress, which is what a low `max_iterations` could not distinguish. Worth revisiting in
  GeometricIntegrators, GeometricProblems and ChargedParticleDynamics.
- **kind:** not verified
- **found:** 2026-08-15
