
"""
GeometricIntegrator

Collects all data structures needed by an integrator:

* `problem`: [`EquationProblem`](@ref) to solve
* `method`: integration method
* `cache`: temprary data structures needed by method
* `solver`: linear or nonlinear solver needed by method
* `iguess`: initial guess for implicit methods
* `projection`: optional projection method

Constructors:

```
GeometricIntegrator(problem::EquationProblem, method::GeometricMethod; solver = default_solver(method), iguess = default_iguess(method), projection = default_projection(method))
```

"""
struct GeometricIntegrator{
    MT <: GeometricMethod,
    PT <: AbstractProblem,
    CT <: CacheDict{PT, MT},
    ST <: Union{<:AbstractSolver, NoSolver},
    IT <: Union{InitialGuess, Extrapolation},
    SST <: AbstractSolverState
} <: AbstractIntegrator
    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
    solverstate::SST
end

function GeometricIntegrator(
        problem::AbstractProblem,
        integratormethod::GeometricMethod,
        solvermethod::SolverMethod,
        iguess::Union{InitialGuess, Extrapolation};
        method = initmethod(integratormethod, problem),
        caches = CacheDict(problem, method),
        options...
)
    solver = initsolver(solvermethod, method, caches; merge(default_options(method, problem), options)...)
    GeometricIntegrator(problem, method, caches, solver, iguess, SolverState(solver))
end

function GeometricIntegrator(
        problem::AbstractProblem,
        integratormethod::GeometricMethod;
        method = initmethod(integratormethod, problem),
        solver = default_solver(method),
        initialguess = default_iguess(method),
        kwargs...
)
    GeometricIntegrator(problem, method, solver, initialguess; kwargs...)
end

GeometricIntegrator(::AbstractProblem, ::Nothing, args...; kwargs...) = missing

problem(int::GeometricIntegrator) = int.problem
method(int::GeometricIntegrator) = int.method
caches(int::GeometricIntegrator) = int.caches
solver(int::GeometricIntegrator) = int.solver
iguess(int::GeometricIntegrator) = int.iguess
initialguess(int::GeometricIntegrator) = int.iguess
solverstate(int::GeometricIntegrator) = int.solverstate

cache(int::GeometricIntegrator, DT) = caches(int)[DT]
cache(int::GeometricIntegrator) = cache(int, datatype(problem(int)))
hasnullvector(int::GeometricIntegrator) = hasnullvector(method(int))
implicit_update(int::GeometricIntegrator) = implicit_update(method(int))
nconstraints(int::GeometricIntegrator) = nconstraints(problem(int))
Base.ndims(int::GeometricIntegrator) = ndims(problem(int))
nlsolution(int::GeometricIntegrator) = nlsolution(cache(int))
nullvector(int::GeometricIntegrator) = nullvector(method(int))
tableau(int::GeometricIntegrator) = tableau(method(int))

equations(int::GeometricIntegrator) = functions(problem(int))
timestep(int::GeometricIntegrator) = timestep(problem(int))

initial_guess!(sol, history, params, ::GeometricIntegrator) = nothing

function integrate(problem::AbstractProblem, method::GeometricMethod; kwargs...)
    integrator = GeometricIntegrator(problem, method; kwargs...)
    integrate(integrator)
end

function integrate(problems::EnsembleProblem, method::GeometricMethod; kwargs...)
    solutions = Solution(problems)

    for (problem, solution) in zip(problems, solutions)
        integrator = GeometricIntegrator(problem, method; kwargs...)
        integrate!(solution, integrator)
    end

    return solutions
end

"""
Performs one integration step of an integrator.
"""
function integrate_step! end

# A method that implements no `integrate_step!` for the problem at hand does not support it.
# Without this fallback the rejection surfaces as whatever exception the first unimplemented
# piece of the machinery happens to raise — a `MethodError` from deep inside the time stepping
# loop, or, where a method's `initial_guess!` is not constrained to the problem types it
# implements, a `FieldError` from the cache — none of which say what is wrong. Methods that do
# use a nonlinear solver are caught earlier, by the same check in `initsolver`, so the wording
# here is that of `src/solvers.jl` verbatim: the rejection reads the same for every method,
# whether or not it solves.
#
# Every `integrate_step!` of an actual method is more specific than this, so nothing that is
# implemented is intercepted.
function integrate_step!(sol, history, params, int::GeometricIntegrator)
    throw(ArgumentError(string(
        nameof(typeof(method(int))), " does not support problems of type ",
        nameof(typeof(equation(problem(int)))), ".")))
end

"""
```julia
internal_variables(::Integrator) = NamedTuple()
```
Returns a `NamedTuple` containing all internal variables of an integrator that
shall be stored in an [`SolutionStep`](@ref). If there is no method for a
specific integrator implemented an empty `NamedTuple()` is returned.
"""
internal_variables(::GeometricIntegrator) = NamedTuple()
internal_variables(::Nothing) = NamedTuple()
