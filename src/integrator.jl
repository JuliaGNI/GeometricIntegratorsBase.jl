
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
    MT<:GeometricMethod,
    PT<:AbstractProblem,
    CT<:CacheDict{PT,MT},
    ST<:Union{NonlinearSolver,SolverMethod},
    IT<:Union{InitialGuess,Extrapolation}
} <: AbstractIntegrator

    problem::PT
    method::MT
    caches::CT
    solver::ST
    iguess::IT
end

function GeometricIntegrator(
    problem::AbstractProblem,
    integratormethod::GeometricMethod,
    solvermethod::SolverMethod,
    iguess::Union{InitialGuess,Extrapolation};
    method=initmethod(integratormethod, problem),
    caches=CacheDict(problem, method),
    options...
)
    solver = initsolver(solvermethod, method, caches; (length(options) == 0 ? default_options(integratormethod) : options)...)
    GeometricIntegrator(problem, method, caches, solver, iguess)
end

function GeometricIntegrator(
    problem::AbstractProblem,
    method::GeometricMethod;
    solver=default_solver(method),
    initialguess=default_iguess(method),
    kwargs...
)
    GeometricIntegrator(problem, method, solver, initialguess; kwargs...)
end

GeometricIntegrator(::AbstractProblem, ::Nothing, args...; kwargs...) = nothing

problem(int::GeometricIntegrator) = int.problem
method(int::GeometricIntegrator) = int.method
caches(int::GeometricIntegrator) = int.caches
solver(int::GeometricIntegrator) = int.solver
iguess(int::GeometricIntegrator) = int.iguess
initialguess(int::GeometricIntegrator) = int.iguess

cache(int::GeometricIntegrator, DT) = caches(int)[DT]
cache(int::GeometricIntegrator) = cache(int, datatype(problem(int)))
hasnullvector(int::GeometricIntegrator) = hasnullvector(method(int))
implicit_update(int::GeometricIntegrator) = implicit_update(method(int))
nconstraints(int::GeometricIntegrator) = nconstraints(problem(int))
Base.ndims(int::GeometricIntegrator) = ndims(problem(int))
nstages(int::GeometricIntegrator) = nstages(tableau(method(int)))
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

"""
function integrate_step! end


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
