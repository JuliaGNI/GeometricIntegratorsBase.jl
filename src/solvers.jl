
default_linesearch(method=nothing) = Backtracking()

default_options(method::GeometricMethod, problem::GeometricProblem) = (
    min_iterations=1,
    f_abstol=max(8, solversize(method, problem)) * eps(datatype(problem))
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
