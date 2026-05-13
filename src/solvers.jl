
default_linesearch(method=nothing) = Backtracking()

default_options(method=nothing) = (
    min_iterations=1,
    f_abstol=8eps(),
    linesearch=default_linesearch(method),
    # verbosity=2,
)

initsolver(::SolverMethod, ::GeometricMethod, ::CacheDict; kwargs...) = NoSolver()

# create nonlinear solver
function initsolver(solvermethod::NonlinearSolverMethod, method::GeometricMethod, caches::CacheDict; kwargs...)
    x = zero(nlsolution(caches))
    y = zero(nlsolution(caches))
    NonlinearSolver(solvermethod, x, residual!, y; kwargs...)
end

# This accounts for the SimpleSolvers interface, expecting a single parameter argument,
# whereas the typical `residual!` methods expect a number of additional arguments.
residual!(y, x, parameters::Union{Tuple,NamedTuple}) = residual!(y, x, parameters...)
