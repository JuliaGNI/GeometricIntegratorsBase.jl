
function Solution(problem::Union{EquationProblem,SubstepProblem}, args...; kwargs...)
    GeometricSolution(problem, args...; kwargs...)
end

function Solution(problem::EnsembleProblem, args...; kwargs...)
    EnsembleSolution(problem, args...; kwargs...)
end


function Base.setindex!(sol::GeometricSolution, solstep::SolutionStep, n)
    sol[n] = current(solstep)
end
