
function Solution(problem::Union{EquationProblem,SubstepProblem}, args...; kwargs...)
    GeometricSolution(problem, args...; kwargs...)
end

function Solution(problem::EnsembleProblem, args...; kwargs...)
    EnsembleSolution(problem, args...; kwargs...)
end
