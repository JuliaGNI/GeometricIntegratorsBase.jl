
const default_extrapolation_stages = 5

# abstract type Extrapolation <: DeterministicMethod end
abstract type Extrapolation end

# default_extrapolation() = missing
default_extrapolation() = MidpointExtrapolation(default_extrapolation_stages)

struct NoExtrapolation <: Extrapolation end


function extrapolate!(newsol, oldsol, ::GeometricProblem, extrap::Union{NoExtrapolation,NoInitialGuess})
    for k in keys(newsol)
        if k != :t
            newsol[k] .= oldsol[k]
        end
    end
    return newsol
end

solutionstep!(sol, history, ::GeometricProblem, extrap::Union{NoExtrapolation,NoInitialGuess}) = sol


# """
# """
# function extrapolate! end



# function extrapolate!(t₀, q₀, q̇₀, t₁, q₁, q̇₁, t₂, q₂, q̇₂, problem::AbstractProblemODE, extrap::Extrapolation; kwargs...)
#     solution = (
#         t = t₀, q = q₀, v = q̇₀,
#     )

#     history = (
#         (t = t₁, q = q₁, v = q̇₁),
#         (t = t₂, q = q₂, v = q̇₂),
#     )

#     solutionstep!(solution, history, problem, extrap; kwargs...)
# end
