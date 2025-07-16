
abstract type InitialGuess end


struct NoInitialGuess <: InitialGuess end

# function initialguess!(t, q, q̇, solstep::SolutionStepODE, ::AbstractProblemODE, ::NoInitialGuess)
#     t  = solstep.t̄
#     q .= solstep.q̄
#     q̇ .= solstep.v̄
#     (t = t, q = q, v = q̇)
# end

# function initialguess!(t, q, p, q̇, ṗ, solstep::SolutionStepPODE, ::Union{AbstractProblemPODE,AbstractProblemIODE}, ::NoInitialGuess)
#     t  = solstep.t̄
#     q .= solstep.q̄
#     p .= solstep.p̄
#     q̇ .= solstep.v̄
#     ṗ .= solstep.f̄
#     (t = t, q = q, p = p, v = q̇, f = ṗ)
# end
