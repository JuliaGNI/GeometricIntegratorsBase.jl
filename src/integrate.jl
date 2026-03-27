#*****************************************************************************#
# General integration functions for all integrators                           #
#*****************************************************************************#


"""
Parts of one integration step that are common to most if not all typical integrators
"""
function integrate!(solstep::SolutionStep, int::AbstractIntegrator)
    # reset solution step
    reset!(solstep, timestep(int))

    # compute initial guess
    initial_guess!(current(solstep), history(solstep), parameters(solstep), int)

    # integrate one initial condition for one time step
    integrate_step!(current(solstep), history(solstep), parameters(solstep), int)

    # copy internal variables from cache to solution step
    copy_internal_variables!(solstep, cache(int))

    # copy solver status to solution step
    # solver_status!(solver(int), solstep.internal[:solver])

    # take care of periodic solutions
    enforce_periodicity!(solstep)

    # update vector field for initial guess
    compute_vectorfields!(current(solstep), problem(int))

    return solstep
end

function integrate!(sol::GeometricSolution, int::AbstractIntegrator, n₁::Int, n₂::Int, solstep, curstate)
    # loop over time steps
    for n in n₁:n₂
        # integrate one step
        integrate!(solstep, int)

        # copy solution from solution step to solution
        copy!(sol, curstate, n)

        # check for NaNs
        if isnan(curstate)
            @warn "Solver encountered NaNs in solution at timestep n=$(n)."
            break
        end
    end
    return sol
end


"""
Solve for time steps n with n₁ ≤ n ≤ n₂.
```julia
integrate!(solution, integrator, n₁, n₂)
```
"""
function integrate!(sol::GeometricSolution, int::AbstractIntegrator, n₁::Int, n₂::Int; kwargs...)
    # check time steps range for consistency
    @assert n₁ ≥ 1
    @assert n₂ ≥ n₁
    @assert n₂ ≤ ntime(sol)

    # copy initial condition from solution to solutionstep and initialize
    solstep = solutionstep(int, sol[n₁-1]; kwargs...)
    curstate = current(solstep)

    integrate!(sol, int, n₁, n₂, solstep, curstate)

    # # loop over time steps
    # for n in n₁:n₂
    #     # integrate one step
    #     integrate!(solstep, int)

    #     # copy solution from solution step to solution
    #     copy!(sol, curstate, n)

    #     # check for NaNs
    #     if isnan(curstate)
    #         @warn "Solver encountered NaNs in solution at timestep n=$(n)."
    #         break
    #     end
    # end

    # try
    #     integrate!(solstep, int)
    #     copy!(sol, current(solstep), n)
    # catch ex
    #     tstr = " in time step " * string(n)

    #     if m₁ ≠ m₂
    #         tstr *= " for initial condition " * string(m)
    #     end

    #     tstr *= "."

    #     if isa(ex, DomainError)
    #         @warn("Domain error" * tstr)
    #     elseif isa(ex, ErrorException)
    #         @warn("Simulation exited early" * tstr)
    #         @warn(ex.msg)
    #     else
    #         @warn(string(typeof(ex)) * tstr)
    #         throw(ex)
    #     end
    # end

    return sol
end

"""
Solve for all time steps n:
```julia
integrate!(solution, integrator)
```
"""
function integrate!(sol::GeometricSolution, int::AbstractIntegrator; kwargs...)
    integrate!(sol, int, 1, ntime(sol); kwargs...)
end


# Apply integrator for ntime time steps and return solution.
function integrate(integrator::AbstractIntegrator; kwargs...)
    solution = Solution(problem(integrator))
    integrate!(solution, integrator; kwargs...)
end
