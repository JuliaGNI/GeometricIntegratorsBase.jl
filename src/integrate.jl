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
    compute_vectorfields!(vectorfield(solstep, 0), solution(solstep, 0), problem(int))

    return solstep
end

"""
Solve for time steps n with n₁ ≤ n ≤ n₂.
```julia
integrate!(solution, integrator, n₁, n₂)
```
"""
function integrate!(sol::GeometricSolution, int::AbstractIntegrator, n₁::Int, n₂::Int)
    # check time steps range for consistency
    @assert n₁ ≥ 1
    @assert n₂ ≥ n₁
    @assert n₂ ≤ ntime(sol)

    # copy initial condition from solution to solutionstep and initialize
    solstep = solutionstep(int, sol[n₁-1])

    # loop over time steps
    for n in n₁:n₂
        # integrate one step and copy solution from cache to solution
        sol[n] = integrate!(solstep, int)

        # try
        #     sol[n] = integrate!(int)
        # catch ex
        #     tstr = " in time step " * string(n)
        #
        #     if m₁ ≠ m₂
        #         tstr *= " for initial condition " * string(m)
        #     end
        #
        #     tstr *= "."
        #
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
    end

    return sol
end

"""
Solve for all time steps n:
```julia
integrate!(solution, integrator)
```
"""
function integrate!(sol::GeometricSolution, int::AbstractIntegrator)
    integrate!(sol, int, 1, ntime(sol))
end


# Apply integrator for ntime time steps and return solution.
function integrate(integrator::AbstractIntegrator; kwargs...)
    solution = Solution(problem(integrator))
    integrate!(solution, integrator; kwargs...)
end
