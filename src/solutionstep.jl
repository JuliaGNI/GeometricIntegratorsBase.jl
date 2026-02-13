"""
Holds the solution of a geometric equation at a single time step.

It stores all the information that is passed from one time step to the next.
This includes the current solution, the vectorfield of the equation, and solution
data from previous time steps.

## Type Parameters

* `equationType`: type of the geometric equation
* `solutionType`: type of the solution tuple
* `vectorfieldType`: type of the vectorfield tuple
* `historyType`: type of the history tuple
* `internalType`: type of the internal variables tuple
* `paramsType`: type of the parameters
* `nHistory`: number of previous time steps to store

## Fields

* `solution`: a `NamedTuple` of `OffsetVector`s holding the solution of the current
  and previous `nHistory` time steps. The indices of the `OffsetVector` are `0...nHistory`.
  `solution[k][0]` is the current solution for variable `k`, `solution[k][1]` the solution
  at the previous time step, and so on.
* `vectorfield`: a `NamedTuple` of `OffsetVector`s holding the vector field of the
  current and previous `nHistory` time steps.
* `history`: a `NamedTuple` that provides convenient access to solution and vectorfield
  of the current and previous time steps.
* `internal`: a `NamedTuple` for integrator-specific internal variables.
* `parameters`: the parameters of the equation.

## Constructors

```julia
SolutionStep{equType}(ics::NamedTuple, parameters::OptionalParameters; nhistory=1, internal=NamedTuple())
SolutionStep(problem::GeometricProblem; kwargs...)
```

The constructor `SolutionStep{equType}(...)` automatically constructs the appropriate
solution step object from the given initial conditions `ics` and `parameters`.
The `internal` field of the solution step is for integrator-specific internal state.
The `solutionstep(integrator, ...)` function is a convenient wrapper to construct
a `SolutionStep` with the correct internal variables for a given integrator.
"""
struct SolutionStep{
    equationType<:GeometricEquation,
    nHistory,
    stateType<:OffsetArray,
    solutionType<:OffsetArray,
    vectorfieldType<:OffsetArray,
    internalType<:NamedTuple,
    paramsType<:OptionalParameters,
    currentType<:State,
    previousType<:State
}

    state::stateType
    solution::solutionType
    vectorfield::vectorfieldType

    internal::internalType
    parameters::paramsType

    current::currentType
    previous::previousType

    function SolutionStep{equType}(ics::NamedTuple, parameters::OptionalParameters=NullParameters(); nhistory=2, internal=NamedTuple()) where {equType}
        @assert nhistory ≥ 1

        # create state vector according to the variables in ics
        states = OffsetVector([State(ics; initialize=false) for _ in 0:nhistory], 0:nhistory)

        # create solution vector from the solution variables in states
        solutions = OffsetVector([solution(st) for st in states], 0:nhistory)

        # create vectorfield vector from the vectorfield variables in states
        vectorfields = OffsetVector([vectorfield(st) for st in states], 0:nhistory)

        # create wrappers for current and previous state
        current = states[0]
        previous = HistoryState(states[1])

        # create solstep
        solstep = new{equType,nhistory,
            typeof(states),
            typeof(solutions),
            typeof(vectorfields),
            typeof(internal),
            typeof(parameters),
            typeof(current),
            typeof(previous)}(
            states, solutions, vectorfields, internal, parameters, current, previous
        )

        # copy initial conditions to current solution of solstep
        copy!(solstep, ics)

        return solstep
    end
end

function SolutionStep(problem::GeometricProblem{superType}; kwargs...) where {superType<:GeometricEquation}
    SolutionStep{superType}(initial_conditions(problem), parameters(problem); kwargs...)
end

# function SolutionStep(problem::GeometricProblem, extrap::Extrapolation = default_extrapolation(); kwargs...)
#     solstep = SolutionStep(initial_conditions(problem), parameters(problem); kwargs...)
#     initialize!(solstep, problem, extrap)
#     return solstep
# end

function solutionstep(int::AbstractIntegrator, sol; extrap::Extrapolation=default_extrapolation(), kwargs...)
    # create solutionstep
    solstep = SolutionStep(problem(int); internal=internal_variables(method(int), problem(int)), kwargs...)

    # copy initial conditions from sol
    copy!(solstep, sol)

    # compute vector fields for initial conditions
    compute_vectorfields!(vectorfield(solstep, 0), solution(solstep, 0), problem(int))

    # initialize solution step history
    initialize!(solstep, problem(int), extrap)

    return solstep
end

@inline hasstatevariable(::SolutionStep{ET,NH,ST,SOLT,VT,IT,PT,CT,HT}, s::Symbol) where {CST,HST,ET,NH,ST,SOLT,VT,IT,PT,CT<:State{CST},HT<:State{HST}} = hasfield(CST, s)
@inline hashistoryvariable(::SolutionStep{ET,NH,ST,SOLT,VT,IT,PT,CT,HT}, s::Symbol) where {CST,HST,ET,NH,ST,SOLT,VT,IT,PT,CT<:State{CST},HT<:State{HST}} = hasfield(HST, s)

@inline function Base.hasproperty(sol::SolutionStep, s::Symbol)
    hasfield(SolutionStep, s) || hasstatevariable(sol, s) || hashistoryvariable(sol, s)
end

@inline function Base.getproperty(sol::SolutionStep, s::Symbol)
    if hasstatevariable(sol, s)
        return state(sol.current)[s]
    elseif hashistoryvariable(sol, s)
        return state(sol.previous)[s]
    else
        return getfield(sol, s)
    end
end

@inline function Base.setproperty!(sol::SolutionStep, s::Symbol, x)
    if hasstatevariable(sol, s)
        return copy!(state(sol.current)[s], x)
    elseif hashistoryvariable(sol, s)
        return copy!(state(sol.previous)[s], x)
    else
        return setfield!(sol, s, x)
    end
end


"""
    keys(solstep::SolutionStep)

Return the keys of the state variables in the solution step.
"""
Base.keys(solstep::SolutionStep) = keys(current(solstep))

"""
    nhistory(solstep::SolutionStep)

Return the number of previous time steps stored in the solution step.
"""
nhistory(::SolutionStep{ET,NT}) where {ET,NT} = NT

"""
    state(solstep::SolutionStep)

Return the state field of the solution step, which contains the state
vectors for all variables at the current and previous time steps.
"""
state(solstep::SolutionStep) = solstep.state

history(solstep::SolutionStep) = state(solstep)

"""
    solution(solstep::SolutionStep)

Return the solution field of the solution step, which contains the solution
vectors for all variables at the current and previous time steps.
"""
solution(solstep::SolutionStep) = solstep.solution

"""
    vectorfield(solstep::SolutionStep)

Return the vectorfield field of the solution step, which contains the
vectorfield evaluations for all variables at the current and previous time steps.
"""
vectorfield(solstep::SolutionStep) = solstep.vectorfield

"""
    internal(solstep::SolutionStep)

Return the internal field of the solution step, which contains integrator-specific
internal variables.
"""
internal(solstep::SolutionStep) = solstep.internal

"""
    parameters(solstep::SolutionStep)

Return the parameters field of the solution step, which contains the parameters
of the geometric equation.
"""
parameters(solstep::SolutionStep) = solstep.parameters

"""
    state(solstep::SolutionStep, i::Int)

Return a [`State`](@ref) with the state at time step `i`, where `i=0` is the current
time step, `i=1` is the previous time step, etc.
"""
state(solstep::SolutionStep, i::Int) = state(solstep)[i]

"""
    solution(solstep::SolutionStep, i::Int)

Return a `NamedTuple` with the solution at time step `i`, where `i=0` is the current
time step, `i=1` is the previous time step, etc.
"""
solution(solstep::SolutionStep, i::Int) = solution(solstep)[i]

"""
    vectorfield(solstep::SolutionStep, i::Int)

Return a NamedTuple with the vectorfield at time step `i`, where `i=0` is the current
time step, `i=1` is the previous time step, etc.
"""
vectorfield(solstep::SolutionStep, i::Int) = vectorfield(solstep)[i]

"""
Returns a NamedTuple with the solution of the current time step.
"""
current(solstep::SolutionStep) = state(solstep, 0)

"""
Returns a NamedTuple with the solution of the previous time step.
"""
previous(solstep::SolutionStep) = state(solstep, 1)

eachsolution(sol::SolutionStep) = 1:nhistory(sol)
backwardhistory(sol::SolutionStep) = nhistory(sol):-1:1


"""
    reset!(solstep::SolutionStep, Δt)

Reset the solution step for the next time step by shifting the solution history
backward and advancing the time by `Δt`.

This function moves the current solution to the previous time step position,
the previous solution to the one before that, and so on. The time variable
is incremented by `Δt`. This prepares the solution step for computing the
next time step.

# Arguments
- `solstep`: the solution step to reset
- `Δt`: the time step size to advance
"""
function GeometricBase.reset!(solstep::SolutionStep, Δt)
    for i in backwardhistory(solstep)
        copy!(state(solstep, i), state(solstep, i - 1))
    end

    # add!(solstep.t, Δt)
    current(solstep).t += Δt

    return solstep
end


"""
    copy!(solstep::SolutionStep, sol::NamedTuple)

Copy the values from a `NamedTuple` `sol` to the current time step of the solution step.

The keys of `sol` must be a subset of the keys of the solution step. Only the
current time step (index 0) of the solution step is modified.

# Arguments
- `solstep`: the solution step to copy into
- `sol`: the named tuple containing the solution values to copy
"""
function Base.copy!(solstep::SolutionStep, sol::NamedTuple)
    copy!(current(solstep), sol)
    return solstep
end

"""
Copy the initial conditions of a `EquationProblem` to the current state of a solution step.
"""
function Base.copy!(solstep::SolutionStep, equ::EquationProblem)
    copy!(solstep, initial_conditions(equ))
end


# The following is a workaround for scalar solution entries that are
# stored as plain Julia numbers instead of AbstractScalarVariable.
# Typically, this concerns the time which should be a TimeVariable.
# _update!(x::OffsetVector{DT, <:AbstractArray{DT}}, y::OffsetVector{DT, <:AbstractArray{DT}}, i, j) where {DT} = x[i] = y[j]
# _update!(x::OffsetVector{<:AbstractScalarVariable}, Δx) = x += Δx

_update!(x::AbstractStateVariable, Δx) = add!(x, Δx)
_update!(x::AlgebraicVariable, y) = copy!(x, y)

"""
    update!(solstep::SolutionStep, Δ::NamedTuple)

Update the current solution in the solution step by adding increments.

This function applies increments to the current solution variables (at index 0)
in the solution step. The increments are added using the appropriate `add!`
method for each variable type, which handles different variable types correctly
(e.g., compensated summation for StateWithError variables).

# Arguments
- `solstep::SolutionStep`: The solution step to update
- `Δ::NamedTuple`: Named tuple containing increments for each variable to update.
  Keys must be a subset of the solution step's variable keys.

# Returns
- `solstep`: The updated solution step (for method chaining)

# Throws
- `AssertionError`: If any key in `Δ` is not present in the solution step

# Examples
```julia
# Create a solution step
solstep = SolutionStep{MyEquation}(initial_conditions, parameters)

# Update with increments
update!(solstep, (q = [0.1, 0.2], p = [0.05, 0.1]))
```
"""
function update!(solstep::SolutionStep, Δ::NamedTuple)
    @assert keys(Δ) ⊆ keys(solstep)

    sol = current(solstep)
    for k in keys(Δ)
        _update!(sol[k], Δ[k])
    end

    return solstep
end


_var(solstep::SolutionStep, s::Symbol, i=0) = state(solstep)[i][s]

"""
    enforce_periodicity!(solstep::SolutionStep)

Enforce periodic boundary conditions on a state variable.

This function checks each component of the state variable for periodic
boundary conditions and adjusts values that fall outside the specified
range by adding or subtracting the range size. The adjustment is applied
to both the current solution and all historical values to maintain
consistency for initial guesses in iterative solvers.

# Arguments
- `solstep::SolutionStep`: The solution step to enforce periodicity on

# Details
This function iterates through all variables in the solution step and
calls `enforce_periodicity!` on each one. Variables that support
periodicity will have their boundary conditions enforced, while others
will be unaffected.

For each component `i` of the state variable:
- If `isperiodic(s[0], i)` is true, check if `s[0][i]` is within `range(s[0], i)`
- If below the range, add the range size until within bounds
- If above the range, subtract the range size until within bounds
- Apply the same adjustment to all historical values `s[j][i]` for consistency
"""
function enforce_periodicity!(solstep::SolutionStep)
    # loop through all solution variables
    for s in keys(current(solstep))
        v = _var(solstep, s, 0)
        # loop through all components of state variable v
        for i in eachindex(v)
            # check if component i is periodic and if so outside of its range
            if isperiodic(v, i)# && !verifyrange(s[0], i)
                # compute range size Δs
                Δs = (range(v, i)[end] - range(v, i)[begin])
                # add Δs as long as value is below range
                while v[i] < range(v, i)[begin]
                    v[i] += Δs
                    # add Δs for all entries of history to allow for consistent initial guesses
                    for j in eachsolution(solstep)
                        _var(solstep, s, j)[i] += Δs
                    end
                end
                # subtract Δs as long as value is above range
                while v[i] > range(v, i)[end]
                    v[i] -= Δs
                    # subtract Δs for all entries of history to allow for consistent initial guesses
                    for j in eachsolution(solstep)
                        _var(solstep, s, j)[i] -= Δs
                    end
                end
            end
        end
    end
end


function initialize!(solstep::SolutionStep, problem::GeometricProblem, extrap::Extrapolation=default_extrapolation())
    for i in eachsolution(solstep)
        state(solstep, i)[:t] .= state(solstep, i - 1).t .- timestep(problem)
        extrapolate!(state(solstep, i), state(solstep, i - 1), problem, extrap)
        compute_vectorfields!(vectorfield(solstep, i), solution(solstep, i), problem)
    end
    return solstep
end

function initialize!(solstep::SolutionStep, problem::DELEProblem, extrap::Extrapolation=default_extrapolation())
    state(solstep, 1)[:t] .= state(solstep, 0).t .- timestep(problem)
    state(solstep, 1)[:q] .= initial_conditions(problem).q̄

    return solstep
end
