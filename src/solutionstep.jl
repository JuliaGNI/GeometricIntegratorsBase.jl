
# q̄ = Symbol(join([Char('q'), Char(0x0304)]))
# q̇ = Symbol(join([Char('q'), Char(0x0307)]))
# q = Symbol(strip(String(:q̄), Char(0x0304)))
# q = Symbol(strip(String(:q̇), Char(0x0307)))

_strip_symbol(s::Symbol, c::Char) = Symbol(strip(normalize(String(s); decompose=true), c))
_strip_bar(s::Symbol) = _strip_symbol(s, Char(0x0304))
_strip_dot(s::Symbol) = _strip_symbol(s, Char(0x0307))

_add_symbol(s::Symbol, c::Char) = Symbol(normalize("$(s)$(c)"))
_add_bar(s::Symbol) = _add_symbol(s, Char(0x0304))
_add_dot(s::Symbol) = _add_symbol(s, Char(0x0307))

_state(x::Number) = zero(x)
_state(x::TimeVariable) = zero(x)
_state(x::StateVariable) = StateWithError(zero(x))
_state(x::VectorfieldVariable) = zero(x)
_state(x::AlgebraicVariable) = zero(x)
_state(x::StateWithError) = zero(x)

_vectorfield(x::Number) = missing
_vectorfield(x::TimeVariable) = missing
_vectorfield(x::StateVariable) = VectorfieldVariable(x)
_vectorfield(x::VectorfieldVariable) = missing
_vectorfield(x::AlgebraicVariable) = missing
_vectorfield(x::StateWithError) = _vectorfield(x.state)


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
    solutionType<:NamedTuple,
    vectorfieldType<:NamedTuple,
    historyType<:NamedTuple,
    internalType<:NamedTuple,
    paramsType<:OptionalParameters,
    nHistory
}

    solution::solutionType
    vectorfield::vectorfieldType
    history::historyType
    internal::internalType
    parameters::paramsType

    function SolutionStep{equType}(ics::NamedTuple, parameters::OptionalParameters; nhistory=1, internal=NamedTuple()) where {equType}
        @assert nhistory ≥ 1

        # create solution vector for all variables in ics
        solution = NamedTuple{keys(ics)}(
            OffsetVector([_state(x) for _ in 0:nhistory], 0:nhistory) for x in ics
        )

        # create vectorfield vector for all state variables in ics
        vectorfield = NamedTuple{keys(ics)}(
            OffsetVector([_vectorfield(x) for _ in 0:nhistory], 0:nhistory) for x in ics
        )

        vectorfield_filtered = NamedTuple{filter(k -> !all(ismissing.(vectorfield[k])), keys(vectorfield))}(vectorfield)

        vectorfield_keys = Tuple(_add_dot(k) for k in keys(vectorfield_filtered))
        vectorfield_dots = NamedTuple{vectorfield_keys}(values(vectorfield_filtered))

        history = merge(solution, vectorfield_dots)

        # create solstep
        solstep = new{equType,typeof(solution),typeof(vectorfield_filtered),typeof(history),typeof(internal),typeof(parameters),nhistory}(solution, vectorfield_filtered, history, internal, parameters)

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

function solutionstep(int::AbstractIntegrator, sol; kwargs...)
    solstep = SolutionStep(problem(int); internal=internal_variables(method(int), problem(int)), kwargs...)
    copy!(solstep, sol)
    solstep
end


@inline function Base.hasproperty(::SolutionStep{ET,ST,VT,HT}, s::Symbol) where {ET,ST,VT,HT}
    hasfield(ST, s) ||
        hasfield(ST, _strip_bar(s)) ||
        hasfield(VT, _strip_dot(s)) ||
        hasfield(SolutionStep, s)
end

@inline function Base.getproperty(sol::SolutionStep{ET,ST,VT,HT}, s::Symbol) where {ET,ST,VT,HT}
    s̄ = _strip_bar(s)
    ṡ = _strip_dot(s)
    if hasfield(ST, s)
        return getfield(sol, :solution)[s][0]
    elseif s̄ ≠ s && hasfield(ST, s̄)
        return getfield(sol, :solution)[s̄][1]
    elseif ṡ ≠ s && hasfield(VT, ṡ)
        return getfield(sol, :vectorfield)[ṡ][0]
    else
        return getfield(sol, s)
    end
end

@inline function Base.setproperty!(sol::SolutionStep{ET,ST,VT,HT}, s::Symbol, v) where {ET,ST,VT,HT}
    s̄ = _strip_bar(s)
    ṡ = _strip_dot(s)
    if hasfield(ST, s)
        return copy!(getfield(sol, :solution)[s][0], v)
    elseif s̄ ≠ s && hasfield(ST, s̄)
        return copy!(getfield(sol, :solution)[s̄][1], v)
    elseif ṡ ≠ s && hasfield(VT, ṡ)
        return copy!(getfield(sol, :vectorfield)[ṡ][0], v)
    else
        return setfield!(sol, s, v)
    end
end

# @inline Base.getindex(sol::SolutionStep, s::Symbol) = getproperty(sol, s)

# GeometricBase.datatype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = dType
# GeometricBase.timetype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = tType
# GeometricBase.arrtype(::SolutionStep{dType, tType, aType}) where {dType, tType, aType} = aType

"""
    nhistory(solstep::SolutionStep)

Return the number of previous time steps stored in the solution step.
"""
nhistory(::SolutionStep{ET,ST,VT,HT,IT,PT,NT}) where {ET,ST,VT,HT,IT,PT,NT} = NT

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
    history(solstep::SolutionStep)

Return the history field of the solution step, which provides convenient access
to both solution and vectorfield data at the current and previous time steps.
"""
history(solstep::SolutionStep) = solstep.history

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
    keys(solstep::SolutionStep)

Return the keys of the solution variables in the solution step.
"""
Base.keys(solstep::SolutionStep) = keys(solution(solstep))

"""
    solution(solstep::SolutionStep, i::Int)

Return a NamedTuple with the solution at time step `i`, where `i=0` is the current
time step, `i=1` is the previous time step, etc.
"""
solution(solstep::SolutionStep, i::Int) = NamedTuple{keys(solution(solstep))}(
    solution(solstep)[k][i] for k in keys(solution(solstep))
)

"""
    vectorfield(solstep::SolutionStep, i::Int)

Return a NamedTuple with the vectorfield at time step `i`, where `i=0` is the current
time step, `i=1` is the previous time step, etc.
"""
vectorfield(solstep::SolutionStep, i::Int) = NamedTuple{keys(vectorfield(solstep))}(
    vectorfield(solstep)[k][i] for k in keys(vectorfield(solstep))
)

"""
    history(solstep::SolutionStep, i::Int)

Return a NamedTuple with the history (solution and vectorfield) at time step `i`,
where `i=0` is the current time step, `i=1` is the previous time step, etc.
"""
history(solstep::SolutionStep, i::Int) = NamedTuple{keys(history(solstep))}(
    history(solstep)[k][i] for k in keys(history(solstep))
)

"""
Returns a NamedTuple with the solution of the current time step.
"""
current(solstep::SolutionStep) = history(solstep, 0)

"""
Returns a NamedTuple with the solution of the previous time step.
"""
previous(solstep::SolutionStep) = history(solstep, 1)

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
    for k in keys(solstep)
        for i in backwardhistory(solstep)
            copy!(solution(solstep)[k][i], solution(solstep)[k][i-1])
        end
    end

    add!(solstep.t, Δt)

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
    @assert keys(sol) ⊆ keys(solstep)

    for k in keys(sol)
        copy!(solution(solstep)[k][0], sol[k])
    end

    # for i in backwardhistory(solstep)
    #     solution(solstep)[k][i] = 0
    # end

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
# _update!(x::OffsetVector{<:AbstractScalarVariable}, Δx) = x[0] += Δx

_update!(x::OffsetVector{<:AbstractStateVariable}, Δx) = add!(x[0], Δx)
_update!(x::OffsetVector{<:AlgebraicVariable}, Δx) = copy!(x[0], Δx)

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

    sol = solution(solstep)
    for k in keys(Δ)
        _update!(sol[k], Δ[k])
    end

    return solstep
end


"""
    enforce_periodicity!(solstep::SolutionStep, ::OffsetVector{<:AbstractVariable})

No-op method for non-periodic variables.

This method is called for variable types that do not support periodicity
constraints and performs no operation.

# Arguments
- `solstep::SolutionStep`: The solution step (unused)
- `::OffsetVector{<:AbstractVariable}`: The variable vector (unused)
"""
enforce_periodicity!(solstep::SolutionStep, ::OffsetVector{<:AbstractVariable}) = nothing

"""
    enforce_periodicity!(solstep::SolutionStep, s::OffsetVector{<:Union{<:StateVariable,<:StateVariableWithError}})

Enforce periodic boundary conditions on a state variable.

This function checks each component of the state variable for periodic
boundary conditions and adjusts values that fall outside the specified
range by adding or subtracting the range size. The adjustment is applied
to both the current solution and all historical values to maintain
consistency for initial guesses in iterative solvers.

# Arguments
- `solstep::SolutionStep`: The solution step containing the variable
- `s::OffsetVector{<:Union{<:StateVariable,<:StateVariableWithError}}`:
  The state variable to enforce periodicity on

# Details
For each component `i` of the state variable:
- If `isperiodic(s[0], i)` is true, check if `s[0][i]` is within `range(s[0], i)`
- If below the range, add the range size until within bounds
- If above the range, subtract the range size until within bounds
- Apply the same adjustment to all historical values `s[j][i]` for consistency
"""
function enforce_periodicity!(solstep::SolutionStep, s::OffsetVector{<:Union{<:StateVariable,<:StateVariableWithError}})
    # loop through all components of state variable s
    for i in eachindex(s[0])
        # check if component i is periodic and if so outside of its range
        if isperiodic(s[0], i)# && !verifyrange(s[0], i)
            # compute range size Δs
            Δs = (range(s[0], i)[end] - range(s[0], i)[begin])
            # add Δs as long as value is below range
            while s[0][i] < range(s[0], i)[begin]
                s[0][i] += Δs
                # add Δs for all entries of history to allow for consistent initial guesses
                for j in eachsolution(solstep)
                    s[j][i] += Δs
                end
            end
            # subtract Δs as long as value is above range
            while s[0][i] > range(s[0], i)[end]
                s[0][i] -= Δs
                # subtract Δs for all entries of history to allow for consistent initial guesses
                for j in eachsolution(solstep)
                    s[j][i] -= Δs
                end
            end
        end
    end
end

"""
    enforce_periodicity!(solstep::SolutionStep)

Enforce periodic boundary conditions on all solution variables.

This function iterates through all variables in the solution step and
calls `enforce_periodicity!` on each one. Variables that support
periodicity will have their boundary conditions enforced, while others
will be unaffected.

# Arguments
- `solstep::SolutionStep`: The solution step to enforce periodicity on
"""
function enforce_periodicity!(solstep::SolutionStep)
    # loop through all solution variables
    for s in solution(solstep)
        enforce_periodicity!(solstep, s)
    end
end
