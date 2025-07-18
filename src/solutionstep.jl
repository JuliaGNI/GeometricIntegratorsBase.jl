
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
Single solution step.

## Constructors:

```julia
SolutionStep(ics, parameters; kwargs...)
SolutionStep(problem; kwargs...)
```

Automatically construct the appropriate atomic solution based on the
given `equation` or `solution` type. If an `integrator` is provided as,
the `internal` field of the solution step is constructed according to
the internal state of the integrator as obtained from the function
`internal_variables`.

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

nhistory(::SolutionStep{ET,ST,VT,HT,IT,PT,NT}) where {ET,ST,VT,HT,IT,PT,NT} = NT

solution(solstep::SolutionStep) = solstep.solution
vectorfield(solstep::SolutionStep) = solstep.vectorfield
history(solstep::SolutionStep) = solstep.history
internal(solstep::SolutionStep) = solstep.internal
parameters(solstep::SolutionStep) = solstep.parameters

Base.keys(solstep::SolutionStep) = keys(solution(solstep))

solution(solstep::SolutionStep, i::Int) = NamedTuple{keys(solution(solstep))}(
    solution(solstep)[k][i] for k in keys(solution(solstep))
)

vectorfield(solstep::SolutionStep, i::Int) = NamedTuple{keys(vectorfield(solstep))}(
    vectorfield(solstep)[k][i] for k in keys(vectorfield(solstep))
)

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


function GeometricBase.reset!(solstep::SolutionStep, Δt)
    for k in keys(solstep)
        for i in backwardhistory(solstep)
            copy!(solution(solstep)[k][i], solution(solstep)[k][i-1])
        end
    end

    add!(solstep.t, Δt)

    return solstep
end


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

function update!(solstep::SolutionStep, Δ::NamedTuple)
    @assert keys(Δ) ⊆ keys(solstep)

    sol = solution(solstep)
    for k in keys(Δ)
        _update!(sol[k], Δ[k])
    end

    return solstep
end


enforce_periodicity!(solstep::SolutionStep, ::OffsetVector{<:AbstractVariable}) = nothing

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

function enforce_periodicity!(solstep::SolutionStep)
    # loop through all solution variables
    for s in solution(solstep)
        enforce_periodicity!(solstep, s)
    end
end
