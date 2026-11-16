"""
`GeometricMethod` is the abstract supertype for all integration methods implemented in GeometricIntegrators.
"""
abstract type GeometricMethod <: AbstractMethod end
abstract type DeterministicMethod <: GeometricMethod end

abstract type ODEMethod <: DeterministicMethod end
abstract type PODEMethod <: DeterministicMethod end
abstract type HODEMethod <: DeterministicMethod end
abstract type IODEMethod <: DeterministicMethod end
abstract type LODEMethod <: DeterministicMethod end
abstract type SODEMethod <: DeterministicMethod end

abstract type DAEMethod <: DeterministicMethod end
abstract type PDAEMethod <: DeterministicMethod end
abstract type HDAEMethod <: DeterministicMethod end
abstract type IDAEMethod <: DeterministicMethod end
abstract type LDAEMethod <: DeterministicMethod end

abstract type DELEMethod <: DeterministicMethod end

initmethod(method::GeometricMethod) = method
initmethod(method::GeometricMethod, ::GeometricProblem) = initmethod(method)
solversize(method::GeometricMethod, ::GeometricProblem) = 0

internal_variables(::GeometricMethod, ::GeometricProblem) = NamedTuple()
nullvector(::GeometricMethod) = nothing
GeometricBase.tableau(::GeometricMethod) = missing

default_solver(::GeometricMethod) = NoSolver()
default_iguess(::GeometricMethod) = NoInitialGuess()
default_projection(::GeometricMethod) = NoProjection()

@doc raw"""
    default_extrapolation(method)

The extrapolation used to fill a `SolutionStep`'s history before the first step.

A `SolutionStep` carries states at ``t_{-1}, t_{-2}, \ldots`` that no initial condition provides,
and `initialize!` manufactures them by running the problem *backwards* from ``t_0``. For a
deterministic problem this is well defined, and `MidpointExtrapolation` does it accurately, which
is why it is the default.

It is dispatched on the method rather than fixed so that a family of methods for which running the
problem backwards is not meaningful can say so. A stochastic problem is the case in point — its
past is not a function of its present, so there is nothing to extrapolate to, and such a method
returns `NoExtrapolation()`, which copies the initial state into the history slots instead.
"""
default_extrapolation(::GeometricMethod) = MidpointExtrapolation(default_extrapolation_stages)

isodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ispodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ishodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isiodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
islodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
issodemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

isdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ispdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ishdaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isidaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
isldaemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false

isodemethod(::Union{ODEMethod, Type{<:ODEMethod}}) = true
ispodemethod(::Union{PODEMethod, Type{<:PODEMethod}}) = true
ishodemethod(::Union{HODEMethod, Type{<:HODEMethod}}) = true
isiodemethod(::Union{IODEMethod, Type{<:IODEMethod}}) = true
islodemethod(::Union{LODEMethod, Type{<:LODEMethod}}) = true
issodemethod(::Union{SODEMethod, Type{<:SODEMethod}}) = true

isdaemethod(::Union{DAEMethod, Type{<:DAEMethod}}) = true

isdelemethod(::Union{GeometricMethod, Type{<:GeometricMethod}}) = false
ispdaemethod(::Union{PDAEMethod, Type{<:PDAEMethod}}) = true
ishdaemethod(::Union{HDAEMethod, Type{<:HDAEMethod}}) = true
isidaemethod(::Union{IDAEMethod, Type{<:IDAEMethod}}) = true
isldaemethod(::Union{LDAEMethod, Type{<:LDAEMethod}}) = true

isdelemethod(::Union{DELEMethod, Type{<:DELEMethod}}) = true

const METHOD_DEFAULTS = (initmethod, nullvector, tableau, default_solver, default_iguess,
    default_projection, default_extrapolation)

const METHOD_FAMILIES = (ODEMethod => isodemethod, PODEMethod => ispodemethod,
    HODEMethod => ishodemethod, IODEMethod => isiodemethod, LODEMethod => islodemethod,
    SODEMethod => issodemethod, DAEMethod => isdaemethod, PDAEMethod => ispdaemethod,
    HDAEMethod => ishdaemethod, IDAEMethod => isidaemethod, LDAEMethod => isldaemethod,
    DELEMethod => isdelemethod)

"""
    interface_violations(method::GeometricMethod; include_missing = false)

The violations of the contract of `GeometricMethod` by `method`, in the form of
`GeometricBase.interface_violations`. The entries of the method for `AbstractMethod`, which checks
the ten properties of GeometricBase and receives `include_missing`, come first. The entries of the
defaults of `GeometricMethod` that take no problem follow, in the order `initmethod`,
`nullvector`, `tableau`, `default_solver`, `default_iguess`, `default_projection` and
`default_extrapolation`, each called on `method`, and then the entry of the family trait:

- `:throws`: the call throws; `detail` names the exception's type. An `InterruptException` is
  thrown again.
- `:type`: `initmethod` answers other than a `GeometricMethod`, `default_solver` other than a
  `GeometricBase.SolverMethod`, `default_iguess` neither an `InitialGuess` nor an `Extrapolation`,
  `default_projection` other than a `ProjectionMethod`, or `default_extrapolation` other than an
  `Extrapolation`; `detail` names the type of the answer. `nullvector` and `tableau` have no type
  rule.
- `:value`: `method` is a subtype of one of the twelve family types, such as `ODEMethod`, and the
  trait of that family, such as `isodemethod`, answers other than `true`; `detail` names the
  family. The other eleven traits have no rule, because a method may be implemented for more kinds
  of problem than its family, and a method under no family gets no entry of this kind.
"""
function GeometricBase.interface_violations(method::GeometricMethod;
        include_missing::Bool = false)
    violations = invoke(GeometricBase.interface_violations, Tuple{AbstractMethod}, method;
        include_missing)
    for f in METHOD_DEFAULTS
        result = call_default!(violations, f, method)
        result === nothing && continue
        answer = something(result)
        default_has_type(f, answer) ||
            push!(violations,
                (; f = nameof(f), kind = :type,
                    detail = "answers a value of type $(typeof(answer))"))
    end
    for (family, trait) in METHOD_FAMILIES
        method isa family || continue
        result = call_default!(violations, trait, method)
        result === nothing && continue
        answer = something(result)
        answer === true ||
            push!(violations,
                (; f = nameof(trait), kind = :value,
                    detail = "answers `$(repr(answer))` on a method of the family `$family`"))
    end
    return violations
end

# Calls `f(method)` and answers `Some(answer)`; where the call throws, pushes a `:throws` entry
# and answers `nothing`. An `InterruptException` is thrown again.
function call_default!(violations, f, method)
    try
        return Some(f(method))
    catch e
        e isa InterruptException && rethrow()
        push!(violations, (; f = nameof(f), kind = :throws, detail = "throws $(typeof(e))"))
        return nothing
    end
end

function default_has_type(f, answer)
    f === initmethod ? answer isa GeometricMethod :
    f === default_solver ? answer isa GeometricBase.SolverMethod :
    f === default_iguess ? answer isa Union{InitialGuess, Extrapolation} :
    f === default_projection ? answer isa ProjectionMethod :
    f === default_extrapolation ? answer isa Extrapolation : true
end

isexplicit(::GeometricMethod) = missing
isimplicit(::GeometricMethod) = missing
issymmetric(::GeometricMethod) = missing
issymplectic(::GeometricMethod) = missing
isenergypreserving(::GeometricMethod) = missing
isstifflyaccurate(::GeometricMethod) = missing

isexplicit(t::Type{<:GeometricMethod}) = applicable(t) ? isexplicit(t()) : missing
isimplicit(t::Type{<:GeometricMethod}) = applicable(t) ? isimplicit(t()) : missing
issymmetric(t::Type{<:GeometricMethod}) = applicable(t) ? issymmetric(t()) : missing
issymplectic(t::Type{<:GeometricMethod}) = applicable(t) ? issymplectic(t()) : missing
function isenergypreserving(t::Type{<:GeometricMethod})
    applicable(t) ? isenergypreserving(t()) : missing
end
function isstifflyaccurate(t::Type{<:GeometricMethod})
    applicable(t) ? isstifflyaccurate(t()) : missing
end

reference(::GeometricMethod) = missing

function print_reference(io, method::GeometricMethod)
    ismissing(reference(method)) || print(io, reference(method))
end

# function check_symplecticity end
function symplecticity_conditions end
