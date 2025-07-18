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
initmethod(method::GeometricMethod, ::AbstractProblem) = initmethod(method)

solversize(problem::AbstractProblemODE, method::GeometricMethod) = 0

internal_variables(::GeometricMethod, ::GeometricProblem) = NamedTuple()
nullvector(::GeometricMethod) = nothing
GeometricBase.tableau(::GeometricMethod) = missing

default_solver(::GeometricMethod) = NoSolver()
default_iguess(::GeometricMethod) = NoInitialGuess()
default_projection(::GeometricMethod) = NoProjection()

isodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
ispodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
ishodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
isiodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
islodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
issodemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false

isdaemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
ispdaemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
ishdaemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
isidaemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false
isldaemethod(::Union{GeometricMethod,Type{<:GeometricMethod}}) = false

isodemethod(::Union{ODEMethod,Type{<:ODEMethod}}) = true
ispodemethod(::Union{PODEMethod,Type{<:PODEMethod}}) = true
ishodemethod(::Union{HODEMethod,Type{<:HODEMethod}}) = true
isiodemethod(::Union{IODEMethod,Type{<:IODEMethod}}) = true
islodemethod(::Union{LODEMethod,Type{<:LODEMethod}}) = true
issodemethod(::Union{SODEMethod,Type{<:SODEMethod}}) = true

isdaemethod(::Union{DAEMethod,Type{<:DAEMethod}}) = true
ispdaemethod(::Union{PDAEMethod,Type{<:PDAEMethod}}) = true
ishdaemethod(::Union{HDAEMethod,Type{<:HDAEMethod}}) = true
isidaemethod(::Union{IDAEMethod,Type{<:IDAEMethod}}) = true
isldaemethod(::Union{LDAEMethod,Type{<:LDAEMethod}}) = true

isdelemethod(::Union{DELEMethod,Type{<:DELEMethod}}) = true

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
isenergypreserving(t::Type{<:GeometricMethod}) = applicable(t) ? isenergypreserving(t()) : missing
isstifflyaccurate(t::Type{<:GeometricMethod}) = applicable(t) ? isstifflyaccurate(t()) : missing

print_reference(io, method::GeometricMethod) =
    try
        ismissing(reference(method)) || print(io, reference(method))
    catch MethodError
        String("")
    end

# function check_symplecticity end
function symplecticity_conditions end
