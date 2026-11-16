using GeometricIntegratorsBase
using Test

import GeometricBase
include("helpers/harmonic_oscillator.jl")
import .HarmonicOscillator: odeproblem

prob = odeproblem()

struct TestMethod <: GeometricMethod end

@test initmethod(TestMethod()) == TestMethod()
@test internal_variables(TestMethod(), prob) == NamedTuple()
@test nullvector(TestMethod()) === nothing
@test ismissing(tableau(TestMethod()))

@test default_solver(TestMethod()) == NoSolver()
@test default_iguess(TestMethod()) == NoInitialGuess()
@test default_projection(TestMethod()) == NoProjection()
@test default_extrapolation(TestMethod()) == MidpointExtrapolation()

@test isodemethod(TestMethod) == isodemethod(TestMethod()) == false
@test ispodemethod(TestMethod) == ispodemethod(TestMethod()) == false
@test ishodemethod(TestMethod) == ishodemethod(TestMethod()) == false
@test isiodemethod(TestMethod) == isiodemethod(TestMethod()) == false
@test islodemethod(TestMethod) == islodemethod(TestMethod()) == false
@test issodemethod(TestMethod) == issodemethod(TestMethod()) == false

@test isdaemethod(TestMethod) == isdaemethod(TestMethod()) == false
@test ispdaemethod(TestMethod) == ispdaemethod(TestMethod()) == false
@test ishdaemethod(TestMethod) == ishdaemethod(TestMethod()) == false
@test isidaemethod(TestMethod) == isidaemethod(TestMethod()) == false
@test isldaemethod(TestMethod) == isldaemethod(TestMethod()) == false

@test ismissing(isexplicit(TestMethod()))
@test ismissing(isimplicit(TestMethod()))
@test ismissing(issymmetric(TestMethod()))
@test ismissing(issymplectic(TestMethod()))
@test ismissing(isenergypreserving(TestMethod()))
@test ismissing(isstifflyaccurate(TestMethod()))

# The method properties are asserted through `GeometricBase` on purpose, not through the bare
# names. The bare names are what this package exports, so a test written with them passes whether
# or not they are the same functions the rest of the ecosystem uses. `RungeKutta` and
# `GeometricIntegrators` attach their properties to `GeometricBase`'s generics, so those are the
# ones that have to answer.
@test issymplectic === GeometricBase.issymplectic
@test GeometricBase.issymplectic(ImplicitMidpoint())
@test GeometricBase.isexplicit(ExplicitEuler())

# `print_reference` reads the reference off `GeometricBase.reference`, where every method in the
# ecosystem defines it.
GeometricBase.reference(::TestMethod) = "n/a"
@test sprint(GeometricIntegratorsBase.print_reference, TestMethod()) == "n/a"

# The fixtures of `interface_violations`. Each fixture has its own type, because a type defined
# again with a changed definition is a new type. `fk` drops the `detail` of each entry.
fk(violations) = [(v.f, v.kind) for v in violations]

struct ThrowingPropertyMethod <: GeometricMethod end
GeometricBase.issymplectic(::ThrowingPropertyMethod) = error("issymplectic")

struct ThrowingPropertyWrongSolverMethod <: GeometricMethod end
GeometricBase.issymplectic(::ThrowingPropertyWrongSolverMethod) = error("issymplectic")
GeometricIntegratorsBase.default_solver(::ThrowingPropertyWrongSolverMethod) = 1

struct NoPropertyMethod <: GeometricMethod end

struct WrongSolverMethod <: GeometricMethod end
GeometricIntegratorsBase.default_solver(::WrongSolverMethod) = 1
struct WrongIguessMethod <: GeometricMethod end
GeometricIntegratorsBase.default_iguess(::WrongIguessMethod) = :none
struct WrongProjectionMethod <: GeometricMethod end
GeometricIntegratorsBase.default_projection(::WrongProjectionMethod) = nothing
struct WrongExtrapolationMethod <: GeometricMethod end
function GeometricIntegratorsBase.default_extrapolation(::WrongExtrapolationMethod)
    NoInitialGuess()
end
struct WrongInitmethodMethod <: GeometricMethod end
GeometricIntegratorsBase.initmethod(::WrongInitmethodMethod) = 0

struct AnyNullvectorTableauMethod <: GeometricMethod end
GeometricIntegratorsBase.nullvector(::AnyNullvectorTableauMethod) = 1.5
GeometricBase.tableau(::AnyNullvectorTableauMethod) = "a tableau"

struct HermiteIguessMethod <: GeometricMethod end
GeometricIntegratorsBase.default_iguess(::HermiteIguessMethod) = HermiteExtrapolation()

struct ThrowingDefaultsMethod <: GeometricMethod end
GeometricIntegratorsBase.nullvector(::ThrowingDefaultsMethod) = error("nullvector")
GeometricIntegratorsBase.default_solver(::ThrowingDefaultsMethod) = error("default_solver")

struct AllDefaultsThrowMethod <: GeometricMethod end
GeometricIntegratorsBase.initmethod(::AllDefaultsThrowMethod) = error("initmethod")
GeometricIntegratorsBase.nullvector(::AllDefaultsThrowMethod) = error("nullvector")
GeometricBase.tableau(::AllDefaultsThrowMethod) = error("tableau")
GeometricIntegratorsBase.default_solver(::AllDefaultsThrowMethod) = error("default_solver")
GeometricIntegratorsBase.default_iguess(::AllDefaultsThrowMethod) = error("default_iguess")
function GeometricIntegratorsBase.default_projection(::AllDefaultsThrowMethod)
    error("default_projection")
end
function GeometricIntegratorsBase.default_extrapolation(::AllDefaultsThrowMethod)
    error("default_extrapolation")
end

struct InterruptingDefaultMethod <: GeometricMethod end
function GeometricIntegratorsBase.default_iguess(::InterruptingDefaultMethod)
    throw(InterruptException())
end

struct InterruptingTraitMethod <: GeometricIntegratorsBase.ODEMethod end
function GeometricIntegratorsBase.isodemethod(::InterruptingTraitMethod)
    throw(InterruptException())
end

struct MyODE <: GeometricIntegratorsBase.ODEMethod end
GeometricIntegratorsBase.ispodemethod(::MyODE) = true
GeometricIntegratorsBase.isiodemethod(::MyODE) = true

struct MyPODE <: GeometricIntegratorsBase.PODEMethod end
GeometricIntegratorsBase.ispodemethod(::MyPODE) = false

struct MissingTraitHODE <: GeometricIntegratorsBase.HODEMethod end
GeometricIntegratorsBase.ishodemethod(::MissingTraitHODE) = missing

struct ThrowingTraitLODE <: GeometricIntegratorsBase.LODEMethod end
GeometricIntegratorsBase.islodemethod(::ThrowingTraitLODE) = error("islodemethod")

# One method per family that denies its own family.
const FAMILY_TRAITS = (
    :ODEMethod => :isodemethod,
    :PODEMethod => :ispodemethod,
    :HODEMethod => :ishodemethod,
    :IODEMethod => :isiodemethod,
    :LODEMethod => :islodemethod,
    :SODEMethod => :issodemethod,
    :DAEMethod => :isdaemethod,
    :PDAEMethod => :ispdaemethod,
    :HDAEMethod => :ishdaemethod,
    :IDAEMethod => :isidaemethod,
    :LDAEMethod => :isldaemethod,
    :DELEMethod => :isdelemethod
)
for (family, trait) in FAMILY_TRAITS
    T = Symbol(:Denying, family)
    @eval struct $T <: GeometricIntegratorsBase.$family end
    @eval GeometricIntegratorsBase.$trait(::$T) = false
end

struct NoFamilyMethod <: GeometricMethod end
struct NoFamilyDeterministicMethod <: GeometricIntegratorsBase.DeterministicMethod end

@testset "the checks of AbstractMethod run first" begin
    v = GeometricBase.interface_violations(ThrowingPropertyMethod())
    @test fk(v) == [(:issymplectic, :throws)]
    v = GeometricBase.interface_violations(ThrowingPropertyWrongSolverMethod())
    @test fk(v) == [(:issymplectic, :throws), (:default_solver, :type)]
end

@testset "include_missing is passed on" begin
    properties = (
        :isexplicit, :isimplicit, :issymmetric, :issymplectic, :isenergypreserving,
        :isstifflyaccurate, :order, :name, :description, :reference)
    v = GeometricBase.interface_violations(NoPropertyMethod(); include_missing = true)
    @test fk(v) == [(f, :missing) for f in properties]
    @test isempty(GeometricBase.interface_violations(NoPropertyMethod()))
    @test isempty(GeometricBase.interface_violations(NoPropertyMethod(); include_missing = false))
end

@testset "the type of a default" begin
    for (M, f) in (
        WrongSolverMethod => :default_solver,
        WrongIguessMethod => :default_iguess,
        WrongProjectionMethod => :default_projection,
        WrongExtrapolationMethod => :default_extrapolation,
        WrongInitmethodMethod => :initmethod
    )
        v = GeometricBase.interface_violations(M())
        @test fk(v) == [(f, :type)]
    end
    @test only(GeometricBase.interface_violations(WrongSolverMethod())).detail ==
          "answers a value of type Int64"
    # `nullvector` and `tableau` have no type rule.
    @test isempty(GeometricBase.interface_violations(AnyNullvectorTableauMethod()))
    @test isempty(GeometricBase.interface_violations(TestMethod()))
end

@testset "an extrapolation is an initial guess" begin
    @test isempty(GeometricBase.interface_violations(HermiteIguessMethod()))
end

@testset "a default that throws" begin
    v = GeometricBase.interface_violations(ThrowingDefaultsMethod())
    @test fk(v) == [(:nullvector, :throws), (:default_solver, :throws)]
    @test v[1].detail == "throws ErrorException"
    defaults = (:initmethod, :nullvector, :tableau, :default_solver, :default_iguess,
        :default_projection, :default_extrapolation)
    v = GeometricBase.interface_violations(AllDefaultsThrowMethod())
    @test fk(v) == [(f, :throws) for f in defaults]
end

@testset "an interrupt is thrown again" begin
    # An `InterruptException` carries no message, and `@test_throws` catches one only when its
    # first argument is the type itself.
    # ast-grep-ignore: test-throws-bare-type
    @test_throws InterruptException GeometricBase.interface_violations(InterruptingDefaultMethod())
    # ast-grep-ignore: test-throws-bare-type
    @test_throws InterruptException GeometricBase.interface_violations(InterruptingTraitMethod())
end

@testset "a method for more kinds of problem" begin
    @test ispodemethod(MyODE()) && isiodemethod(MyODE())
    @test isempty(GeometricBase.interface_violations(MyODE()))
    @test isempty(GeometricBase.interface_violations(ImplicitMidpoint()))
end

@testset "a family trait that does not match" begin
    v = GeometricBase.interface_violations(MyPODE())
    @test fk(v) == [(:ispodemethod, :value)]
    @test occursin("PODEMethod", only(v).detail)
    for (family, trait) in FAMILY_TRAITS
        M = getfield(@__MODULE__, Symbol(:Denying, family))
        v = GeometricBase.interface_violations(M())
        @test fk(v) == [(trait, :value)]
        @test occursin(string(family), only(v).detail)
    end
    @test fk(GeometricBase.interface_violations(MissingTraitHODE())) ==
          [(:ishodemethod, :value)]
    @test fk(GeometricBase.interface_violations(ThrowingTraitLODE())) ==
          [(:islodemethod, :throws)]
end

@testset "a method under no family" begin
    @test isempty(GeometricBase.interface_violations(NoFamilyMethod()))
    @test isempty(GeometricBase.interface_violations(NoFamilyDeterministicMethod()))
    @test isempty(GeometricBase.interface_violations(NoProjection()))
end

@testset "the methods of this package conform" begin
    for m in (ExplicitEuler(), ImplicitEuler(), ImplicitMidpoint(), CrankNicolson(),
        SymplecticEulerA(), SymplecticEulerB(), NoProjection())
        @test isempty(GeometricBase.interface_violations(m))
    end
end
