using GeometricIntegratorsBase
using Test

import GeometricBase
import ..HarmonicOscillator: odeproblem

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
@test GeometricBase.issymplectic(ImplicitMidpoint()) == true
@test GeometricBase.isexplicit(ExplicitEuler()) == true

# `print_reference` reads the reference off `GeometricBase.reference`, where every method in the
# ecosystem defines it.
GeometricBase.reference(::TestMethod) = "n/a"
@test sprint(GeometricIntegratorsBase.print_reference, TestMethod()) == "n/a"
