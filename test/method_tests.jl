using GeometricIntegratorsBase
using Test

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
