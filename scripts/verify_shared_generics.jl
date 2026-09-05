# Verifies the claims the 0.6.7 CHANGELOG entry makes about the interface generics: that the
# eleven names this package defines methods on *are* the upstream generics rather than private
# functions of the same name, that none of the eleven imports is stale, and that
# `GeometricBase.reference` gains exactly one method from this package.
#
#     julia --startup-file=no scripts/verify_shared_generics.jl
#
# The measurement needs `GeometricIntegrators` loaded, because that is what pulls in `RungeKutta`
# and with it the reference strings on the tableau types. `GeometricIntegrators` is downstream of
# this package, so the script builds a temporary environment with this tree developed into it
# rather than using the repository's own — a probe environment inside the repository would resolve
# the registered version instead of the working tree.

using Pkg

Pkg.activate(; temp = true)
Pkg.develop(; path = dirname(@__DIR__))
# The three upstreams are added explicitly and not left to resolve transitively: only a direct
# dependency of the active project can be reached with `using`, and reading a generic off
# `GeometricBase` is the entire point of the script.
Pkg.add(["GeometricIntegrators", "GeometricBase", "GeometricEquations", "SimpleSolvers"])

using GeometricBase
using GeometricEquations
using GeometricIntegrators
using GeometricIntegratorsBase
using SimpleSolvers

const GIB = GeometricIntegratorsBase

const SHARED = [(:isexplicit, GeometricBase),
    (:isimplicit, GeometricBase),
    (:issymmetric, GeometricBase),
    (:issymplectic, GeometricBase),
    (:isenergypreserving, GeometricBase),
    (:isstifflyaccurate, GeometricBase),
    (:reference, GeometricBase),
    (:problem, GeometricEquations),
    (:cache, SimpleSolvers),
    (:initialize!, SimpleSolvers),
    (:method, SimpleSolvers)]

failures = String[]

for (name, upstream) in SHARED
    ours = getglobal(GIB, name)
    theirs = getglobal(upstream, name)
    if ours !== theirs
        push!(failures, "$name is not $upstream's function")
        continue
    end
    if !any(m -> m.module === GIB, methods(theirs))
        push!(failures, "$name carries no method from this package, so the import is stale")
    end
end

# The absolute counts are reported and not asserted: they move with any release of `RungeKutta` or
# `GeometricIntegrators` that adds a property, which is not this package's business. What is
# asserted is the decomposition behind the entry's "47 methods to 48" — this package contributes
# exactly one method to `GeometricBase.reference`, so striking it leaves the 47 that were already
# there, which is what the pre-fix state was.
reference_methods = methods(GeometricBase.reference)
reference_ours = count(m -> m.module === GIB, reference_methods)

if reference_ours != 1
    push!(failures,
        "GeometricBase.reference carries $reference_ours methods from this " *
        "package, expected exactly 1")
end

isexplicit_methods = methods(GeometricBase.isexplicit)

println()
println("GeometricBase.reference  : ", length(reference_methods), " methods (",
    reference_ours, " here, ", length(reference_methods) - reference_ours, " upstream)")
println("GeometricBase.isexplicit : ", length(isexplicit_methods), " methods (",
    count(m -> m.module === GIB, isexplicit_methods), " here)")
println("shared generics checked  : ", length(SHARED))

if isempty(failures)
    println("\nOK")
else
    println("\nFAILED")
    foreach(f -> println("  - ", f), failures)
    exit(1)
end
