using GeometricIntegratorsBase
using Test

import GeometricBase
import GeometricEquations
import GeometricSolutions
import SimpleSolvers
import Unicode

# A name this package *owns* while one of its dependencies owns a different function of the same
# name is a re-definition rather than an extension. The two generics never meet: a property set
# here is invisible to a caller that asks the dependency's function, and dispatch that looks like
# it should reach the dependency's methods raises a `MethodError` instead. Julia reports nothing:
# an unexported upstream name is never in scope to clash with, and a definition here wins over an
# exported one anyway — which is what the import block at the top of
# `src/GeometricIntegratorsBase.jl` is for.
#
# The dependency list is derived so that a new dependency is covered without editing this test.
# `identify_package` is what answers "is this module a direct dependency of `mod`", which is the
# right scope: a transitively loaded package is not one this module can shadow by accident. Testing
# instead whether `mod` *binds* the module's name would miss a dependency reached as
# `using Dep: name`, which is how this package reaches `Unicode`.
#
# Two things are outside the scan, both empty for this package and neither of them silent if that
# changes: a function declared in a *submodule*, because `names` lists the submodule itself and not
# its contents, and one defined in a package *extension*, because `identify_package` does not
# resolve an extension to a dependency. This package has neither, so the machinery to cover them
# would be guarding nothing.

function upstream_modules(mod::Module)
    ups = Module[Base, Core]
    for m in values(Base.loaded_modules)
        (m === mod || m === Base || m === Core) && continue
        Base.identify_package(mod, String(nameof(m))) === nothing && continue
        push!(ups, m)
    end
    unique(ups)
end

# Ownership is the whole test, and `parentmodule` is what decides it. Functions, macros and types
# all count, because a name is a name to `using` whatever kind of thing it names.
#
# A union has to be unwrapped before it is asked, and then skipped: `parentmodule` throws for one,
# and this package has several — `AbstractProblemIODE` and friends are `Union`s of parametric
# `EquationProblem`s, so they are `UnionAll`s whose body is a `Union` and the throw is reached
# through the `UnionAll` branch rather than the obvious one. Skipping them loses nothing. A union
# belongs to no single module, so there is no ownership to compare; and an alias of upstream types
# is correctly *not* owned here, which is what the unwrap already establishes for the parametric
# structs this package does own.
function is_owned(mod::Module, x)
    x isa Function && return parentmodule(x) === mod
    x isa Type || return false
    unwrapped = Base.unwrap_unionall(x)
    unwrapped isa DataType || return false
    parentmodule(unwrapped) === mod
end

function shadowed_generics(mod::Module)
    ups = upstream_modules(mod)
    scanned = Symbol[]
    shadows = Tuple{Symbol, Module}[]
    for n in names(mod; all = true)
        startswith(String(n), "#") && continue
        isdefined(mod, n) || continue
        is_owned(mod, getglobal(mod, n)) || continue
        push!(scanned, n)
        # Every upstream owner is reported rather than the first one found. A name can collide with
        # two dependencies at once — `issymmetric` collides with both `LinearAlgebra` and
        # `GeometricBase`, which own genuinely different functions — and this iterates a `Dict`, so
        # stopping at the first would name one of the two arbitrarily and a reader who imported it
        # would fix nothing. There is no same-function case to skip either: a binding owned by
        # `mod` and one owned by any other module are never the same object.
        for up in ups
            isdefined(up, n) && is_owned(up, getglobal(up, n)) && push!(shadows, (n, up))
        end
    end
    (scanned = sort!(scanned; by = String),
        shadows = sort!(shadows; by = t -> (String(t[1]), String(nameof(t[2])))))
end

# Both preconditions are asserted, because either one failing quietly would make the check below
# pass without looking at anything: `names(mod; all = true)` lists only what the module declares
# itself, and it excludes every `using`- and `import`-introduced binding, so the upstreams cannot
# be recovered from it.
#
# `Unicode` is in the list because it is reached as `using Unicode: normalize`, the binding form an
# earlier version of this check missed entirely. `Base` and `Core` are deliberately not asserted:
# `upstream_modules` seeds those itself, so asserting them would pin the seeding rather than the
# derivation.
@test upstream_modules(GeometricIntegratorsBase) ⊇
      [GeometricBase, GeometricEquations, GeometricSolutions, SimpleSolvers, Unicode]

result = shadowed_generics(GeometricIntegratorsBase)

# Named rather than counted. `scanned > 0` is satisfied by a scan that reached one binding, and a
# bare count would need editing every time the package gains a function. These five are declared in
# this package and imported from nowhere, so their absence means the scan missed what it is for.
@test result.scanned ⊇
      [:default_solver, :initmethod, :internal_variables, :nullvector, :print_reference]
@test result.shadows == Tuple{Symbol, Module}[]
