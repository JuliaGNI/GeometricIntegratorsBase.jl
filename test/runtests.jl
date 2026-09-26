using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Interface" include("integration/interface.jl")
    @safetestset "Method" include("method.jl")
    @safetestset "Integrator cache" include("cache.jl")
    @safetestset "Solution step" include("solutionstep.jl")
    @safetestset "Solvers" include("solvers.jl")
    @safetestset "Extrapolation" include("extrapolation/extrapolation.jl")
    @safetestset "Integrator" include("integrator.jl")
    @safetestset "Explicit Euler" include("integrators/explicit_euler.jl")
    @safetestset "Implicit Euler" include("integrators/implicit_euler.jl")
    @safetestset "Symplectic Euler" include("integrators/symplectic_euler.jl")
    @safetestset "Implicit midpoint" include("integrators/implicit_midpoint.jl")
    @safetestset "Crank-Nicolson" include("integrators/crank_nicolson.jl")
end
if "slow" in GROUPS
    @safetestset "Common integrator properties" include("integrators/common.jl")
end
