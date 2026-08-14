using SafeTestsets

# the test problems, included at top level so that every test set below can reach the modules
# they define as `..HarmonicOscillator` and so on
include("examples/harmonic_oscillator.jl")
include("examples/nonautonomous.jl")
include("examples/nonlinear.jl")

@safetestset "Initial Guess Tests                                                             " begin
    include("initialguess_tests.jl")
end
@safetestset "Method Tests                                                                    " begin
    include("method_tests.jl")
end
@safetestset "Integrator Cache Tests                                                          " begin
    include("cache_tests.jl")
end
@safetestset "Solution Step Tests                                                             " begin
    include("solutionstep_tests.jl")
end
@safetestset "Solver Tests                                                                    " begin
    include("solver_tests.jl")
end
@safetestset "Extrapolation Tests                                                             " begin
    include("extrapolation_tests.jl")
end
@safetestset "Integrator Tests                                                                " begin
    include("integrator_tests.jl")
end
@safetestset "Projection Tests                                                                " begin
    include("projection_tests.jl")
end
@safetestset "Common Integrator Tests                                                         " begin
    include("integrators/common_tests.jl")
end
@safetestset "ExplicitEuler Tests                                                             " begin
    include("integrators/explicit_euler_tests.jl")
end
@safetestset "ImplicitEuler Tests                                                             " begin
    include("integrators/implicit_euler_tests.jl")
end
@safetestset "SymplecticEuler Tests                                                           " begin
    include("integrators/symplectic_euler_tests.jl")
end
@safetestset "ImplicitMidpoint Tests                                                          " begin
    include("integrators/implicit_midpoint_tests.jl")
end
@safetestset "CrankNicolson Tests                                                             " begin
    include("integrators/crank_nicolson_tests.jl")
end
