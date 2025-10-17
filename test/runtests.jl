using SafeTestsets

include("harmonic_oscillator.jl")

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
@safetestset "Euler Tests                                                                     " begin
    include("euler_tests.jl")
end
@safetestset "ImplicitEuler Tests                                                             " begin
    include("implicit_euler_tests.jl")
end
