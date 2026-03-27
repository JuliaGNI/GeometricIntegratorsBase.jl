using Profile
# using ProfileView
using PProf
using GeometricIntegratorsBase
using GeometricIntegratorsBase: Solution, solutionstep, enforce_periodicity!, ntime

include("harmonic_oscillator.jl")

using ..HarmonicOscillator

const Δt = 0.1
const nt = 100_000
const t₀ = 0.0
const t₁ = nt * Δt

# method = ExplicitEuler()
# method = ImplicitEuler()


function timeit(method)
    ode = odeproblem(timespan=(t₀, t₁), timestep=Δt)
    int = GeometricIntegrator(ode, method)
    sol = Solution(ode)

    @time integrate!(sol, int)
    @time integrate!(sol, int)
end

# timeit(ExplicitEuler())
# timeit(ImplicitEuler())


function profile(method)
    ode = odeproblem(timespan=(t₀, t₁), timestep=Δt)
    int = GeometricIntegrator(ode, method)
    sol = Solution(ode)

    @time integrate!(sol, int)
    @time integrate!(sol, int)

    Profile.clear()
    Profile.clear_malloc_data()

    Profile.Allocs.@profile integrate!(sol, int)
end

# profile(ExplicitEuler())
# profile(ImplicitEuler())


function profilesomemore(method)
    ode = odeproblem(timespan=(t₀, t₁), timestep=Δt)
    int = GeometricIntegrator(ode, method)
    sol = Solution(ode)
    solstep = SolutionStep(ode)
    curstate = current(solstep)

    integrate!(sol, int)
    @time integrate!(sol, int)

    solstep = solutionstep(int, sol[0])
    @time integrate!(solstep, int)

    solstep = solutionstep(int, sol[0])
    integrate!(sol, int, 1, ntime(sol), solstep, curstate)
    solstep = solutionstep(int, sol[0])
    @time integrate!(sol, int, 1, ntime(sol), solstep, curstate)

    copy!(sol, curstate, 1)
    @time copy!(sol, curstate, 1)

    copy!(sol, current(solstep), 1)
    @time copy!(sol, current(solstep), 1)

    isnan(curstate)
    @time isnan(curstate)

    isnan(current(solstep))
    @time isnan(current(solstep))

    @time enforce_periodicity!(solstep)

    solstep = solutionstep(int, sol[0])

    Profile.clear()
    Profile.clear_malloc_data()

    # Profile.Allocs.@profile integrate!(sol, int)
    # Profile.Allocs.@profile integrate!(solstep, int)
    # Profile.Allocs.@profile integrate!(sol, int, 1, ntime(sol), solstep, curstate);
    # Profile.Allocs.@profile copy!(sol, curstate, 1)
    Profile.Allocs.@profile copy!(sol, current(solstep), 1)
    # Profile.Allocs.@profile isnan(current(solstep))

    # prof = Profile.Allocs.fetch()
    # PProf.Allocs.pprof(prof; from_c=false)
end

profilesomemore(ExplicitEuler())
profilesomemore(ImplicitEuler())


function profileview(method)
    ode = odeproblem(timespan=(t₀, t₁), timestep=Δt)
    int = GeometricIntegrator(ode, method)
    sol = Solution(ode)

    integrate!(sol, int)

    Profile.clear()
    Profile.clear_malloc_data()

    Profile.Allocs.@profile integrate!(sol, int)
end

# @profview profileview(ExplicitEuler())
# @profview profileview(ImplicitEuler())




# results = Profile.Allocs.fetch()
# allocs = results.allocs
# t = 0
# for r in allocs
#     st = r.stacktrace
#     for s in st
#         fn = string(s.file)
#         println("$(basename(fn)):$(s.line) $(r.size)")
#         t += r.size
#         break
#     end
# end
# println("total: $t")
