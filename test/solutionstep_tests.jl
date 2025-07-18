using GeometricIntegratorsBase
using GeometricEquations
using Test

using GeometricBase: StateWithError, TimeVariable, VectorfieldVariable
using GeometricBase: periodic, value
using GeometricIntegratorsBase: enforce_periodicity!, internal, nhistory
using GeometricIntegratorsBase: _strip_symbol, _strip_bar, _strip_dot, _add_symbol, _add_bar, _add_dot, _state, _vectorfield

using ..HarmonicOscillator


Δt = 0.1
t0 = 0.0
t1 = t0 + Δt
x0 = rand(2)
q0 = rand(1)
p0 = q0 .^ 2
λ0 = rand(1)
λ1 = rand(1)
μ0 = rand(1)
μ1 = rand(1)
Δx = rand(2)
Δq = rand(1)
Δp = rand(1)
ΔW = rand(3)
ΔZ = rand(3)

ode = odeproblem()
dae = daeproblem()
pode = podeproblem()
pdae = pdaeproblem()
hode = hodeproblem()
hdae = hdaeproblem()
iode = iodeproblem()
idae = idaeproblem()
lode = lodeproblem()
ldae = ldaeproblem()


@testset "$(rpad("Solution Step Helper Functions",80))" begin
    x = rand(3)

    @test _strip_symbol(:q̄, Char(0x0304)) == :q
    @test _strip_symbol(:q̇, Char(0x0307)) == :q

    @test _strip_bar(:q̄) == :q
    @test _strip_dot(:q̇) == :q

    @test _add_symbol(:q, Char(0x0304)) == :q̄
    @test _add_symbol(:q, Char(0x0307)) == :q̇

    @test _add_bar(:q) == :q̄
    @test _add_dot(:q) == :q̇

    @test _state(1) == 0
    @test _state(StateVariable(x)) == StateWithError(StateVariable(zero(x)))

    @test ismissing(_vectorfield(1))
    @test ismissing(_vectorfield(TimeVariable(1.0)))
    @test ismissing(_vectorfield(AlgebraicVariable(x)))
    @test ismissing(_vectorfield(VectorfieldVariable(StateVariable(x))))

    @test _vectorfield(StateVariable(x)) == VectorfieldVariable(zero(x))
    @test _vectorfield(StateWithError(StateVariable(x))) == VectorfieldVariable(zero(x))

end


@testset "$(rpad("Solution Step Constructors",80))" begin
    @test typeof(SolutionStep(ode)) <: SolutionStep{ODE}
    @test typeof(SolutionStep(dae)) <: SolutionStep{DAE}
    @test typeof(SolutionStep(pode)) <: SolutionStep{PODE}
    @test typeof(SolutionStep(hode)) <: SolutionStep{HODE}
    @test typeof(SolutionStep(iode)) <: SolutionStep{IODE}
    @test typeof(SolutionStep(lode)) <: SolutionStep{LODE}
    @test typeof(SolutionStep(pdae)) <: SolutionStep{PDAE}
    @test typeof(SolutionStep(hdae)) <: SolutionStep{HDAE}
    @test typeof(SolutionStep(idae)) <: SolutionStep{IDAE}
    @test typeof(SolutionStep(ldae)) <: SolutionStep{LDAE}
end


@testset "$(rpad("SolutionStep Accessor Functions",80))" begin

    # Test with ODE solution step
    solstep_ode = SolutionStep(ode; nhistory=3)

    # Test nhistory
    @test nhistory(solstep_ode) == 3

    # Test basic accessor functions
    @test solution(solstep_ode) == solstep_ode.solution
    @test vectorfield(solstep_ode) == solstep_ode.vectorfield
    @test history(solstep_ode) == solstep_ode.history
    @test internal(solstep_ode) == solstep_ode.internal
    @test parameters(solstep_ode) == solstep_ode.parameters

    # Test keys
    @test keys(solstep_ode) == (:t, :q)

    # Test indexed accessor functions
    sol_current = solution(solstep_ode, 0)
    sol_previous = solution(solstep_ode, 1)

    @test haskey(sol_current, :t)
    @test haskey(sol_current, :q)
    @test haskey(sol_previous, :t)
    @test haskey(sol_previous, :q)

    # Test vectorfield accessor with index
    vf_current = vectorfield(solstep_ode, 0)
    vf_previous = vectorfield(solstep_ode, 1)

    @test haskey(vf_current, :q)
    @test haskey(vf_previous, :q)

    # Test history accessor with index
    hist_current = history(solstep_ode, 0)
    hist_previous = history(solstep_ode, 1)

    @test haskey(hist_current, :t)
    @test haskey(hist_current, :q)
    @test haskey(hist_current, :q̇)
    @test haskey(hist_previous, :t)
    @test haskey(hist_previous, :q)
    @test haskey(hist_previous, :q̇)

    # Test that indexed accessors return NamedTuples with correct structure
    @test typeof(sol_current) <: NamedTuple
    @test typeof(vf_current) <: NamedTuple
    @test typeof(hist_current) <: NamedTuple

    # Test that the returned values have the same keys as the original
    @test keys(sol_current) == keys(solution(solstep_ode))
    @test keys(vf_current) == keys(vectorfield(solstep_ode))
    @test keys(hist_current) == keys(history(solstep_ode))

    # Test with PODE solution step
    solstep_pode = SolutionStep(pode; nhistory=2)

    @test nhistory(solstep_pode) == 2
    @test keys(solstep_pode) == (:t, :q, :p)

    sol_pode_current = solution(solstep_pode, 0)
    @test haskey(sol_pode_current, :t)
    @test haskey(sol_pode_current, :q)
    @test haskey(sol_pode_current, :p)

    # Test with DAE solution step
    solstep_dae = SolutionStep(dae; nhistory=1)

    @test nhistory(solstep_dae) == 1
    @test keys(solstep_dae) == (:t, :q, :λ, :μ)

    sol_dae_current = solution(solstep_dae, 0)
    @test haskey(sol_dae_current, :t)
    @test haskey(sol_dae_current, :q)
    @test haskey(sol_dae_current, :λ)
    @test haskey(sol_dae_current, :μ)

end


@testset "$(rpad("ODE Solution Step",80))" begin

    solstep = SolutionStep(ode)

    @test hasproperty(solstep, :t)
    @test hasproperty(solstep, :q)
    @test hasproperty(solstep, :t̄)
    @test hasproperty(solstep, :q̄)
    @test hasproperty(solstep, :q̇)
    @test !hasproperty(solstep, :ṫ)

    # @test solstep == SolutionStep(initial_conditions(ode), parameters(ode))

    @test solstep.t == solution(solstep).t[0] == solution(solstep, 0).t == current(solstep).t
    @test solstep.q == solution(solstep).q[0] == solution(solstep, 0).q == current(solstep).q

    @test solstep.t̄ == solution(solstep).t[1] == solution(solstep, 1).t == previous(solstep).t
    @test solstep.q̄ == solution(solstep).q[1] == solution(solstep, 1).q == previous(solstep).q

    @test solstep.q̇ == vectorfield(solstep).q[0] == vectorfield(solstep, 0).q == current(solstep).q̇

    @test solstep.t == initial_conditions(ode).t
    @test solstep.q == initial_conditions(ode).q

    @test solstep.t̄ == zero(initial_conditions(ode).t)
    @test solstep.q̄ == zero(initial_conditions(ode).q)

    solstep.t = initial_conditions(ode).t
    solstep.q .= initial_conditions(ode).q

    solstep.t̄ = zero(initial_conditions(ode).t)
    solstep.q̄ .= zero(initial_conditions(ode).q)

    @test current(solstep) == (
        t=initial_conditions(ode).t,
        q=initial_conditions(ode).q,
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )
    @test previous(solstep) == (
        t=zero(initial_conditions(ode).t),
        q=zero(initial_conditions(ode).q),
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )

    reset!(solstep, Δt)

    @test solstep.t == initial_conditions(ode).t + Δt
    @test solstep.q == initial_conditions(ode).q
    @test solstep.t̄ == initial_conditions(ode).t
    @test solstep.q̄ == initial_conditions(ode).q

    @test current(solstep) == (
        t=initial_conditions(ode).t + Δt,
        q=initial_conditions(ode).q,
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )
    @test previous(solstep) == (
        t=initial_conditions(ode).t,
        q=initial_conditions(ode).q,
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )

    update!(solstep, (q=Δx,))

    @test solstep.t == initial_conditions(ode).t + Δt
    @test solstep.q == initial_conditions(ode).q .+ Δx
    @test solstep.t̄ == initial_conditions(ode).t
    @test solstep.q̄ == initial_conditions(ode).q

    @test current(solstep) == (
        t=initial_conditions(ode).t + Δt,
        q=initial_conditions(ode).q .+ Δx,
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )
    @test previous(solstep) == (
        t=initial_conditions(ode).t,
        q=initial_conditions(ode).q,
        q̇=VectorfieldVariable(initial_conditions(ode).q),
    )

    # test for periodicity treatment

    ode_periodic = ODEProblem(functions(ode).v, timespan(ode), timestep(ode), initial_conditions(ode).q; invariants=invariants(ode), parameters=parameters(ode), periodicity=([0, -Inf], [2π, +Inf]))

    solstep = SolutionStep(ode_periodic)

    @test solstep.q ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol = 1E-14
    @test solstep.q̄ ≈ [0.0, 0.0] atol = 1E-14

    reset!(solstep, Δt)

    @test solstep.q ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol = 1E-14
    @test solstep.q̄ ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol = 1E-14

    update!(solstep, (q=[3π, 0.0],))

    @test solstep.q ≈ [initial_conditions(ode).q[1] + 3π, initial_conditions(ode).q[2]] atol = 1E-14
    @test solstep.q̄ ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol = 1E-14

    enforce_periodicity!(solstep)

    @test solstep.q == [initial_conditions(ode).q[1] + π, initial_conditions(ode).q[2]]
    @test solstep.q̄ == [initial_conditions(ode).q[1] - 2π, initial_conditions(ode).q[2]]



    #     solstep = SolutionStep(ode)

    #     k = HarmonicOscillator.k
    #     ω = HarmonicOscillator.ω

    #     A = sqrt(solstep.q[2]^2 / k + solstep.q[1]^2)
    #     ϕ = asin(solstep.q[1] / A)

    #     @test solstep.t == initial_conditions(ode).t
    #     @test solstep.q == initial_conditions(ode).q

    #     @test solstep.t̄ ≈ - timestep(ode)
    #     @test solstep.q̄ ≈ [A * sin(- ω * timestep(ode) + ϕ), ω * A * cos(- ω * timestep(ode) + ϕ)]

end


@testset "$(rpad("PODE Solution Step",80))" begin

    solstep = SolutionStep(pode)

    @test hasproperty(solstep, :t)
    @test hasproperty(solstep, :q)
    @test hasproperty(solstep, :p)
    @test hasproperty(solstep, :t̄)
    @test hasproperty(solstep, :q̄)
    @test hasproperty(solstep, :p̄)
    @test hasproperty(solstep, :q̇)
    @test hasproperty(solstep, :ṗ)
    @test !hasproperty(solstep, :ṫ)

    # @test solstep == SolutionStep(initial_conditions(pode), parameters(pode))

    @test solstep.t == solution(solstep).t[0] == solution(solstep, 0).t == current(solstep).t
    @test solstep.q == solution(solstep).q[0] == solution(solstep, 0).q == current(solstep).q
    @test solstep.p == solution(solstep).p[0] == solution(solstep, 0).p == current(solstep).p

    @test solstep.t̄ == solution(solstep).t[1] == solution(solstep, 1).t == previous(solstep).t
    @test solstep.q̄ == solution(solstep).q[1] == solution(solstep, 1).q == previous(solstep).q
    @test solstep.p̄ == solution(solstep).p[1] == solution(solstep, 1).p == previous(solstep).p

    @test solstep.q̇ == vectorfield(solstep).q[0] == vectorfield(solstep, 0).q == current(solstep).q̇
    @test solstep.ṗ == vectorfield(solstep).p[0] == vectorfield(solstep, 0).p == current(solstep).ṗ

    @test solstep.t == initial_conditions(pode).t
    @test solstep.q == initial_conditions(pode).q
    @test solstep.p == initial_conditions(pode).p

    @test solstep.t̄ == zero(initial_conditions(pode).t)
    @test solstep.q̄ == zero(initial_conditions(pode).q)
    @test solstep.p̄ == zero(initial_conditions(pode).p)

    solstep.t = initial_conditions(pode).t
    solstep.q .= initial_conditions(pode).q
    solstep.p .= initial_conditions(pode).p

    solstep.t̄ = zero(initial_conditions(pode).t)
    solstep.q̄ .= zero(initial_conditions(pode).q)
    solstep.p̄ .= zero(initial_conditions(pode).p)

    @test current(solstep) == (
        t=initial_conditions(pode).t,
        q=initial_conditions(pode).q,
        p=initial_conditions(pode).p,
        q̇=VectorfieldVariable(initial_conditions(pode).q),
        ṗ=VectorfieldVariable(initial_conditions(pode).p),
    )
    @test previous(solstep) == (
        t=zero(initial_conditions(pode).t),
        q=zero(initial_conditions(pode).q),
        p=zero(initial_conditions(pode).p),
        q̇=VectorfieldVariable(initial_conditions(pode).q),
        ṗ=VectorfieldVariable(initial_conditions(pode).p),
    )

    reset!(solstep, Δt)
    @test current(solstep) == (
        t=initial_conditions(pode).t + Δt,
        q=initial_conditions(pode).q,
        p=initial_conditions(pode).p,
        q̇=VectorfieldVariable(initial_conditions(pode).q),
        ṗ=VectorfieldVariable(initial_conditions(pode).p),
    )
    @test previous(solstep) == (
        t=initial_conditions(pode).t,
        q=initial_conditions(pode).q,
        p=initial_conditions(pode).p,
        q̇=VectorfieldVariable(initial_conditions(pode).q),
        ṗ=VectorfieldVariable(initial_conditions(pode).p),
    )

    #     update!(solstep, Δq, Δp)
    #     @test solstep.t == t0  + Δt
    #     @test solstep.q == q0 .+ Δq
    #     @test solstep.p == p0 .+ Δp

    #     copy!(solstep, (t = t1, q = [2π], p = [2π]))
    #     @test solstep.t == t1
    #     @test solstep.q == [2π]
    #     @test solstep.p == [2π]

    #     cut_periodic_solution!(solstep, (q = [2π],))
    #     @test solstep.q == [0.]
    #     @test solstep.p == [2π]


    #     solstep = SolutionStep(pode)

    #     k = HarmonicOscillator.k
    #     ω = HarmonicOscillator.ω

    #     A = sqrt(solstep.p[1]^2 / k + solstep.q[1]^2)
    #     ϕ = asin(solstep.q[1] / A)

    #     @test solstep.t == initial_conditions(pode).t
    #     @test solstep.q == initial_conditions(pode).q
    #     @test solstep.p == initial_conditions(pode).p

    #     @test solstep.t̄ ≈ - timestep(pode)
    #     @test solstep.q̄ ≈ [A * sin(- ω * timestep(pode) + ϕ)]
    #     @test solstep.p̄ ≈ [ω * A * cos(- ω * timestep(pode) + ϕ)]

end


@testset "$(rpad("DAE Solution Step",80))" begin

    solstep = SolutionStep(dae)

    @test hasproperty(solstep, :t)
    @test hasproperty(solstep, :q)
    @test hasproperty(solstep, :λ)
    @test hasproperty(solstep, :t̄)
    @test hasproperty(solstep, :q̄)
    @test hasproperty(solstep, :λ̄)
    @test hasproperty(solstep, :q̇)
    @test !hasproperty(solstep, :ṫ)
    @test !hasproperty(solstep, :λ̇)
    @test !hasproperty(solstep, :μ̇)

    # @test solstep == SolutionStep(initial_conditions(ode), parameters(ode))

    @test solstep.t == solution(solstep).t[0] == solution(solstep, 0).t == current(solstep).t
    @test solstep.q == solution(solstep).q[0] == solution(solstep, 0).q == current(solstep).q
    @test solstep.λ == solution(solstep).λ[0] == solution(solstep, 0).λ == current(solstep).λ
    @test solstep.μ == solution(solstep).μ[0] == solution(solstep, 0).μ == current(solstep).μ

    @test solstep.t̄ == solution(solstep).t[1] == solution(solstep, 1).t == previous(solstep).t
    @test solstep.q̄ == solution(solstep).q[1] == solution(solstep, 1).q == previous(solstep).q
    @test solstep.λ̄ == solution(solstep).λ[1] == solution(solstep, 1).λ == previous(solstep).λ
    @test solstep.μ̄ == solution(solstep).μ[1] == solution(solstep, 1).μ == previous(solstep).μ

    @test solstep.q̇ == vectorfield(solstep).q[0] == vectorfield(solstep, 0).q == current(solstep).q̇

    @test solstep.t == initial_conditions(dae).t
    @test solstep.q == initial_conditions(dae).q
    @test solstep.λ == initial_conditions(dae).λ
    @test solstep.μ == initial_conditions(dae).μ

    @test solstep.t̄ == zero(initial_conditions(dae).t)
    @test solstep.q̄ == zero(initial_conditions(dae).q)
    @test solstep.λ̄ == zero(initial_conditions(dae).λ)
    @test solstep.μ̄ == zero(initial_conditions(dae).μ)

    solstep.t = initial_conditions(dae).t
    solstep.q .= initial_conditions(dae).q
    solstep.λ .= initial_conditions(dae).λ
    solstep.μ .= initial_conditions(dae).μ

    solstep.t̄ = zero(initial_conditions(dae).t)
    solstep.q̄ .= zero(initial_conditions(dae).q)
    solstep.λ̄ .= zero(initial_conditions(dae).λ)
    solstep.μ̄ .= zero(initial_conditions(dae).μ)

    @test current(solstep) == (
        t=initial_conditions(dae).t,
        q=initial_conditions(dae).q,
        λ=initial_conditions(dae).λ,
        μ=zero(initial_conditions(dae).μ),
        q̇=VectorfieldVariable(initial_conditions(dae).q),
    )
    @test previous(solstep) == (
        t=zero(initial_conditions(dae).t),
        q=zero(initial_conditions(dae).q),
        λ=zero(initial_conditions(dae).λ),
        μ=zero(initial_conditions(dae).μ),
        q̇=VectorfieldVariable(initial_conditions(dae).q),
    )

    reset!(solstep, Δt)
    @test current(solstep) == (
        t=initial_conditions(dae).t + Δt,
        q=initial_conditions(dae).q,
        λ=initial_conditions(dae).λ,
        μ=zero(initial_conditions(dae).μ),
        q̇=VectorfieldVariable(initial_conditions(dae).q),
    )
    @test previous(solstep) == (
        t=initial_conditions(dae).t,
        q=initial_conditions(dae).q,
        λ=initial_conditions(dae).λ,
        μ=zero(initial_conditions(dae).μ),
        q̇=VectorfieldVariable(initial_conditions(dae).q),
    )

    #     update!(solstep, Δx, AlgebraicVariable(λ1))
    #     @test solstep.t == t0  + Δt
    #     @test solstep.q == x0 .+ Δx
    #     @test solstep.λ == λ1
    #     @test solstep.μ == μ0

    #     update!(solstep, Δx, λ1, μ1)
    #     @test solstep.t == t0  + Δt
    #     @test solstep.q == x0 .+ Δx .+ Δx
    #     @test solstep.λ == λ1
    #     @test solstep.μ == μ1

    #     copy!(solstep, (t = t1, q = [-2π,2π], λ = λ1, μ = μ1))
    #     @test solstep.t == t1
    #     @test solstep.q == [-2π,2π]
    #     @test solstep.λ == λ1
    #     @test solstep.μ == μ1
    #     cut_periodic_solution!(solstep, (q = [2π, 0.],))
    #     @test solstep.q == [0., 2π]


    #     solstep = SolutionStep(dae)

    #     k = HarmonicOscillator.k
    #     ω = HarmonicOscillator.ω

    #     A = sqrt(solstep.q[2]^2 / k + solstep.q[1]^2)
    #     ϕ = asin(solstep.q[1] / A)

    #     @test solstep.t == initial_conditions(dae).t
    #     @test solstep.q == initial_conditions(dae).q

    #     @test solstep.t̄ ≈ - timestep(dae)
    #     @test solstep.q̄[1:2] ≈ [A * sin(- ω * timestep(dae) + ϕ), ω * A * cos(- ω * timestep(dae) + ϕ)]

end


# @testset "$(rpad("PDAE Solution Step",80))" begin

#     solstep = SolutionStepPDAE(t0, StateVariable(q0), StateVariable(p0), AlgebraicVariable(λ0), AlgebraicVariable(μ0), parameters(pdae))

#     @test solstep.t == solution(solstep).t[0] == current(solstep).t
#     @test solstep.q == solution(solstep).q[0] == current(solstep).q
#     @test solstep.p == solution(solstep).p[0] == current(solstep).p
#     @test solstep.λ == solution(solstep).λ[0] == current(solstep).λ
#     @test solstep.μ == solution(solstep).μ[0] == current(solstep).μ

#     @test solstep.t̄ == solution(solstep).t[1] == previous(solstep).t
#     @test solstep.q̄ == solution(solstep).q[1] == previous(solstep).q
#     @test solstep.p̄ == solution(solstep).p[1] == previous(solstep).p
#     @test solstep.λ̄ == solution(solstep).λ[1] == previous(solstep).λ
#     @test solstep.μ̄ == solution(solstep).μ[1] == previous(solstep).μ

#     solstep.t  = t0
#     solstep.q .= q0
#     solstep.p .= p0
#     solstep.λ .= λ0
#     solstep.μ .= μ0

#     solstep.t̄  = zero(t0)
#     solstep.q̄ .= zero(q0)
#     solstep.p̄ .= zero(p0)
#     solstep.λ̄ .= zero(λ0)
#     solstep.μ̄ .= zero(μ0)

#     @test current(solstep) == (t = t0, q = q0, p = p0, λ = λ0, μ = μ0)
#     @test previous(solstep) == (t = zero(t0), q = zero(q0), p = zero(p0), λ = zero(λ0), μ = zero(μ0))

#     reset!(solstep, Δt)
#     @test current(solstep) == (t = t0 + Δt, q = q0, p = p0, λ = λ0, μ = μ0)
#     @test previous(solstep) == (t = t0, q = q0, p = p0, λ = λ0, μ = μ0)

#     update!(solstep, Δq, Δp, λ1)
#     @test solstep.t == t0  + Δt
#     @test solstep.q == q0 .+ Δq
#     @test solstep.p == p0 .+ Δp
#     @test solstep.λ == λ1
#     @test solstep.μ == μ0

#     update!(solstep, Δq, Δp, λ1, μ1)
#     @test solstep.t == t0  + Δt
#     @test solstep.q == q0 .+ Δq .+ Δq
#     @test solstep.p == p0 .+ Δp .+ Δp
#     @test solstep.λ == λ1
#     @test solstep.μ == μ1

#     copy!(solstep, (t = t1, q = [2π], p = [2π], λ = λ1, μ = μ1))
#     @test solstep.t == t1
#     @test solstep.q == [2π]
#     @test solstep.p == [2π]
#     @test solstep.λ == λ1
#     cut_periodic_solution!(solstep, (q = [2π],))
#     @test solstep.q == [0.]
#     @test solstep.p == [2π]


#     solstep = SolutionStep(pdae)

#     k = HarmonicOscillator.k
#     ω = HarmonicOscillator.ω

#     A = sqrt(solstep.p[1]^2 / k + solstep.q[1]^2)
#     ϕ = asin(solstep.q[1] / A)

#     @test solstep.t == initial_conditions(pdae).t
#     @test solstep.q == initial_conditions(pdae).q
#     @test solstep.p == initial_conditions(pdae).p

#     @test solstep.t̄ ≈ - timestep(pdae)
#     @test solstep.q̄[1] ≈ A * sin(- ω * timestep(pdae) + ϕ)
#     @test solstep.p̄[1] ≈ ω * A * cos(- ω * timestep(pdae) + ϕ)

# end


# Create a simple equation type for testing
struct TestEquation <: GeometricEquation{Nothing,Nothing,Nothing} end

# Helper function to create a basic SolutionStep for testing
function create_test_solutionstep(; nhistory=2)
    # Create initial conditions with both regular and periodic variables
    q_val = [1.0, 2.0, 3.0]
    p_val = [0.5, 1.5, 2.5]

    # Create periodic ranges for q (component 1 is periodic from 0 to 2π)
    q_range_min = [0.0, -Inf, -Inf]
    q_range_max = [2π, Inf, Inf]
    q_periodic = BitArray([true, false, false])

    # Create state variables
    t = TimeVariable(0.0)
    q = StateVariable(q_val, (q_range_min, q_range_max), q_periodic)
    p = StateVariable(p_val, (fill(-Inf, 3), fill(Inf, 3)), BitArray(fill(false, 3)))

    ics = (t=t, q=q, p=p)

    return SolutionStep{TestEquation}(ics, NullParameters(); nhistory=nhistory)
end


@testset "SolutionStep Update and Periodicity Functions" begin
    @testset "update! function" begin
        @testset "Basic functionality" begin
            solstep = create_test_solutionstep()

            # Store original values
            orig_t = value(solstep.t)
            orig_q = copy(solstep.q)
            orig_p = copy(solstep.p)

            # Define increments
            Δt = 0.1
            Δq = [0.1, 0.2, 0.3]
            Δp = [0.05, 0.1, 0.15]

            # Apply reset
            reset!(solstep, Δt)

            # Check that time value was updated correctly
            @test value(solstep.t) ≈ orig_t + Δt

            # Apply update
            result = update!(solstep, (q=Δq, p=Δp))

            # Check that function returns the solution step
            @test result === solstep

            # Check that values were updated correctly
            @test solstep.q ≈ orig_q + Δq
            @test solstep.p ≈ orig_p + Δp

            # Check that history wasn't affected
            @test solution(solstep)[:q][1] ≈ orig_q
            @test solution(solstep)[:p][1] ≈ orig_p
        end

        @testset "Partial updates" begin
            solstep = create_test_solutionstep()

            orig_q = copy(solstep.q)
            orig_p = copy(solstep.p)

            # Update only q
            Δq = [0.1, 0.2, 0.3]
            update!(solstep, (q=Δq,))

            @test solstep.q ≈ orig_q + Δq
            @test solstep.p ≈ orig_p  # p should remain unchanged
        end

        @testset "Empty update" begin
            solstep = create_test_solutionstep()

            orig_q = copy(solstep.q)
            orig_p = copy(solstep.p)

            # Empty update should not change anything
            update!(solstep, NamedTuple())

            @test solstep.q ≈ orig_q
            @test solstep.p ≈ orig_p
        end

        @testset "Error conditions" begin
            solstep = create_test_solutionstep()

            # Test with invalid key
            @test_throws AssertionError update!(solstep, (invalid_key=[1.0, 2.0],))

            # Test with subset of invalid keys
            @test_throws AssertionError update!(solstep, (q=[0.1, 0.2, 0.3], invalid_key=[1.0, 2.0]))
        end
    end

    @testset "enforce_periodicity! functions" begin

        @testset "Non-periodic variables (no-op)" begin
            solstep = create_test_solutionstep()

            # Get a non-periodic variable (p has no periodic components)
            p_vector = solution(solstep)[:p]
            orig_p = copy(p_vector[0])

            # This should do nothing
            enforce_periodicity!(solstep, p_vector)

            @test p_vector[0] ≈ orig_p
        end

        @testset "Periodic variable enforcement" begin
            solstep = create_test_solutionstep()

            # Test with value below range
            @testset "Value below range" begin
                solstep = create_test_solutionstep()
                q_vector = solution(solstep)[:q]

                # Set first component (periodic) to below range
                q_vector[0][1] = -0.5  # Below 0

                # Store history values
                hist_val = copy(q_vector[1][1])

                enforce_periodicity!(solstep, q_vector)

                # Should be adjusted by adding 2π
                @test q_vector[0][1] ≈ -0.5 + 2π
                @test q_vector[1][1] ≈ hist_val + 2π  # History should be adjusted too

                # Non-periodic components should be unchanged
                @test q_vector[0][2] ≈ 2.0
                @test q_vector[0][3] ≈ 3.0
            end

            @testset "Value above range" begin
                solstep = create_test_solutionstep()
                q_vector = solution(solstep)[:q]

                # Set first component (periodic) to above range
                q_vector[0][1] = 7.0  # Above 2π ≈ 6.28

                # Store history values
                hist_val = copy(q_vector[1][1])

                enforce_periodicity!(solstep, q_vector)

                # Should be adjusted by subtracting 2π
                @test q_vector[0][1] ≈ 7.0 - 2π
                @test q_vector[1][1] ≈ hist_val - 2π  # History should be adjusted too
            end

            @testset "Multiple period adjustments" begin
                solstep = create_test_solutionstep()
                q_vector = solution(solstep)[:q]

                # Set value multiple periods outside range
                q_vector[0][1] = -4π  # Two periods below

                enforce_periodicity!(solstep, q_vector)

                # Should be adjusted by adding 4π to bring it into [0, 2π]
                @test q_vector[0][1] ≈ -4π + 4π
                @test 0 ≤ q_vector[0][1] ≤ 2π
            end

            @testset "Value already in range" begin
                solstep = create_test_solutionstep()
                q_vector = solution(solstep)[:q]

                # Set value within range
                original_val = π  # Within [0, 2π]
                q_vector[0][1] = original_val

                enforce_periodicity!(solstep, q_vector)

                # Should remain unchanged
                @test q_vector[0][1] ≈ original_val
            end
        end

        @testset "Full solution step periodicity" begin
            solstep = create_test_solutionstep()

            # Set periodic component out of range
            solution(solstep)[:q][0][1] = -0.5

            # Store original non-periodic values
            orig_q2 = copy(solution(solstep)[:q][0][2])
            orig_q3 = copy(solution(solstep)[:q][0][3])
            orig_p = copy(solution(solstep)[:p][0])

            # Apply periodicity to entire solution step
            enforce_periodicity!(solstep)

            # Periodic component should be adjusted
            @test solution(solstep)[:q][0][1] ≈ -0.5 + 2π

            # Non-periodic components should be unchanged
            @test solution(solstep)[:q][0][2] ≈ orig_q2
            @test solution(solstep)[:q][0][3] ≈ orig_q3
            @test solution(solstep)[:p][0] ≈ orig_p
        end

        @testset "StateVariableWithError periodicity" begin
            # Create a solution step with StateVariableWithError
            q_val = [1.0, 2.0, 3.0]
            q_range_min = [0.0, -Inf, -Inf]
            q_range_max = [2π, Inf, Inf]
            q_periodic = BitArray([true, false, false])

            q_state = StateVariable(q_val, (q_range_min, q_range_max), q_periodic)
            q_with_error = StateWithError(q_state)

            t = TimeVariable(0.0)
            ics = (q=q_with_error, t=t)

            solstep = SolutionStep{TestEquation}(ics, NullParameters(); nhistory=2)

            # Set periodic component out of range
            solution(solstep)[:q][0][1] = -0.5

            enforce_periodicity!(solstep)

            # Should be adjusted for StateVariableWithError too
            @test solution(solstep)[:q][0][1] ≈ -0.5 + 2π
        end

        @testset "Edge cases" begin
            solstep = create_test_solutionstep()

            @testset "Boundary values" begin
                q_vector = solution(solstep)[:q]

                # Test exact boundary values
                q_vector[0][1] = 0.0  # At lower bound
                enforce_periodicity!(solstep, q_vector)
                @test q_vector[0][1] ≈ 0.0

                q_vector[0][1] = 2π  # At upper bound
                enforce_periodicity!(solstep, q_vector)
                @test q_vector[0][1] ≈ 2π
            end

            @testset "Very small adjustments" begin
                q_vector = solution(solstep)[:q]

                # Test with value just slightly outside range
                q_vector[0][1] = 2π + 1e-14
                enforce_periodicity!(solstep, q_vector)
                @test 0 ≤ q_vector[0][1] ≤ 2π

                q_vector[0][1] = -1e-14
                enforce_periodicity!(solstep, q_vector)
                @test 0 ≤ q_vector[0][1] ≤ 2π
            end
        end
    end

    @testset "Integration tests" begin
        @testset "Update followed by periodicity enforcement" begin
            solstep = create_test_solutionstep()

            # Update with increment that takes periodic component out of range
            # Starting at q[1] = 1.0, add 6.0 to get 7.0 > 2π
            update!(solstep, (q=[6.0, 0.0, 0.0],))

            # Before periodicity enforcement
            @test solution(solstep)[:q][0][1] ≈ 7.0

            # Enforce periodicity
            enforce_periodicity!(solstep)

            # After periodicity enforcement
            @test solution(solstep)[:q][0][1] ≈ 7.0 - 2π
            @test 0 ≤ solution(solstep)[:q][0][1] ≤ 2π
        end

        @testset "Multiple update and periodicity cycles" begin
            solstep = create_test_solutionstep()

            # Simulate multiple time steps with updates and periodicity
            for i in 1:5
                update!(solstep, (q=[0.0, 1.5, 0.0],))  # Add 1.5 to periodic component
                enforce_periodicity!(solstep)

                # Should always be in range
                @test 0 ≤ solution(solstep)[:q][0][1] ≤ 2π
            end
        end
    end
end
