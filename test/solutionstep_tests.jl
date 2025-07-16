using GeometricIntegratorsBase
using GeometricEquations
using Test

using GeometricBase: VectorfieldVariable, periodic
using GeometricIntegratorsBase: enforce_periodicity!

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

    @test solstep.q ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol=1E-14
    @test solstep.q̄ ≈ [0.0, 0.0] atol=1E-14

    reset!(solstep, Δt)

    @test solstep.q ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol=1E-14
    @test solstep.q̄ ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol=1E-14

    update!(solstep, (q=[3π, 0.0],))

    @test solstep.q ≈ [initial_conditions(ode).q[1] + 3π, initial_conditions(ode).q[2]] atol=1E-14
    @test solstep.q̄ ≈ [initial_conditions(ode).q[1], initial_conditions(ode).q[2]] atol=1E-14

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
