@doc raw"""
# Nonlinear Test Problems

Every problem in [`HarmonicOscillator`](@ref) and [`NonautonomousProblems`](@ref) is linear, and
every partitioned one is separable on top of that, that is ``v = v(t,p)`` and ``f = f(t,q)``.
Three things go untested on such problems: the nonlinear solve converges in a single Newton step,
methods that coincide on a linear problem are indistinguishable, and the off-diagonal blocks
``\partial v / \partial q`` and ``\partial f / \partial p`` of the residual Jacobian vanish, so
the two halves of the solver solution vector of a partitioned integrator are never actually
coupled. The problems collected here exercise all three.
"""
module NonlinearProblems

using GeometricEquations

export riccati_problem, riccati_error
export nonlinear_podeproblem, nonlinear_hodeproblem, nonlinear_pode_odeproblem


const t₀ = 0.0
const Δt = 0.1
const nt = 10
const timespan = (t₀, Δt * nt)

const q₀ = [0.5]
const p₀ = [0.3]


@doc raw"""
The Riccati equation
```math
\dot{q} = - q^{2} ,
```
whose exact solution is ``q(t) = 1 / (1 + t)`` for ``q(0) = 1``.

The ordinary counterpart of the partitioned problem below. It tells the midpoint rule and the
trapezoidal rule apart, which are the same map on a linear problem, and it makes the nonlinear
solve take more than the single Newton step a linear problem needs.
"""
riccati_v(v, t, q, params) = (v[1] = -q[1]^2; nothing)

riccati_problem(; timestep=Δt) = ODEProblem(riccati_v, (0.0, 1.0), timestep, [1.0])

"""
Maximum absolute error of a solution of [`riccati_problem`](@ref) against its exact solution.
"""
riccati_error(sol) = maximum(abs(sol.q[n][1] - 1 / (1 + sol.t[n])) for n in axes(sol.q, 1))


@doc raw"""
A nonlinear, non-separable partitioned system,
```math
\dot{q} = p + \tfrac{3}{10} \, q p ,
\qquad
\dot{p} = - \sin (q) - \tfrac{1}{5} \, p^{2} ,
```
in which the velocity depends on the position and the force depends on the momentum, so that the
residual Jacobian of a partitioned integrator couples the two components.

It is *not* a Hamiltonian system — the Hamiltonian formulation below exists only so that the
`HODEProblem` code path sees the same equation, and no test asserts an energy behaviour on it.
"""
function nonlinear_pode_v(v, t, q, p, params)
    v[1] = p[1] + 0.3 * q[1] * p[1]
    nothing
end

function nonlinear_pode_f(f, t, q, p, params)
    f[1] = -sin(q[1]) - 0.2 * p[1]^2
    nothing
end

# only needed to construct the `HODEProblem`, never evaluated by a test
nonlinear_hamiltonian(t, q, p, params) = p[1]^2 / 2 + (1 - cos(q[1]))

function nonlinear_podeproblem(q₀=q₀, p₀=p₀; timespan=timespan, timestep=Δt)
    PODEProblem(nonlinear_pode_v, nonlinear_pode_f, timespan, timestep, q₀, p₀)
end

function nonlinear_hodeproblem(q₀=q₀, p₀=p₀; timespan=timespan, timestep=Δt)
    HODEProblem(nonlinear_pode_v, nonlinear_pode_f, nonlinear_hamiltonian,
        timespan, timestep, q₀, p₀)
end


@doc raw"""
The first order system ``x = (q, p)`` that [`nonlinear_podeproblem`](@ref) is identical to,
```math
\dot{x}_{1} = x_{2} + \tfrac{3}{10} \, x_{1} x_{2} ,
\qquad
\dot{x}_{2} = - \sin (x_{1}) - \tfrac{1}{5} \, x_{2}^{2} .
```
A method that is implemented for partitioned and for ordinary differential equations alike has to
produce the very same map on the two. In contrast to the linear problems, agreement here also
pins down the coupled blocks of the residual and the iterates of the nonlinear solve.
"""
function nonlinear_pode_ode_v(v, t, x, params)
    v[1] = x[2] + 0.3 * x[1] * x[2]
    v[2] = -sin(x[1]) - 0.2 * x[2]^2
    nothing
end

function nonlinear_pode_odeproblem(q₀=q₀, p₀=p₀; timespan=timespan, timestep=Δt)
    ODEProblem(nonlinear_pode_ode_v, timespan, timestep, vcat(q₀, p₀))
end

end
