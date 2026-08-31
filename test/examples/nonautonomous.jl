@doc raw"""
# Non-Autonomous Test Problems

Every problem in [`HarmonicOscillator`](@ref) is autonomous, so none of them can tell the
stage times of an integrator apart: evaluating a vector field at ``t_{n}`` or at ``t_{n+1}``
gives the same result there. The problems collected here have explicitly time dependent
vector fields and pin those times down.
"""
module NonautonomousProblems

using GeometricEquations
using GeometricSolutions

export nonautonomous_odeproblem, nonautonomous_ode_solution
export nonautonomous_podeproblem, nonautonomous_hodeproblem
export nonautonomous_pode_odeproblem
export nonautonomous_iodeproblem, nonautonomous_lodeproblem
export nonautonomous_iode_odeproblem
export nonautonomous_hamiltonian, nonautonomous_lagrangian

const t₀ = 0.0
const Δt = 0.1
const nt = 10
const timespan = (t₀, Δt * nt)

const x₀ = [1.0]
const q₀ = [1.0]
const p₀ = [0.5]

@doc raw"""
The linear, non-autonomous scalar equation
```math
\dot{q} = \cos (t) \, q ,
```
whose exact solution is ``q(t) = q_{0} \exp ( \sin t - \sin t_{0} )``.

A method that evaluates the vector field at the wrong stage time loses an order of accuracy
on this problem, so the convergence order tests discriminate directly.
"""
function nonautonomous_ode_v(v, t, q, params)
    v[1] = cos(t) * q[1]
    nothing
end

function nonautonomous_odeproblem(x₀ = x₀; timespan = timespan, timestep = Δt)
    ODEProblem(nonautonomous_ode_v, timespan, timestep, x₀)
end

nonautonomous_ode_solution(t, q₀, t₀) = q₀ .* exp(sin(t) - sin(t₀))

function nonautonomous_ode_solution(prob::ODEProblem)
    sol = GeometricSolution(prob)
    for n in eachtimestep(sol)
        sol[n].q .= nonautonomous_ode_solution(sol[n].t, sol[0].q, sol[0].t)
    end
    return sol
end

@doc raw"""
A separable, non-autonomous Hamiltonian system with
```math
H (t, q, p) = \frac{p^{2}}{2} + (1 + t) \, \frac{q^{2}}{2} ,
```
so that ``v (t, p) = p`` and ``f (t, q) = - (1 + t) \, q``.

Separability makes it a legitimate problem for the symplectic Euler methods, while the time
dependence of ``f`` distinguishes the stage times of the A and B variants.
"""
function nonautonomous_pode_v(v, t, q, p, params)
    v[1] = p[1]
    nothing
end

function nonautonomous_pode_f(f, t, q, p, params)
    f[1] = -(1 + t) * q[1]
    nothing
end

function nonautonomous_hamiltonian(t, q, p, params)
    p[1]^2 / 2 + (1 + t) * q[1]^2 / 2
end

function nonautonomous_podeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    PODEProblem(nonautonomous_pode_v, nonautonomous_pode_f, timespan, timestep, q₀, p₀;
        invariants = (h = nonautonomous_hamiltonian,))
end

function nonautonomous_hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    HODEProblem(nonautonomous_pode_v, nonautonomous_pode_f, nonautonomous_hamiltonian,
        timespan, timestep, q₀, p₀)
end

@doc raw"""
The first order system ``x = (q, p)`` that [`nonautonomous_podeproblem`](@ref) is identical to,
```math
\dot{x}_{1} = x_{2} , \qquad \dot{x}_{2} = - (1 + t) \, x_{1} .
```
A method that is implemented for partitioned and for ordinary differential equations alike has to
produce the very same map on the two, which pins down every stage time of the partitioned variant
exactly, not just to leading order.
"""
function nonautonomous_pode_ode_v(v, t, x, params)
    v[1] = x[2]
    v[2] = -(1 + t) * x[1]
    nothing
end

function nonautonomous_pode_odeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    ODEProblem(nonautonomous_pode_ode_v, timespan, timestep, vcat(q₀, p₀))
end

@doc raw"""
The non-autonomous, regular Lagrangian
```math
L (t, q, v) = (1 + t) \, \frac{v^{2} - q^{2}}{2} ,
```
whose momentum map and force,
```math
\vartheta (t, q, v) = (1 + t) \, v , \qquad f (t, q, v) = - (1 + t) \, q ,
```
are *both* explicitly time dependent, so that a wrong stage time in either of them is detected.

Since ``\vartheta`` is regular, the implicit equation is equivalent to the partitioned system
``\dot{q} = p / (1+t)``, ``\dot{p} = - (1+t) \, q``, which is available as
[`nonautonomous_iode_odeproblem`](@ref) and serves as a reference solution.
"""
function nonautonomous_iode_ϑ(p, t, q, v, params)
    p[1] = (1 + t) * v[1]
    nothing
end

function nonautonomous_iode_f(f, t, q, v, params)
    f[1] = -(1 + t) * q[1]
    nothing
end

function nonautonomous_iode_g(g, t, q, v, λ, params)
    g[1] = λ[1]
    nothing
end

function nonautonomous_iode_v(v, t, q, p, params)
    v[1] = p[1] / (1 + t)
    nothing
end

function nonautonomous_lagrangian(t, q, v, params)
    (1 + t) * (v[1]^2 - q[1]^2) / 2
end

function nonautonomous_iode_ω(ω, t, q, params)
    ω[1, 1] = 0
    ω[1, 2] = -1
    ω[2, 1] = +1
    ω[2, 2] = 0
    nothing
end

function nonautonomous_iodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    IODEProblem(nonautonomous_iode_ϑ, nonautonomous_iode_f, nonautonomous_iode_g,
        timespan, timestep, q₀, p₀; v̄ = nonautonomous_iode_v)
end

function nonautonomous_lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    LODEProblem(nonautonomous_iode_ϑ, nonautonomous_iode_f, nonautonomous_iode_g,
        nonautonomous_iode_ω, nonautonomous_lagrangian,
        timespan, timestep, q₀, p₀; v̄ = nonautonomous_iode_v)
end

@doc raw"""
The first order system ``x = (q, p)`` that [`nonautonomous_iodeproblem`](@ref) is equivalent to,
```math
\dot{x}_{1} = \frac{x_{2}}{1 + t} , \qquad \dot{x}_{2} = - (1 + t) \, x_{1} .
```
For a regular momentum map the implicit midpoint and Crank-Nicolson methods applied to the
implicit equation are the same map as applied to this system, so integrating both and comparing
pins down every stage time of the implicit variants exactly, not just to leading order.
"""
function nonautonomous_iode_ode_v(v, t, x, params)
    v[1] = x[2] / (1 + t)
    v[2] = -(1 + t) * x[1]
    nothing
end

function nonautonomous_iode_odeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = Δt)
    ODEProblem(nonautonomous_iode_ode_v, timespan, timestep, vcat(q₀, p₀))
end

end
