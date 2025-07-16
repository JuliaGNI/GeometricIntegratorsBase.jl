"""
A `ProjectionMethod` is an algorithm that is applied together with a geometric integrator
to enforce constraints which are not automatically satisfied by the integrator.
Examples include conservation of invariants like energy or the Dirac constraint in [`IODE`](@ref)s.
"""
abstract type ProjectionMethod <: GeometricMethod end

struct NoProjection <: ProjectionMethod end

projection(::GeometricMethod) = NoProjection()
