"""
Compute p(x) where p is the unique polynomial of degree length(xi),
such that p(x[i]) = y[i]) for all i.
Call with
```julia
aitken_neville!(x::AbstractVector, t::Real, ti::AbstractVector, xi::AbstractMatrix)
```
where
* `x`:  evaluation value
* `t`:  evaluation point
* `ti`: interpolation nodes
* `xi`: interpolation values
"""
function aitken_neville!(x::AbstractArray, t::TT, ti::AbstractVector{TT}, xi::AbstractVector) where {TT}
    length(ti) == length(xi) || throw(ArgumentError("ti and xi must have the same length"))

    for _xi in xi
        axes(x) == axes(_xi) || throw(ArgumentError("All xi entries must have the same axes as x"))
    end

    for j in eachindex(ti)
        for i in eachindex(ti)[begin:end-j]
            @. xi[i] = xi[i+1] + (xi[i] - xi[i+1]) * (ti[i+j] - t) / (ti[i+j] - ti[i])
        end
    end

    copyto!(x, xi[1])
end
