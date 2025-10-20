# Extrapolation Methods

The extrapolation routines are exclusively used for computing
initial guesses and are usually not called directly by the user.

```@autodocs
Modules = [GeometricIntegratorsBase]
Order   = [:constant, :type, :macro, :function]
Pages   = ["extrapolation/extrapolation.jl",
           "extrapolation/aitken_neville.jl",
           "extrapolation/euler.jl",
           "extrapolation/hermite.jl",
           "extrapolation/midpoint.jl"]
```
