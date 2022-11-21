module ShearFunctions

using StaticArrays
using ..PhysicalVectors
using ..CrystalDefects
using ..PBCImageFields

export b1g_shear, b2g_shear, bxg_shears, Δsplitting

@doc raw"""
    b1g_shear( eval_r::Vector2D, bob::Vector2D )
    b1g_shear( position::Vector2D, dis::Dislocation )

Calculate the ``B_{1g}`` from an edge dislocation for a given Burgers vector `bob` at site `eval_r`. 

# Additional Information
* The edge dislocation is aligned along the ``z`` axis at the origin.
* This is in units of ``1/(1-σ)`` for a Poisson ratio ``σ``.
* This quantity is obtained for the 2-Fe crystallographic unit cell, however, it is computed within the 1-Fe basis.

```math
ε_{B_{1g}}^{(2 - {\rm Fe})} = \partial_y u_x + \partial_x u_y = ε_{xy}^{(1-{\rm Fe})}.
```

# Examples
```jldocstest
julia> b1g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
-0.07957747154594767
```
"""
@inline function b1g_shear( eval_r::Vector2D, bob::Vector2D )
    return 1/(2π) * ( ( xcomponent(eval_r)^2 - ycomponent(eval_r)^2 ) / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end
b1g_shear(position, dis::Dislocation) = b1g_shear( position - dislocationorigin(dis), burgersvector(dis) )

@doc raw"""
    b2g_shear( eval_r::Vector2D, bob::Vector2D )
    b2g_shear( position::Vector2D, dis::Dislocation )

Calculate the ``B_{2g}`` from an edge dislocation for a given Burgers vector `bob` at site `eval_r`. 

# Additional Information
* The edge dislocation is aligned along the ``z`` axis at the origin.
* This is in units of ``1/(1-σ)`` for a Poisson ratio ``σ``.
* This quantity is obtained for the 2-Fe crystallographic unit cell, however, it is computed within the 1-Fe basis.

```math
ε_{B_{2g}}^{(2 - {\rm Fe})} = \frac{1}{2} \left(  \partial_x u_x - \partial_y u_y  \right) = \frac{1}{2} \left( ε_{xx}^{(1-{\rm Fe})} - ε_{yy}^{(1-{\rm Fe})} \right).
```

# Examples
```jldocstest
julia> b2g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
0.15915494309189535
```
"""
@inline function b2g_shear( eval_r::Vector2D, bob::Vector2D )
    return -1/(2π) * ( xcomponent(eval_r) * ycomponent(eval_r) / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end
b2g_shear(position, dis::Dislocation) = b2g_shear( position - dislocationorigin(dis), burgersvector(dis) )

"""
    bxg_shears(eval_r::Vector2D, bob::Vector2D; diff::Function, source_r::Vector2D = Vector2D(0.,0.)) -> Tuple{T, T} where {T <: Real}
    bxg_shears(position::Vector2D, dis::Dislocation; diff::Function) -> Tuple{T, T} where {T <: Real}

Calculate the two shears `(b1g, b2g)` from the ``B_{1g}`` and ``B_{2g}`` channels from and edge dislocation situated at the `source_r`
location and aligned along the ``z`` axis with a Burger's vector `bob`. 

# Additional Information
* The keyword argument `diff <: Function` should be used in cases where `Vector2D` subtraction has a different definition than the normal one expected, for example with periodic boundary conditions.

# Examples
```jldocstest
julia> using StaticArrays

julia> A = MVector{2}(0., 0.)
2-element MVector{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> bxg_shears!( A, Vector2D(1.,1.), Vector2D(1.,0.) )

julia> A
2-element MVector{2, Float64} with indices SOneTo(2):
 -0.07957747154594767
  0.15915494309189535
```
"""
function bxg_shears( eval_r::Vector2D, bob::Vector2D; diff::Function = Base.:(-), source_r::Vector2D = Vector2D(0.,0.) )
    disp = diff(eval_r, source_r)
    b1g = b1g_shear(disp, bob)
    b2g = b2g_shear(disp, bob)
    return (b1g, b2g)
end
bxg_shears(position, dis::Dislocation; diff = Base.:(-)) = bxg_shears( position, burgersvector(dis); diff = diff, source_r = dislocationorigin(dis) )

"""
    Δsplitting(ε₁, ε₂)

Calculate the level-splitting Δ as the quadrature sum of the strains.
"""
@inline Δsplitting(ε₁, ε₂) = hypot(ε₁, ε₂)
# Δsplitting(ε₁, ε₂) = sqrt(ε₁^2 + ε₂^2)
    
end # module ShearFunctions