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

```math
 =  =  .
\begin{aligned}
ε_{B_{1g}} &= ε_{xx} - ε_{yy}
\\
&= \partial_x u_x - \partial_y u_y
\\
&= -\frac{ (\mathbf{b} \cdot \mathbf{r}) \cdot 2xy }{2πr^4}.
\end{aligned}
```

# Examples
```jldocstest
julia> b1g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
-0.07957747154594767
```
"""
@inline function b1g_shear( eval_r::Vector2D, bob::Vector2D )
    return -(bob ⋅ eval_r) * ( 2 * xcomponent(eval_r) * ycomponent(eval_r) ) / ( 2π * magnitude2(eval_r)^2 )
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
\begin{aligned}
ε_{B_{2g}} &= 2ε_{xy}
\\
&= \partial_y u_x + \partial_x u_y
\\
&= \frac{(\mathbf{b} \cdot \mathbf{r}) \cdot (x^2 - y^2) }{2πr^4}.
\end{aligned}
```

# Examples
```jldocstest
julia> b2g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
0.15915494309189535
```
"""
@inline function b2g_shear( eval_r::Vector2D, bob::Vector2D )
    return (bob ⋅ eval_r) * ( xcomponent(eval_r)^2 - ycomponent(eval_r)^2 ) / ( 2π * magnitude2(eval_r)^2 )
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
    Δsplitting(ε₁, ε₂, [ratio = oneunit(ε₂)])

Calculate the level-splitting Δ as the quadrature sum of the strains.

!!! note 
    In principle, there is no reason for the two strain fields to have a 
    equal couplings. The optional paramter `ratio = λ₂ / λ₁` represents the 
    relative strength of the elastic couplings `λ₁` and `λ₂`.  
"""
@inline Δsplitting(ε₁, ε₂, ratio = oneunit(ε₂)) = hypot(ε₁, ratio * ε₂)
    
end # module ShearFunctions