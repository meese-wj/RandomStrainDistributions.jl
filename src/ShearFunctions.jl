module ShearFunctions

using StaticArrays
using ..PhysicalVectors
using ..CrystalDefects
using ..PBCImageFields

export b1g_shear, b2g_shear, bxg_shears, Δsplitting, b1gshear_PBC, b2g_shear_PBC, b2g_PBC_manual

@doc raw"""
    b1g_shear( eval_r::Vector2D, bob::Vector2D )
    b1g_shear( position::Vector2D, dis::Dislocation )

Calculate the ``B_{1g}`` from an edge dislocation for a given Burgers vector `bob` at site `eval_r`. 

# Additional Information
* The edge dislocation is aligned along the ``z`` axis at the origin.
* This is in units of ``1/(1-σ)`` for a Poisson ratio ``σ``.
* This quantity is obtained for the 2-Fe crystallographic unit cell, however, it is computed within the 1-Fe basis.

```math
ε_{B_{1g}}^{(2 - {\rm Fe})} = \frac{1}{2} \left(  \partial_y u_x + \partial_x u_y  \right) = ε_{xy}^{(1-{\rm Fe})}.
```

# Examples
```jldocstest
julia> b1g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
-0.07957747154594767
```
"""
function b1g_shear( eval_r::Vector2D, bob::Vector2D )
    return 1/(4π) * ( ( xcomponent(eval_r)^2 - ycomponent(eval_r)^2 ) / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end
b1g_shear(position, dis::Dislocation) = b1g_shear( position - DislocationOrigin(dis), BurgersVector(dis) )

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
function b2g_shear( eval_r::Vector2D, bob::Vector2D )
    return -1/(2π) * ( xcomponent(eval_r) * ycomponent(eval_r) / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end
b2g_shear(position, dis::Dislocation) = b2g_shear( position - DislocationOrigin(dis), BurgersVector(dis) )

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
bxg_shears(position, dis::Dislocation; diff = Base.:(-)) = bxg_shears( position, BurgersVector(dis); diff = diff, source_r = DislocationOrigin(dis) )

"""
    Δsplitting(ε₁, ε₂)

Calculate the level-splitting Δ as the quadrature sum of the strains.
"""
Δsplitting(ε₁, ε₂) = sqrt(ε₁^2 + ε₂^2)

function b1g_shear_PBC(position, dis::Dislocation, Lx, Ly = Lx, tolerance = sqrt(eps()))
    ϕfunc(r, origin) = let bob = BurgersVector(dis) 
        b1g_shear(r - origin, bob)
    end
    return PBCField(ϕfunc, position, DislocationOrigin(dis), Lx, Ly, tolerance)
end

function b2g_shear_PBC(position, dis::Dislocation, Lx, Ly = Lx, tolerance = sqrt(eps()))
    ϕfunc(r, origin) = let bob = BurgersVector(dis) 
        b2g_shear(r - origin, bob)
    end
    return PBCField(ϕfunc, position, DislocationOrigin(dis), Lx, Ly, tolerance)
end

PBC_cell_shift( nx, ny, Lx, Ly ) = Vector2D( nx * Lx, ny * Ly )

function _left_edge( ϕfunc::Function, square_index, position::Vector2D{T}, dis::Dislocation, Lx, Ly ) where T <: AbstractFloat
    edge_total::T = zero(T)
    left_xdx::Int = -square_index
    for ydx ∈ Iterators.reverse( UnitRange(-square_index, square_index) )
        edge_total += ϕfunc( position, shift_origin( dis, PBC_cell_shift(left_xdx, ydx, Lx, Ly) ) )
    end
    return edge_total
end

function _bottom_edge( ϕfunc::Function, square_index, position::Vector2D{T}, dis::Dislocation, Lx, Ly ) where T <: AbstractFloat
    edge_total::T = zero(T)
    bottom_ydx::Int = -square_index
    for xdx ∈ UnitRange(-square_index + one(Int), square_index)
        edge_total += ϕfunc( position, shift_origin( dis, PBC_cell_shift(xdx, bottom_ydx, Lx, Ly) ) )
    end
    return edge_total
end

function _right_edge( ϕfunc::Function, square_index, position::Vector2D{T}, dis::Dislocation, Lx, Ly ) where T <: AbstractFloat
    edge_total::T = zero(T)
    right_xdx::Int = square_index
    for ydx ∈ UnitRange(-square_index + one(Int), square_index)
        edge_total += ϕfunc( position, shift_origin( dis, PBC_cell_shift(right_xdx, ydx, Lx, Ly) ) )
    end
    return edge_total
end

function _top_edge( ϕfunc::Function, square_index, position::Vector2D{T}, dis::Dislocation, Lx, Ly ) where T <: AbstractFloat
    edge_total::T = zero(T)
    top_ydx::Int = square_index
    for xdx ∈ Iterators.reverse( UnitRange(-square_index + one(Int), square_index - one(Int)) )
        edge_total += ϕfunc( position, shift_origin( dis, PBC_cell_shift(xdx, top_ydx, Lx, Ly) ) )
    end
    return edge_total
end

function _cycle_square( ϕfunc::Function, square_index, position::Vector2D{T}, dis::Dislocation, Lx, Ly ) where T <: AbstractFloat
    square_total::T = zero(T)
    square_total += _left_edge(ϕfunc, square_index, position, dis, Lx, Ly)
    square_total += _bottom_edge(ϕfunc, square_index, position, dis, Lx, Ly)
    square_total += _right_edge(ϕfunc, square_index, position, dis, Lx, Ly)
    square_total += _top_edge(ϕfunc, square_index, position, dis, Lx, Ly)
    return square_total
end

function b2g_PBC_manual( ϕfunc::Function, position::Vector2D{T}, dis::Dislocation, Lx, Ly, tolerance::T = sqrt(eps()) ) where T <: AbstractFloat
    # First evaluate the field within its own cell
    # ϕfunc(r, d) = b2g_shear(r, d)
    output::T = ϕfunc(position, dis)
    
    # Now start computing along square trajectories
    square_index = zero(Int)
    not_complete::Bool = true
    while not_complete
        square_index += one(square_index)
        square_total::T = _cycle_square(ϕfunc, square_index, position, dis, Lx, Ly)   

        # Now check that the contribution along the square is smaller
        # than the output * tolerance (in absolute scale)
        # @show square_index
        # @show square_total
        output += square_total
        # @show output
        complete = abs(square_total) / abs(output) < tolerance

        not_complete = !complete
    end

    return output, square_index
end
    
end # module ShearFunctions