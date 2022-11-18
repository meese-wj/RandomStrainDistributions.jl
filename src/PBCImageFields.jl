"""
    module PBCImageFields

A submodule to calculate values of a long-range
field on a periodic square lattice. 
"""
module PBCImageFields

using ..PhysicalVectors
using ..CrystalDefects

export PBCField, image_origin, _left_edge

# """
#     image_origin(origin, nx, ny, Lx, Ly)

# Compute the origin for the image dislocation in the (nx, ny) image cell.
# """
# image_origin( origin, nx, ny, Lx, Ly, xhat = Vector2D(1.0, 0.0), yhat = Vector2D(0.0, 1.0) ) = origin + nx * Lx * xhat + ny * Ly * yhat

# image_contribution( ϕfunc::Function, position::Vector2D{T}, _image_origin ) where T = ϕfunc( position, _image_origin )

# function _left_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origin, Lx, Ly ) where T
#     edge_total::T = zero(T) 
#     left_xdx::Int = -square_index
#     for ydx ∈ Iterators.reverse( UnitRange(-square_index, square_index) )
#         edge_total += image_contribution( ϕfunc, position, image_origin(origin, left_xdx, ydx, Lx, Ly) )
#     end
#     return edge_total
# end

# function _bottom_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origin, Lx, Ly ) where T
#     edge_total::T = zero(T) 
#     bottom_ydx::Int = -square_index
#     @inbounds for xdx ∈ UnitRange(-square_index + one(square_index), square_index)
#         edge_total += image_contribution( ϕfunc, position, image_origin(origin, xdx, bottom_ydx, Lx, Ly) )
#     end
#     return edge_total
# end

# function _right_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origin, Lx, Ly ) where T
#     edge_total::T = zero(T) 
#     right_xdx::Int = square_index
#     @inbounds for ydx ∈ UnitRange(-square_index + one(square_index), square_index)
#         edge_total += image_contribution(ϕfunc, position, image_origin(origin, right_xdx, ydx, Lx, Ly))
#     end
#     return edge_total
# end

# function _top_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origin, Lx, Ly ) where T
#     edge_total::T = zero(T) 
#     top_ydx::Int = square_index
#     @inbounds for xdx ∈ Iterators.reverse( UnitRange(-square_index + one(square_index), square_index - one(square_index)) )
#         edge_total += image_contribution(ϕfunc, position, image_origin(origin, xdx, top_ydx, Lx, Ly))
#     end
#     return edge_total
# end

# function _cycle_square(ϕfunc::Function, square_index, position::Vector2D{T}, origin, Lx, Ly) where T
#     local square_total::T = zero(T)
#     # Start in top left and move down at constant x 
#     square_total += _left_edge(ϕfunc, square_index, position, origin, Lx, Ly)

#     # Now move along the bottom edge at constant y
#     square_total += _bottom_edge(ϕfunc, square_index, position, origin, Lx, Ly)
    
#     # Next move along the right edge at constant x
#     square_total += _right_edge(ϕfunc, square_index, position, origin, Lx, Ly)
    
#     # Finally move along the top edge at constant y
#     square_total += _top_edge(ϕfunc, square_index, position, origin, Lx, Ly)
#     return square_total
# end

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

"""
    PBCField(ϕfunc::Function, position, origin, Lx, Ly, [tolerance])

Calculate the `ϕfunc`tion on a square periodic lattice of
size Lx, Ly. 

# Additional Information

* By default, the `tolerance` variable has the value of `sqrt( eps() )`.

!!! note
    The image contribution is summed along square trajectories.
    If the `ϕfunc` field has an angular dependence, then this 
    square trajectory should converge. If one supplies a divergent
    field, however, then this sum will always clearly diverge.
"""
function PBCField( ϕfunc::Function, position::Vector2D{T}, dis::Dislocation, Lx, Ly, tolerance::T = sqrt(eps()) ) where T <: AbstractFloat
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
# function PBCField(ϕfunc::Function, position::Vector2D{T}, origin, Lx, Ly, tolerance::T = sqrt(eps()) ) where T <: AbstractFloat
#     # First evaluate the field within its own cell
#     output::T = ϕfunc(position, origin)
    
#     # Now start computing along square trajectories
#     square_index = zero(Int)
#     not_complete::Bool = true
#     while not_complete
#         square_index += one(square_index)
#         square_total::T = _cycle_square(ϕfunc, square_index, position, origin, Lx, Ly)

#         # Now check that the contribution along the square is smaller
#         # than the output * tolerance (in absolute scale)
#         output += square_total
#         complete = abs(square_total) / abs(output) < tolerance

#         not_complete = !complete
#     end

#     return output, square_index
# end


    
end # module PBCImageFields