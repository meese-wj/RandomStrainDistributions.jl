"""
    module PBCImageFields

A submodule to calculate values of a long-range
field on a periodic square lattice. 
"""
module PBCImageFields

using ..PhysicalVectors

export PBCField, image_origin

"""
    image_origin(origin, nx, ny, Lx, Ly)

Compute the origin for the image dislocation in the (nx, ny) image cell.
"""
image_origin( origin, nx, ny, Lx, Ly, xhat = Vector2D(1.0, 0.0), yhat = Vector2D(0.0, 1.0) ) = origin + nx * Lx * xhat + ny * Ly * yhat

image_contribution( ϕfunc::Function, position::Vector2D{T}, _image_origin ) where T = ( ϕfunc( position, _image_origin ) )::T

function _left_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origins, Lx, Ly ) where T
    edge_total::T = zero(T) 
    left_xdx::Int = -square_index
    @inbounds for ydx ∈ UnitRange(square_index, -square_index)
        edge_total += image_contribution( ϕfunc, position, image_origin.(origins, left_xdx, ydx, Lx, Ly) )
    end
    return edge_total
end

function _bottom_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origins, Lx, Ly ) where T
    edge_total::T = zero(T) 
    bottom_ydx = -square_index
    @inbounds for xdx ∈ UnitRange(-square_index, square_index)
        edge_total += image_contribution( ϕfunc, position, image_origin.(origins, xdx, bottom_ydx, Lx, Ly) )
    end
    return edge_total
end

function _right_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origins, Lx, Ly ) where T
    edge_total::T = zero(T) 
    right_xdx = square_index
    @inbounds for ydx ∈ UnitRange(-square_index, square_index)
        edge_total += image_contribution(ϕfunc, position, image_origin.(origins, right_xdx, ydx, Lx, Ly))
    end
    return edge_total
end

function _top_edge( ϕfunc::Function, square_index, position::Vector2D{T}, origins, Lx, Ly ) where T
    edge_total::T = zero(T) 
    top_ydx = square_index
    @inbounds for xdx ∈ UnitRange(square_index, -square_index)
        edge_total += image_contribution(ϕfunc, position, image_origin.(origins, xdx, top_ydx, Lx, Ly))
    end
    return edge_total
end

"""
    PBCField(ϕfunc::Function, position, origins, Lx, Ly, [tolerance])

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
function PBCField(ϕfunc::Function, position::Vector2D{T}, origins, Lx, Ly, tolerance::T = sqrt(eps()) ) where T <: AbstractFloat
    # First evaluate the field within its own cell
    output::T = ϕfunc(position, origins)
    
    # Now start computing along square trajectories
    square_index = zero(Int)
    not_complete::Bool = true
    while not_complete
        square_index::Int += one(square_index)
        square_total::T = zero(output)

        # Start in top left and move down at constant x 
        square_total += _left_edge(ϕfunc, square_index, position, origins, Lx, Ly)

        # Now move along the bottom edge at constant y
        square_total += _bottom_edge(ϕfunc, square_index, position, origins, Lx, Ly)
        
        # Next move along the right edge at constant x
        square_total += _right_edge(ϕfunc, square_index, position, origins, Lx, Ly)
        
        # Finally move along the top edge at constant y
        square_total += _top_edge(ϕfunc, square_index, position, origins, Lx, Ly)

        # Now check that the contribution along the square is smaller
        # than the output * tolerance (in absolute scale)
        output += square_total
        complete = abs(square_total) / abs(output) < tolerance

        not_complete = !complete
    end

    return output, square_index
end


    
end # module PBCImageFields