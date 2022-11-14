"""
    module PBCImageFields

A submodule to calculate values of a long-range
field on a periodic square lattice. 
"""
module PBCImageFields

using ..PhysicalVectors

export PBCField

"""
    PBCField(ϕfunc::Function, position, Lx, Ly, [tolerance])

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
function PBCField(ϕfunc::Function, position, origin::Vector2D, Lx, Ly, tolerance = sqrt(eps()) )
    # First evaluate the field within its own grid
    output = ϕfunc(position, origin)
    
    # Now start computing along square trajectories
    old_output = output
    square_index = zero(Int)
    not_complete = true
    while ( not_complete || square_index == zero(Int) )
        square_index += one(square_index)
        
        # max_abs_out = zero(output)
        square_total = zero(output)

        @show square_index, output
        # Start in top left and move down at constant x 
        left_xdx = -square_index
        for ydx ∈ range(-square_index, square_index, step = 1)
            image_origin = origin + left_xdx * Lx * Vector2D( 1.0, 0.0 ) + ydx * Ly * Vector2D( 0.0, 1.0 )
            new_out = ϕfunc(position, image_origin)
            square_total += new_out
            # max_abs_out = abs(new_out) > max_abs_out ? abs(new_out) : max_abs_out
        end

        # Now move along the bottom edge at constant y
        bottom_ydx = -square_index
        for xdx ∈ range(-square_index, square_index, step = 1)
            image_origin = origin + xdx * Lx * Vector2D( 1.0, 0.0 ) + bottom_ydx * Ly * Vector2D( 0.0, 1.0 )
            new_out = ϕfunc(position, image_origin)
            square_total += new_out
            # max_abs_out = abs(new_out) > max_abs_out ? abs(new_out) : max_abs_out
        end
        
        # Next move along the right edge at constant x
        right_xdx = square_index
        for ydx ∈ range(square_index, -square_index, step = 1)
            image_origin = origin + right_xdx * Lx * Vector2D( 1.0, 0.0 ) + ydx * Ly * Vector2D( 0.0, 1.0 )
            new_out = ϕfunc(position, image_origin)
            square_total += new_out
            # max_abs_out = abs(new_out) > max_abs_out ? abs(new_out) : max_abs_out
        end
        
        # Finally move along the top edge at constant y
        top_ydx = square_index
        for xdx ∈ range(square_index, -square_index, step = 1)
            image_origin = origin + xdx * Lx * Vector2D( 1.0, 0.0 ) + top_ydx * Ly * Vector2D( 0.0, 1.0 )
            new_out = ϕfunc(position, image_origin)
            square_total += new_out
            # max_abs_out = abs(new_out) > max_abs_out ? abs(new_out) : max_abs_out
        end

        # Now check to see if the cutoff (to within a tolerance)
        # was reached. This is done by multiplying the perimeter by
        # the maximum bound. If this is smaller than the tolerance * output
        # then the next iteration will definitely be convergent and we end
        # the summation.
        # @show perimeter = 4 * ( 2 * square_index + 1 )
        # @show max_bound = perimeter * max_abs_out
        # @show complete = max_bound / abs(output) < tolerance
        @show square_total
        @show output += square_total
        @show complete = abs(square_total) / abs(output) < tolerance

        not_complete = !complete
    end

    output, square_index
end


    
end # module PBCImageFields