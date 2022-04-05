export vector_diff!, subtract_PBC, subtract_PBC!

"""
    vector_diff!( A::PhysicalVector, B::PhysicalVector; diff!::Function )

Wrapper around the in-place difference function `diff!`. 

# Additional Information
To be used when geometry, etc., may change the definition of the `PhysicalVector` `subtract!`.
"""
vector_diff!( a::PhysicalVector, b::PhysicalVector; diff::Function = subtract! ) = diff(a, b)

"""
    subtract_PBC( x1::Real, x0::Real; axis_size )

Subtraction defined on a _ring_ with periodic boundary conditions of size `axis_size`.

# Examples
```jldocstest
julia> subtract_PBC( 2., 1.; axis_size = 10 )
1.0

julia> subtract_PBC( 7., 1.; axis_size = 10 )
-4.0

julia> subtract_PBC( 1., 8.; axis_size = 10 )
3.0
```
"""
function subtract_PBC( x1::Real, x0::Real; axis_size )
    value = x1 - x0 
    while value > axis_size/2
        value -= axis_size
    end
    while value < -axis_size/2
        value += axis_size
    end
    return value
end

"""
    subtract_PBC!( A::Vector2D, B::Vector2D; Lx, Ly = Lx )

In-place `Vector2D` subtraction with periodic boundary conditions.

# Examples
```jldocstest
julia> A = Vector2D(1., 2.)
Vector2D{Float64}([1.0, 2.0])

julia> B = Vector2D(7., 9.)
Vector2D{Float64}([7.0, 9.0])

julia> subtract_PBC!(A, B; Lx = 10)

julia> A
Vector2D{Float64}([4.0, 3.0])
```
"""
function subtract_PBC!( A::Vector2D, B::Vector2D; Lx, Ly = Lx )
    A.vec[1] = subtract_PBC( A.vec[1], B.vec[1]; axis_size = Lx )
    A.vec[2] = subtract_PBC( A.vec[2], B.vec[2]; axis_size = Ly )
    return nothing
end