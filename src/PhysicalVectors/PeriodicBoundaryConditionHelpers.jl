export subtract_PBC, subtract_PBC!

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
    broadcast!( (a, b, ax) -> subtract_PBC(a, b; axis_size=ax), A.vec, A.vec, B.vec, (Lx, Ly) )
    return nothing
end