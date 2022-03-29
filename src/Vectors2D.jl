using StaticArrays

# Exports from Vector2D implementation
export Vector2D, equal, add!, multiply!

"""
Wrapper around StaticArrays.MVector{2,T} for simpler use in 
physics applications.
"""
mutable struct Vector2D{T <: Real} <: PhysicalVector
    vec::MVector{2,T}
end

#= ========================================================================================= =#
#  Struct constructors
#= ========================================================================================= =#
"""
Vector2D constructor from an AbstractVector
"""
Vector2D(arr::AbstractArray{T,1}) where {T <: Real} = Vector2D{T}( MVector{2,T}(arr[1], arr[2]) )
"""
Vector2D constructor from a pair x and y of the same type T 
"""
Vector2D(x::T, y::T) where {T <: Real} = Vector2D{T}( MVector{2,T}(x,y) )
"""
Vector2D constructor from a pair x and y of potentially different types.
This function makes a call to Vector2D(::Tuple{S, T}) where {S <: Real, T <: Real}.
"""
Vector2D(x::Real, y::Real) = Vector2D((x, y))
"""
Vector2D constructor from a tuple of reals (x, y). In the case that the
elements are of different types, they are first passed through promote(x, y)
and the common-type values are stored in the output.
"""
function Vector2D(tup::Tuple{S, T}) where {S <: Real, T <: Real}
    vals = promote(tup...)
    return Vector2D{eltype(vals)}( MVector{2, eltype(vals)}(vals) )
end

#= ========================================================================================= =#
#  Interface requirement definitions.
#= ========================================================================================= =#
"""
2D Vector Equality based on value equality in the vectors
"""
equal(A::Vector2D, B::Vector2D) = prod(A.vec .== B.vec)
"""
2D Vector Addition
"""
Base.:+(A::Vector2D, B::Vector2D) = Vector2D(A.vec[1] + B.vec[1], A.vec[2] + B.vec[2])
"""
2D Vector Scalar Multiplication
"""
Base.:*(λ, A::Vector2D) = Vector2D( λ * A.vec[1], λ * A.vec[2] )
"""
2D Vector Scalar Product
"""
⋅(A::Vector2D, B::Vector2D) = A.vec[1] * B.vec[1] + A.vec[2] * B.vec[2]

"""
In-place Vector2D addition. 

Examples:

    A = [1., 2.]; B = [3., 4.]
    add!(A, B)
    A = [4., 6.]
"""
add!(A::Vector2D, B::Vector2D) = A.vec .+= B.vec

"""
In-place Vector2D addition. 

Examples:

    A = [1., 2.]; B = [3., 4.]
    subract!(A, B)
    A = [-2., -2.]
"""
subtract!(A::Vector2D, B::Vector2D) = A.vec .-= B.vec

"""
In-place Vector2D scalar multiplication. 

Examples:

    A = [1., 2.]; λ = 2.
    multiply!(λ, A)
    A = [2., 4.]
"""
multiply!(λ, A::Vector2D) = A.vec .*= λ