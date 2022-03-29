module PhysicalVectors

# Exports from PhysicalVector interface
export PhysicalVector, equal, add!, multiply!, ⋅, subtract!, divide!, dot, magnitude, magnitude2, normalize, normalize!, unit, unit!
# Exports from Vector2D implementation
export Vector2D, equal, add!, multiply!

"""
This is an interface for vectors used in physics.

The only functions that the user must define for a
vector <: PhysicalVector are:

    -- Equality:                       equal(vector1, vector2)
    -- Addition:                       Base.:+(vector1, vector 2)
    -- In-place Addition:              add!(vector1, vector 2)
    -- In-place Subtraction:           subract!(vector1, vector 2)
    -- Scalar Multiplication:          Base.:*(scalar, vector)
    -- In-place Scalar Multiplication: multiply!(scalar, vector)
    -- Scalar Product:                 ⋅(vector1, vector2) #\\cdot [tab]

The other exported functions are standard vector 
operations which can be written in terms of those
above:

    -- Scalar Multiplication:          Base.:*(vector, scalar)
    -- In-place Scalar Multiplication: multiply!(vector, scalar)
    -- Subtraction:                    Base.:-(vector1, vector2)
    -- Divison by Scalar:              Base.:/(vector, scalar)
    -- In-place Divison by Scalar:     divide!(vector, scalar)
    -- Scalar Product:                 dot(vector1, vector2)
    -- Square Magnitude:               magnitude2(vector)
    -- Magnitude:                      magnitude(vector)
    -- Normalization:                  normalize(vector)
    -- In-place Normalization:         normalize!(vector)
    -- Unit:                           unit(vector)
    -- In-place Unit:                  unit!(vector)
"""
abstract type PhysicalVector end

"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Check for equality between two PhysicalVectors. Note that Base.:(==) is not used 
for issues with immutability and equality.
"""
equal(A::T, B::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Add two PhysicalVectors together.
"""
Base.:+(A::T, B::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Add two PhysicalVectors together in-place. (Here A is modified, not B).
"""
add!(A::T, B::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Scalar multiplication of a PhysicalVector.
"""
Base.:*(λ, A::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

In-place scalar multiplication of a PhysicalVector.
"""
multiply!(λ, A::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Scalar (dot) product of two PhysicalVectors.
"""
⋅(A::T, B::T) where {T <: PhysicalVector} = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

# These are functions that follow the interface but are not
# necessary to define for subtypes of PhysicalVector
Base.:*(A::T, λ) where {T <: PhysicalVector} = λ * A
multiply!(A::T, λ) where {T <: PhysicalVector} = multiply!(λ, A) 
Base.:-(A::T, B::T) where {T <: PhysicalVector} = A + (-1 * B)
Base.:/(A::T, λ) where {T <: PhysicalVector} = A * (1 / λ)
divide!(A::T, λ) where {T <: PhysicalVector} = multiply!(A, 1 / λ)
dot(A::T, B::T) where {T <: PhysicalVector} = A ⋅ B 
magnitude2(A::T) where {T <: PhysicalVector} = A ⋅ A 
magnitude(A::T) where {T <: PhysicalVector} = sqrt(magnitude2(A))
normalize(A::T) where {T <: PhysicalVector} = A / magnitude(A) 
normalize!(A::T) where {T <: PhysicalVector} = divide!(A, magnitude(A)) 
unit(A::T) where {T <: PhysicalVector} = normalize(A) 
unit!(A::T) where {T <: PhysicalVector} = normalize!(A) 

using StaticArrays

"""
Wrapper around StaticArrays.MVector{2,T} for simpler use in 
physics applications.
"""
mutable struct Vector2D{T <: Real} <: PhysicalVector
    vec::MVector{2,T}
end

Vector2D(x::T, y::T) where {T <: Real} = Vector2D{T}( MVector{2,T}(x,y) )

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

end # module PhysicalVectors