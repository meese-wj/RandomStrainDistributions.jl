module PhysicalVectors

# Exports from PhysicalVector interface
export PhysicalVector, equal, add!, multiply!, ⋅, subtract!, divide!, dot, magnitude, magnitude2, normalize, normalize!, unit, unit!

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
equal(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Add two PhysicalVectors together.
"""
Base.:+(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Add two PhysicalVectors together in-place. (Here A is modified, not B).
"""
add!(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Scalar multiplication of a PhysicalVector.
"""
Base.:*(λ, A::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

In-place scalar multiplication of a PhysicalVector.
"""
multiply!(λ, A::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")
"""
#= REQUIRED FOR PhysicalVector INTERFACE =#

Scalar (dot) product of two PhysicalVectors.
"""
⋅(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

#= ========================================================================================= =#
#  These are functions that follow the interface but are not
#  necessary to define for subtypes of PhysicalVector
#= ========================================================================================= =#
"""
Scalar multiplication of a PhysicalVector. Alias for Base.:*(λ, ::PhysicalVector)
"""
Base.:*(A::PhysicalVector, λ) = λ * A
multiply!(A::PhysicalVector, λ) = multiply!(λ, A) 
Base.:-(A::PhysicalVector, B::PhysicalVector) = A + (-1 * B)
Base.:/(A::PhysicalVector, λ) = A * (1 / λ)
divide!(A::PhysicalVector, λ) = multiply!(A, 1 / λ)
dot(A::PhysicalVector, B::PhysicalVector) = A ⋅ B 
magnitude2(A::PhysicalVector) = A ⋅ A 
magnitude(A::PhysicalVector) = sqrt(magnitude2(A))
normalize(A::PhysicalVector) = A / magnitude(A) 
normalize!(A::PhysicalVector) = divide!(A, magnitude(A)) 
unit(A::PhysicalVector) = normalize(A) 
unit!(A::PhysicalVector) = normalize!(A) 

#= ========================================================================================= =#
#  Include implementations of the PhysicalVector interface.
#  Each implementation file is included and contains it's own export list.
#= ========================================================================================= =#
include("Vectors2D.jl")

end # module PhysicalVectors