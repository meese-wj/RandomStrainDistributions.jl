"""
Module to house an interface of vectors commonly used in physics.
"""
module PhysicalVectors

import Base

# Exports from PhysicalVector interface
export PhysicalVector, isequal, add!, multiply!, ⋅, subtract!, divide!, dot, magnitude, magnitude2, normalize, normalize!, unit, unit!

"""
This is an interface for vectors used in physics.

# Additional Information 

The only functions that the user must define for a
`VectorType <: PhysicalVector` are:

* Equality:                       `Base.isequal(vector1, vector2)`
* Addition:                       `Base.:+(vector1, vector2)`
* In-place Addition:              `add!(vector1, vector2)`
* In-place Subtraction:           `subract!(vector1, vector2)`
* Scalar Multiplication:          `Base.:*(scalar, vector)`
* In-place Scalar Multiplication: `multiply!(vector, scalar)`
* Scalar Product:                 `⋅(vector1, vector2) #\\cdot [tab]`

The other exported functions are standard vector 
operations which can be written in terms of those
above:

* Scalar Multiplication:          `Base.:*(vector, scalar)`
* Subtraction:                    `Base.:-(vector1, vector2)`
* Divison by Scalar:              `Base.:/(vector, scalar)`
* In-place Divison by Scalar:     `divide!(vector, scalar)`
* Scalar Product:                 `dot(vector1, vector2)`
* Square Magnitude:               `magnitude2(vector)`
* Magnitude:                      `magnitude(vector)`
* Normalization:                  `normalize(vector)`
* In-place Normalization:         `normalize!(vector)`
* Unit:                           `unit(vector)`
* In-place Unit:                  `unit!(vector)`
"""
abstract type PhysicalVector end

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    Base.isequal(A::PhysicalVector, B::PhysicalVector)

Check for equality between two `PhysicalVector`s. 
    
# Additional Information 
* Note that Base.:(==) is not used for issues with mutability and equality versus egality.
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
Base.isequal(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    +(A::PhysicalVector, B::PhysicalVector)

Addtion of two `PhysicalVector`s. 
    
# Additional Information 
* A call to a copy constructor is implied.
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
Base.:+(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    add!(A::PhysicalVector, B::PhysicalVector)

In-place addition of two `PhysicalVector`s. 
    
# Additional Information 
* It is implied that the first argument, `A`, is modified in-place.
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
add!(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    *(λ, A::PhysicalVector)

Multiplication of a `PhysicalVector` by the scalar `λ`. 
    
# Additional Information 
* A call to a copy constructor is implied.
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
Base.:*(λ, A::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    multiply!(A::PhysicalVector, λ)

In-place multiplication of a `PhysicalVector` by the scalar `λ`. 
    
# Additional Information 
* It is implied that the first argument, `A`, is modified in-place.
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
multiply!(A::PhysicalVector, λ) = error("No implementation defined for vectors of type $(typeof(λ)) and $(typeof(A)).")

"""
# REQUIRED FOR `PhysicalVector` INTERFACE 

    ⋅(A::PhysicalVector, B::PhysicalVector) #\\cdot [tab]

Scalar product of two `PhysicalVector`s. This is also referred to as the `dot` product. 
    
# Additional Information 
* This function `error`s out in the default implementation to enforce its definition in subtypes of `PhysicalVector`.
"""
⋅(A::PhysicalVector, B::PhysicalVector) = error("No implementation defined for vectors of type $(typeof(A)) and $(typeof(B)).")

#= ========================================================================================= =#
#  These are functions that follow the interface but are not
#  necessary to define for subtypes of PhysicalVector
#= ========================================================================================= =#
"""
    *(A::PhysicalVector, λ)

Multiplication of a `PhysicalVector` by the scalar `λ`. 
    
# Additional Information 
* This is an alias for `*(λ, A::PhysicalVector)`.
"""
Base.:*(A::PhysicalVector, λ) = λ * A

"""
    -(A::PhysicalVector, B::PhysicalVector)

Subtraction of two `PhysicalVector`s. 
    
# Additional Information 
* A copy constructor call is implied.
"""
Base.:-(A::PhysicalVector, B::PhysicalVector) = A + (-1 * B)

"""
    /(A::PhysicalVector, λ)

Division of a `PhysicalVector` by the scalar `λ`. 
    
# Additional Information 
* A copy constructor call is implied.
"""
Base.:/(A::PhysicalVector, λ) = A * (1 / λ)

"""
    divide!(A::PhysicalVector, λ)

Division of a `PhysicalVector` by the scalar `λ`. 
    
# Additional Information 
* It is implied that the first argument, `A`, is modified in-place.
"""
divide!(A::PhysicalVector, λ) = multiply!(A, 1 / λ)

"""
    dot(A::PhysicalVector, B::PhysicalVector)

Scalar product of two `PhysicalVector`s. 
    
# Additional Information 
* This is an alias call to `⋅(A::PhysicalVector, B::PhysicalVector) #\\cdot [tab]`
"""
dot(A::PhysicalVector, B::PhysicalVector) = A ⋅ B 

"""
    magnitude2(A::PhysicalVector)

The square magnitude of a `PhysicalVector`. 
    
# Additional Information 
* This is an alias call to `⋅(A::PhysicalVector, A::PhysicalVector) #\\cdot [tab]`
"""
magnitude2(A::PhysicalVector) = A ⋅ A 

"""
    magnitude(A::PhysicalVector)

The magnitude of a `PhysicalVector`. 
    
# Additional Information 
* This is an alias call to `sqrt( magnitude2(A::PhysicalVector) )`
"""
magnitude(A::PhysicalVector) = sqrt(magnitude2(A))

"""
    normalize(A::PhysicalVector)

Divide a `PhysicalVector` by its length. 
    
# Additional Information 
* A call to a copy constructor is implied. 
"""
normalize(A::PhysicalVector) = A / magnitude(A) 

"""
    normalize!(A::PhysicalVector)

Divide a `PhysicalVector` by its length. 
    
# Additional Information 
* It is implied that the argument, `A`, is modified in-place.
"""
normalize!(A::PhysicalVector) = divide!(A, magnitude(A)) 

"""
    unit(A::PhysicalVector)

Divide a `PhysicalVector` by its length.
    
# Additional Information 
* This is an alias for `normalize(A::PhysicalVector)`.
* A call to a copy constructor is implied. 
"""
unit(A::PhysicalVector) = normalize(A) 

"""
    unit!(A::PhysicalVector)

Divide a `PhysicalVector` by its length. 
    
# Additional Information 
* This is an alias for `normalize!(A::PhysicalVector)`.
* It is implied that the argument, `A`, is modified in-place.
"""
unit!(A::PhysicalVector) = normalize!(A) 

#= ========================================================================================= =#
#  Include implementations of the PhysicalVector interface.
#  Each implementation file is included and contains it's own export list.
#= ========================================================================================= =#
include("Vectors2D.jl")

end # module PhysicalVectors