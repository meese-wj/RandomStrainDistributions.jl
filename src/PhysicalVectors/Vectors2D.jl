import Base

# Exports from Vector2D implementation
export Vector2D

"""
    struct Vector2D{T <: Real} <: PhysicalVector

Wrapper around `NTuple{2, T}` for convenience in physics applications.
"""
struct Vector2D{T <: Real} <: PhysicalVector
    vec::NTuple{2, T}
end

Base.show(io::IO, vec::Vector2D{T}) where T = Base.print(io, "Vector2D{$T}($(vec.vec[1]), $(vec.vec[2]))")
Base.show(vec::Vector2D) = Base.show(stdout, vec)

#= ========================================================================================= =#
#  Struct constructors
#= ========================================================================================= =#
"""
    Vector2D(arr::AbstractArray{Real, 1})

Copy constructor from an `AbstractVector`

# Examples

```jldoctest
julia> Vector2D([1., 2])
Vector2D{Float64}(1.0, 2.0)
```
"""
Vector2D(arr::AbstractArray{T,1}) where {T <: Real} = Vector2D{T}( (arr[1], arr[2]) )

"""
    Vector2D(vec::Vector2D)

Copy constructor from an `Vector2D`

# Examples

```jldoctest
julia> Vector2D( Vector2D([1., 2]) )
Vector2D{Float64}(1.0, 2.0)
```
"""
Vector2D(vec::Vector2D{T}) where {T <: Real} = Vector2D( vec.vec )

"""
    Vector2D(x::Real, y::Real)

Copy constructor from a pair x and y of potentially different `Real` types.

# Additional Information
This function makes a call to `Vector2D(::Tuple{S, T})` where `{S <: Real, T <: Real}`.

# Examples

```jldoctest
julia> Vector2D(1, 2.)
Vector2D{Float64}(1.0, 2.0)
```
"""
Vector2D(x::Real, y::Real) = Vector2D((x, y))

"""
    Vector2D(tup::Tuple{S, T}) where {S <: Real, T <: Real}

Copy constructor from a tuple of reals `tup === (x::Real, y::Real)`. 

# Additional Information
In the case that the elements are of different types, they are first
passed through `promote(x, y)` and the common-type values are stored
in the output.

# Examples

```jldoctest
julia> Vector2D((1, 2.))
Vector2D{Float64}(1.0, 2.0)
```
"""
function Vector2D(tup::Tuple{S, T}) where {S <: Real, T <: Real}
    vals = promote(tup...)
    return Vector2D{eltype(vals)}( vals )
end

Vector2D{T}(::UndefInitializer) where T = Vector2D( T(undef), T(undef) )   

#= ========================================================================================= =#
#  Interface requirement definitions.
#= ========================================================================================= =#
"""
    Base.isequal(A::Vector2D, B::Vector2D)

Vector2D _equality_ based on component-wise equality.

# Additional Information
Use this function for equality instead of `==`.

# Examples

```jldoctest
julia> A = Vector2D(1., 2)
Vector2D{Float64}(1.0, 2.0)

julia> B = Vector2D(1, 2.)
Vector2D{Float64}(1.0, 2.0)

julia> isequal(A, B)
true
```
"""
Base.isequal(A::Vector2D{S}, B::Vector2D{T}) where {S, T} = prod(A.vec .== B.vec)

"""
    +(A::Vector2D, B::Vector2D) -> Vector2D

Addition for `Vector2D`. Includes a copy constructor.

```jldoctest
julia> A = Vector2D(1., 2)
Vector2D{Float64}(1.0, 2.0)

julia> B = Vector2D(1, 2.)
Vector2D{Float64}(1.0, 2.0)

julia> A + B
Vector2D{Float64}(2.0, 4.0)
```
"""
Base.:+(A::Vector2D, B::Vector2D) = Vector2D( broadcast(+, A.vec, B.vec) )

"""
    *(λ, A::Vector2D) -> Vector2D

Multiplication by a scalar `λ` for `Vector2D`. Includes a copy constructor.

```jldoctest
julia> A = Vector2D(1., 2)
Vector2D{Float64}(1.0, 2.0)

julia> λ = 2.
2.0

julia> λ * A
Vector2D{Float64}(2.0, 4.0)
```
"""
Base.:*(λ, A::Vector2D) = Vector2D( broadcast(*, A.vec, λ) )

"""
    ⋅(A::Vector2D, B::Vector2D)

Scalar product, i.e. dot product, for `Vector2D`. 

```jldoctest
julia> A = Vector2D(1., 2)
Vector2D{Float64}([1.0, 2.0])

julia> B = Vector2D(1, 2.)
Vector2D{Float64}([1.0, 2.0])

julia> A ⋅ B
5.0
```
"""
function ⋅(A::Vector2D, B::Vector2D)
    ans = zero(eltype(A.vec))
    for idx in eachindex(A.vec, B.vec)
        ans += A.vec[idx] * B.vec[idx]
    end
    return ans
end

#= ========================================================================================= =#
#  Includes for the Vector2D struct.
#  Each file **must** contain it's own export list.
#= ========================================================================================= =#
include("PeriodicBoundaryConditionHelpers.jl")
