module RandomStrainDistributions
using Reexport

include("PhysicalVectors/PhysicalVectors.jl")
@reexport using .PhysicalVectors

# ===================================================================== #
#  RandomStrainDistributions module code is below
# ===================================================================== #
using StaticArrays

export b1g_shear, b2g_shear, bxg_shears!

"""
    b1g_shear( eval_r::Vector2D, bob::Vector2D )

Calculate the ``B_{1g}`` from an edge dislocation for a given Burgers vector `bob` at site `eval_r`. 

# Additional Information
* The edge dislocation is aligned along the ``z`` axis at the origin.
* This is in units of ``1/(1-σ)`` for a Poisson ratio ``σ``.

# Examples
```jldocstest
julia> b1g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
-0.07957747154594767
```
"""
function b1g_shear( eval_r::Vector2D, bob::Vector2D )
    return -1/π * ( eval_r.vec[1] * eval_r.vec[2] / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end

"""
    b2g_shear( eval_r::Vector2D, bob::Vector2D )

Calculate the ``B_{2g}`` from an edge dislocation for a given Burgers vector `bob` at site `eval_r`. 

# Additional Information
* The edge dislocation is aligned along the ``z`` axis at the origin.
* This is in units of ``1/(1-σ)`` for a Poisson ratio ``σ``.

# Examples
```jldocstest
julia> b2g_shear( Vector2D(1.,1.), Vector2D(1.,0.) )
0.15915494309189535
```
"""
function b2g_shear( bob::Vector2D, eval_r::Vector2D )
    return 1/(2*π) * ( ( eval_r.vec[1]^2 - eval_r.vec[2]^2 ) / magnitude2(eval_r) ) * (bob ⋅ eval_r) / magnitude2(eval_r)
end

"""
    bxg_shears!(shears::MVector{2}, eval_r::Vector2D, bob::Vector2D; diff::Function, source_r::Vector2D = Vector2D(0.,0.)) -> nothing 

Calculate the two `shears` from the ``B_{1g}`` and ``B_{2g}`` channels from and edge dislocation situated at the `source_r`
location and aligned along the ``z`` axis with a Burger's vector `bob`.

# Additional Information
* This function mutates both `shears` and `eval_r`. It should be used in cases where `eval_r` is a temporary vector.

# Examples
```jldocstest
julia> using StaticArrays

julia> A = MVector{2}(0., 0.)
2-element MVector{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> bxg_shears!( A, Vector2D(1.,1.), Vector2D(1.,0.) )

julia> A
2-element MVector{2, Float64} with indices SOneTo(2):
 -0.07957747154594767
  0.15915494309189535
```
"""
function bxg_shears!( shears::MVector{2}, eval_r::Vector2D, bob::Vector2D; diff::Function = Base.:-, source_r::Vector2D = Vector2D(0.,0.) )
    vector_diff!(eval_r, source_r; diff = diff)
    broadcast!( (x, v) -> v, shears, shears, ( b1g_shear(eval_r, bob), b2g_shear(eval_r, bob) ) )
    return nothing
end

end # module RandomStrainDistributions
