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
"""
function bxg_shears!( shears::MVector{2}, eval_r::Vector2D, bob::Vector2D; diff::Function = subtract!, source_r::Vector2D = Vector2D(0.,0.) )
    vector_diff!(eval_r, source_r; diff = diff)
    shears[1] = b1g_shear( eval_r, bob )
    shears[2] = b2g_shear( eval_r, bob )
    return nothing
end

end # module RandomStrainDistributions
