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
    bxg_shears!(eval_r::Vector2D, bob::Vector2D; diff::Function, source_r::Vector2D = Vector2D(0.,0.)) -> Tuple{T, T} where {T <: Real}

Calculate the two shears `(b1g, b2g)` from the ``B_{1g}`` and ``B_{2g}`` channels from and edge dislocation situated at the `source_r`
location and aligned along the ``z`` axis with a Burger's vector `bob`. 

# Additional Information
* This function mutates both and `eval_r`. It should be used in cases where `eval_r` is a temporary vector.
* The keyword argument `diff <: Function` should be used in cases where `Vector2D` subtraction has a different definition than the normal one expected, for example with periodic boundary conditions.

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
function bxg_shears!( eval_r::Vector2D, bob::Vector2D; diff::Function = subtract!, source_r::Vector2D = Vector2D(0.,0.) )
    diff(eval_r, source_r)
    b1g = b1g_shear(eval_r, bob)
    b2g = b2g_shear(eval_r, bob)
    return (b1g, b2g)
end

using Distributions
export RandomDislocation, UniformBurgersVector, RandomStrainDistribution, RandomDislocationDistribution

"""
Interface type for the joint-distribution for the dislocations.

# Additional Information
For any given dislocation, it has a random Burgers vector and random location in the lattice. 
This type should define both additional distributions for the dislocations.
"""
abstract type RandomDislocation end

@enum DislocationProperties begin 
    BurgersVector
    DislocationOrigin
end

"""
    rand_burger_vector( rbv::RandomDislocation )

Interface function for the `RandomDislocation` type. 
"""
rand_burger_vector( rbv::RandomDislocation ) = error("No implementation defined for $(typeof(rbv)).")

"""
    rand_dislocation_source( rbv::RandomDislocation )

Interface function for the `RandomDislocation` type. 
"""
rand_dislocation_source( rbv::RandomDislocation ) = error("No implementation defined for $(typeof(rbv)).")

"""
    rand_dislocation( rbv::UniformBurgersVector ) -> Tuple{Vector2D, Vector2D}

Generate a random dislocation by defining its Burgers vector and origin. Returns a tuple of `( Burgers_Vector, Source_Location )`.
"""
rand_dislocation( rbv::RandomDislocation ) = ( rand_burger_vector(rbv), rand_dislocation_source(rbv) )

"""
    origin_already_present( idx, origin, all_dislocations )

Determine if the `origin isa Vector2D` is already present within the `all_dislocations` array up though the `idx` element using `Vectors2d.equal`.
"""
function origin_already_present( idx, origin, all_dislocations )
    if idx == 1
        return false
    end
    already_here = true
    for dis_idx ∈ 1:idx
        already_here = equal( origin, all_dislocations[dis_idx, Int(DislocationOrigin)] )
    end
    return true
end

"""
    collect_dislocations( rbv::RandomDislocation, num_dislocations ) -> Array{Vector2D, 2, num_dislocations}
"""
function collect_dislocations( rbv::RandomDislocation, num_dislocations )
    all_dislocations = Array{Vector2D{Float64}, 2}(undef, 2, num_dislocations)
    for idx ∈ eachindex(1:num_dislocations)
        (bob, origin) = rand_dislocation( rbv )
        # Be sure that the same origin is not used multiple times
        while origin_already_present( idx, origin, all_dislocations )
            origin = rand_dislocation( rbv )[ Int(DislocationOrigin) ]
        end
        all_dislocations[idx, Int(BurgersVector)]     = Vector2D( bob )
        all_dislocations[idx, Int(DislocationOrigin)] = Vector2D( origin )
    end
    return all_dislocations
end

"""
Struct containing uniform distribution for the Burger's vectors and their locations 

# Struct Members
* `axes::Tuple{Int}`: a tuple containing `axes[1] = Lx` and `axes[2] = Ly`.
* `burgers_vectors::AbstractArray`: the container with possible Burger's vectors 
* `random_location::Rand_Loc`: the `Distribution` for generating random locations 
* `random_burger_vector::Rand_BV`: the `Distribution` for generating random Burgers vectors 
"""
struct UniformBurgersVector {Rand_Loc <: Distribution, Rand_BV <: Distribution} <: RandomDislocation
    axes::Tuple{Int, Int}
    burgers_vectors::AbstractArray{Vector2D}
    random_location::Rand_Loc
    random_burger_vector::Rand_BV
end

"""
    UniformBurgersVector(; axes, burgers_vectors, random_location, random_burger_vector )

Keyword constructor for the `UniformBurgersVector` struct.
"""
UniformBurgersVector(; axes, burgers_vectors, random_location, random_burger_vector ) = UniformBurgersVector( axes, burgers_vectors, random_location, random_burger_vector ) 

"""
    rand_burger_vector( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` for a Burger's vector according to the `ubv <: UniformBurgersVector` allowed Burgers vectors.
"""
function rand_burger_vector( ubv::UniformBurgersVector )
    index = rand( ubv.random_burger_vector )
    return burgers_vectors[index]    
end

"""
    rand_dislocation_source( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` that denotes the center of a plaquette for a dislocation origin.
"""
rand_dislocation_source( ubv::UniformBurgersVector ) = Vector2D( ( 0.5 + rand(ubv.random_location), 0.5 + rand(ubv.random_location) ) )

"""
Interface type for the random strain distributions
"""
abstract type RandomStrainDistribution end

"""
Struct containing all information required for random strains generated by edge dislocations.

# Struct Members
* `concentration::AbstractFloat`: the concentration of dislocations in a 2D square lattice
* `vector_diff::Function`: the difference function used to for `Vector2D` types depending on boundary conditions
* `rand_num_dis::Dis`: the `Dis <: Distribution` that returns a random number of dislocations per disorder configuration
* `burgers_vectors::RBVD`: the distribution `RBVD <: RandomDislocation` of random Burger's vectors and locations in the `axes`
"""
struct RandomDislocationDistribution{Dis <: Distribution, RBVD <: RandomDislocation} <: RandomStrainDistribution
    concentration::AbstractFloat
    vector_diff::Function
    rand_num_dis::Dis
    burgers_vectors::RBVD
end

"""
    RandomDislocationDistribution(; axes, concentration, vector_diff, rand_num_dis, burgers_vectors )

Keyword constructor for the `RandomDislocationDistribution` struct.
"""
RandomDislocationDistribution(; concentration, vector_diff, rand_num_dis, burgers_vectors ) = RandomDislocationDistribution( concentration, vector_diff, rand_num_dis, burgers_vectors )

end # module RandomStrainDistributions
