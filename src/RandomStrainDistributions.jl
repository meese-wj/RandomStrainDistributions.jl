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
export RandomDislocation, UniformBurgersVector, RandomStrainDistribution, RandomDislocationDistribution, collect_dislocations

"""
Interface type for the joint-distribution for the dislocations.

# Additional Information
For any given dislocation, it has a random Burgers vector and random location in the lattice. 
This type should define both additional distributions for the dislocations.
"""
abstract type RandomDislocation end

@enum DislocationProperties begin 
    BurgersVector = 1
    DislocationOrigin
end

"""
    system_size( rbv::RandomDislocation )

Interface function for the `RandomDislocation` type.
"""
system_size( rbv::RandomDislocation ) = error("No implementation defined for $(typeof(rbv)).")

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

Determine if the `origin isa Vector2D` is already present within the `all_dislocations` array up until the `idx` element using `Vectors2d.equal`.
"""
function origin_already_present( idx, origin, all_dislocations )
    if idx == 1
        return false
    end
    already_here = false
    for dis_idx ∈ 1:(idx-1)
        already_here = equal( origin, all_dislocations[Int(DislocationOrigin), dis_idx] )
    end
    return already_here
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
        all_dislocations[Int(BurgersVector), idx]     = Vector2D( bob )
        all_dislocations[Int(DislocationOrigin), idx] = Vector2D( origin )
    end
    return all_dislocations
end

"""
Struct containing uniform distribution for the Burger's vectors and their locations 

# Struct Members
* `axes::Tuple{Int}`: a tuple containing `axes[1] = Lx` and `axes[2] = Ly`.
* `burgers_vectors::Vector`: the container with possible Burger's vectors 
* `random_position::Tuple{Rand_Loc, Rand_Loc}`: a `Tuple` of `Distribution`s for generating random dislocation positions 
* `random_burger_vector::Rand_BV`: the `Distribution` for generating random Burgers vectors 
"""
struct UniformBurgersVector{Rand_Loc <: Distribution, Rand_BV <: Distribution} <: RandomDislocation
    axes::Tuple{Int, Int}
    burgers_vectors::Vector{Vector2D{Float64}}
    random_position::Tuple{Rand_Loc, Rand_Loc}
    random_burger_vector::Rand_BV
end

"""
    UniformBurgersVector(; axes, burgers_vectors, random_position, random_burger_vector )

Keyword constructor for the `UniformBurgersVector` struct.
"""
# UniformBurgersVector(; axes, burgers_vectors, random_position, random_burger_vector ) = UniformBurgersVector( axes, burgers_vectors, random_position, random_burger_vector ) 

"""
    UniformBurgersVector(; Lx, Ly = Lx, burgers_vectors )

Convenient keyword constructor for the `UniformBurgersVector` struct which makes sure the distributions agree with other members.
"""
function UniformBurgersVector(; Lx, Ly = Lx, burgers_vectors )
    axes = (Lx, Ly)
    rand_pos = ( DiscreteUniform( 1, Lx ), DiscreteUniform(1, Ly) )
    rand_bv = DiscreteUniform(1, length(burgers_vectors))
    return UniformBurgersVector( axes, 
                                 burgers_vectors,
                                 rand_pos,
                                 rand_bv )
end

"""
    rand_burger_vector( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` for a Burger's vector according to the `ubv <: UniformBurgersVector` allowed Burgers vectors.
"""
function rand_burger_vector( ubv::UniformBurgersVector )
    index = rand( ubv.random_burger_vector )
    return ubv.burgers_vectors[index]    
end

"""
    rand_dislocation_source( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` that denotes the center of a plaquette for a dislocation origin.
"""
rand_dislocation_source( ubv::UniformBurgersVector ) = Vector2D( ( 0.5 + rand(ubv.random_position[1]), 0.5 + rand(ubv.random_position[2]) ) )

"""
    system_size( ubv::RandomDislocationDistribution )

Return the `prod`uct over the `axes` member `Tuple`.
"""
system_size( ubv::UniformBurgersVector ) = prod(ubv.axes)

"""
Interface type for the random strain distributions
"""
abstract type RandomStrainDistribution end

"""
    collect_dislocations( rsd::RandomStrainDistribution )

Infterface function to generate dislocations from a `RandomStrainDistribution`.
"""
collect_dislocations( rsd::RandomStrainDistribution ) = error("No implementation for $(typeof(rsd)) has been defined.")

"""
    system_size( rsd::RandomStrainDistribution )

Interface function which returns the system size of a `RandomStrainDistribution`.
"""
system_size( rsd::RandomStrainDistribution ) = error("No implementation for $(typeof(rsd)) has been defined.")

# function calculate_strains( eval_r::Vector2D, rsd::RandomStrainDistribution )
#     all_dislocations = collect_dislocations( rsd )
#     axes = system_size( rsd )
#     num_dis = length( all_dislocations[:, Int(BurgersVector)] )
#     b1g = zero(eltype(eval_r.vec))
#     b2g = zero(eltype(eval_r.vec))
#     temp_eval_r = Vector2D( eval_r )
#     for dis_idx ∈ 1:num_dis
#         Vector2D!( temp_eval_r, eval_r )
#         new_strains = bxg_shears!( temp_eval_r, all_dislocations[dis_idx, Int(BurgersVector)];
#                                    source_r = all_dislocations[dis_idx, Int(DislocationOrigin) ],
#                                    diff = (A, B) -> subtract_PBC!( A, B; Lx = axes[1], Ly = axes[2] ) )
#         b1g += new_strains[1]
#         b2g += new_strains[2]
#     end
#     return ( b1g, b2g )
# end

"""
Struct containing all information required for random strains generated by edge dislocations.

# Struct Members
* `concentration::AbstractFloat`: the concentration of dislocations in a 2D square lattice
* `vector_diff::Function`: the difference function used to for `Vector2D` types depending on boundary conditions
* `rand_num_dis::Dis`: the `Dis <: Distribution` that returns a random number of dislocations per disorder configuration
* `burgers_vector_dist::RBVD`: the distribution `RBVD <: RandomDislocation` of random Burger's vectors and locations in the `axes`
"""
struct RandomDislocationDistribution{Dis <: Distribution, RBVD <: RandomDislocation} <: RandomStrainDistribution
    concentration::AbstractFloat
    vector_diff::Function
    rand_num_dis::Dis
    burgers_vector_dist::RBVD
end

"""
    RandomDislocationDistribution(; axes, concentration, vector_diff, rand_num_dis, burgers_vectors )

Keyword constructor for the `RandomDislocationDistribution` struct.
"""
# RandomDislocationDistribution(; concentration::AbstractFloat, vector_diff, rand_num_dis, burgers_vector_dist ) = RandomDislocationDistribution( concentration, vector_diff, rand_num_dis, burgers_vector_dist )

"""
    RandomDislocationDistribution(; expected_num_dislocations::Int, Lx, Ly = Lx, vector_diff = subtract_PBC!, burgers_vectors )

Convenient keyword constructor for the `RandomDislocationDistribution` struct which makes sure the distributions agree with other members.

# Additional Information
* This keyword constructor assumes that the `RandomDislocation` is a `UniformBurgersVector` and that the `rand_num_dis` is a `Binomial`.
"""
function RandomDislocationDistribution(; expected_num_dislocations::Int, Lx, Ly = Lx, vector_diff = subtract_PBC!, burgers_vectors )
    num_plaquettes = Lx * Ly 
    concentration = convert(Float64, expected_num_dislocations / num_plaquettes)
    diff = (A, B) -> vector_diff( A, B; Lx = Lx, Ly = Ly )
    if vector_diff != subtract_PBC!
        diff = vector_diff
    end 
    rand_num = Binomial( num_plaquettes, concentration )
    bv_dist = UniformBurgersVector(; Lx = Lx, Ly = Ly, burgers_vectors = burgers_vectors )
    return RandomDislocationDistribution(  concentration,
                                           diff,
                                           rand_num,
                                           bv_dist )
end

"""
    system_size( rdd::RandomDislocationDistribution )

Alias call to the `RandomDislocation` `system_size` function.
"""
system_size( rdd::RandomDislocationDistribution ) = system_size( rdd.burgers_vector_dist )

"""
    collect_dislocations( rdd::RandomDislocationDistribution )

Alias call to the `RandomDislocation` `collect_dislocations` function.
"""
function collect_dislocations( rdd::RandomDislocationDistribution )
    num_dis = rand(rdd.rand_num_dis)
    while num_dis == zero(typeof(num_dis)) || num_dis > system_size(rdd)
        num_dis = rand(rdd.rand_num_dis)
    end
    return collect_dislocations( rdd.burgers_vector_dist, num_dis )
end

end # module RandomStrainDistributions
