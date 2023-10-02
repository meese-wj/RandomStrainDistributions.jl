using Distributions
using ..PhysicalVectors
using ..CrystalDefects

export RandomDislocation, UniformBurgersVector, tetragonal_burgers_vectors, diagonal_burgers_vectors, combined_burgers_vectors

"""
Interface type for the joint-distribution for the dislocations.

# Additional Information
For any given dislocation, it has a random Burgers vector and random location in the lattice. 
This type should define both additional distributions for the dislocations.
"""
abstract type RandomDislocation end

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
    rand_dislocation( rbv::UniformBurgersVector, [::Type{T <: Dislocation} = Dislocation2D{Float64}] )

Generate a random dislocation by defining its Burgers vector and origin. Returns a [`Dislocation2D`](@ref).
"""
rand_dislocation(rng,  rbv::RandomDislocation, ::Type{T} = Dislocation2D{Float64} ) where T <: Dislocation = T( rand_burger_vector(rng, rbv), rand_dislocation_source(rng, rbv) )

"""
    origin_already_present( origin, all_dislocations )

Determine if the `origin isa Vector2D` is already present within the `all_dislocations` array using `∈`.
"""
function origin_already_present( origin, all_dislocations )
    return origin ∈ dislocationorigin.(all_dislocations)
end

function unique_origin(origin, all_dislocations, rbv, rng, ::Type{T} = Dislocation2D{Float64}) where T <: Dislocation
    # Be sure that the same origin is not used multiple times
    new_origin::typeof(origin) = origin
    while origin_already_present( new_origin, all_dislocations )
        new_origin = dislocationorigin( rand_dislocation(rng, rbv, T ) )
    end
    return new_origin
end

function set_dislocation!(all_dislocations, idx, dis::T, rbv::RandomDislocation, rng) where T <: Dislocation
    origin = unique_origin(dislocationorigin(dis), all_dislocations, rbv, rng, T)
    all_dislocations[idx] = T( burgersvector(dis), dislocationorigin(dis) )
    return nothing
end

function generate_dislocation!(all_dislocations, idx, rbv::RandomDislocation, rng, ::Type{T} = Dislocation2D{Float64}) where T <: Dislocation
    dis = rand_dislocation( rng, rbv, T )
    set_dislocation!(all_dislocations, idx, dis, rbv, rng)    
    return nothing
end

"""
    collect_dislocations( rbv::RandomDislocation, num_dislocations ) -> Vector{Dislocation}(num_dislocations}

Generate the set of dislocations. Only the first `num_dislocations ÷ 2` are *unique*. The last 
`num_dislocations ÷ 2` always have opposite Burgers vectors to the first set such that the 
total topological charge is zero.
"""
function collect_dislocations( rng, rbv::RandomDislocation, num_dislocations, ::Type{T} = Dislocation2D{Float64}) where T <: Dislocation
    all_dislocations = Vector{T}(undef, num_dislocations)
    unique_dislocations = num_dislocations ÷ 2
    for idx ∈ eachindex(1:unique_dislocations)
        generate_dislocation!(all_dislocations, idx, rbv, rng, T)
    end

    for (old, new) ∈ zip( 1:unique_dislocations, (unique_dislocations + one(Int):num_dislocations) )
        old_burger = burgersvector(all_dislocations[old])
        new_burger = -1 * old_burger
        new_origin = unique_origin( dislocationorigin(all_dislocations[old]), all_dislocations, rbv, rng, T )
        set_dislocation!(all_dislocations, new, T(new_burger, new_origin), rbv, rng)
    end
    return all_dislocations
end

const tetragonal_burgers_vectors = ( Vector2D(1., 0), Vector2D(-1., -0.), Vector2D(0, 1.), Vector2D(-0., -1.) )
const diagonal_burgers_vectors = (
	Vector2D(1/sqrt(2), 1/sqrt(2)),
	Vector2D(-1/sqrt(2), -1/sqrt(2)),
	Vector2D(-1/sqrt(2), 1/sqrt(2)),
	Vector2D(1/sqrt(2), -1/sqrt(2)),
)

const combined_burgers_vectors = (tetragonal_burgers_vectors..., diagonal_burgers_vectors...)

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
    UniformBurgersVector(; Lx, Ly = Lx, burgers_vectors )

Convenient keyword constructor for the `UniformBurgersVector` struct which makes sure the distributions agree with other members.
"""
function UniformBurgersVector(; Lx, Ly = Lx, burgers_vectors )
    axes = (Lx, Ly)
    bv = [ bob for bob ∈ burgers_vectors ]
    rand_pos = ( DiscreteUniform( 1, Lx ), DiscreteUniform(1, Ly) )
    rand_bv = DiscreteUniform(1, length(burgers_vectors))
    return UniformBurgersVector( axes, 
                                 bv,
                                 rand_pos,
                                 rand_bv )
end

"""
    rand_burger_vector( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` for a Burger's vector according to the `ubv <: UniformBurgersVector` allowed Burgers vectors.
"""
function rand_burger_vector(rng, ubv::UniformBurgersVector )
    index = rand(rng, ubv.random_burger_vector )
    return ubv.burgers_vectors[index]    
end

"""
    rand_dislocation_source( ubv::UniformBurgersVector ) -> Vector2D{Float64}

Returns a random `Vector2D` that denotes the center of a plaquette for a dislocation origin.
"""
rand_dislocation_source(rng, ubv::UniformBurgersVector ) = Vector2D( map( x -> 0.5 + rand(rng, x), ubv.random_position ) )

"""
    system_size( ubv::RandomDislocationDistribution )

Return the `prod`uct over the `axes` member `Tuple`.
"""
system_size( ubv::UniformBurgersVector ) = prod(ubv.axes)