module RandomStrainDistributions
using Reexport

include("PhysicalVectors/PhysicalVectors.jl")
@reexport using .PhysicalVectors

include("CrystalDefects.jl")
@reexport using .CrystalDefects

include("PBCImageFields.jl")
@reexport using .PBCImageFields

include("ShearFunctions.jl")
@reexport using .ShearFunctions

include("RandomDislocations.jl")

include("DisorderConfigurations.jl")
@reexport using .DisorderConfigurations

# ===================================================================== #
#  RandomStrainDistributions module code is below
# ===================================================================== #

using Distributions
export RandomStrainDistribution, RandomDislocationDistribution, collect_dislocations

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
* `truncated::Bool`: whether one should check that the number of dislocations lie on the interval ``[1, system_size]``.
* `concentration::AbstractFloat`: the concentration of dislocations in a 2D square lattice
* `vector_diff::Function`: the difference function used to for `Vector2D` types depending on boundary conditions
* `rand_num_dis::Dis`: the `Dis <: Distribution` that returns a random number of dislocations per disorder configuration
* `burgers_vector_dist::RBVD`: the distribution `RBVD <: RandomDislocation` of random Burger's vectors and locations in the `axes`
"""
struct RandomDislocationDistribution{Dis <: Distribution, RBVD <: RandomDislocation} <: RandomStrainDistribution
    truncated::Bool
    concentration::AbstractFloat
    vector_diff::Function
    rand_num_dis::Dis
    burgers_vector_dist::RBVD
end

"""
    RandomDislocationDistribution(; expected_num_dislocations::Int, Lx, Ly = Lx, vector_diff = subtract_PBC, burgers_vectors = tetragonal_burgers_vectors )

Convenient keyword constructor for the `RandomDislocationDistribution` struct which makes sure the distributions agree with other members.

# Additional Information
* This keyword constructor assumes that the `RandomDislocation` is a `UniformBurgersVector`.
* This does not use `StatsBase.Truncated` because it's about an order of magnitude slower than my version in `collect_dislocations`.
"""
function RandomDislocationDistribution(; expected_num_dislocations::Int, Lx, Ly = Lx, vector_diff = subtract_PBC, burgers_vectors = tetragonal_burgers_vectors, truncated = true )
    diff = (A, B) -> vector_diff( A, B, Lx, Ly )
    if vector_diff != subtract_PBC
        diff = vector_diff
    end 
    bv_dist = UniformBurgersVector(; Lx = Lx, Ly = Ly, burgers_vectors = burgers_vectors )

    num_plaquettes = system_size(bv_dist) 
    concentration = convert(Float64, expected_num_dislocations / num_plaquettes)
    rand_num = Binomial( num_plaquettes, concentration )

    return RandomDislocationDistribution(  truncated,
                                           concentration,
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

# Additional Information
* Note that this distribution will not allow for configurations with zero dislocations or for a number of dislocations larger than the `system_size`.
* The number of dislocations is constrained to be even only so as to have a net-zero topological charge.
* This keyword constructor assumes that the `rand_num_dis` is a truncated `Binomial` `Distribution` on the interval ``{2, 4, ..., Lx \\times Ly}`` when `truncated == true`.
"""
function collect_dislocations( rdd::RandomDislocationDistribution )
    num_dis = rand(rdd.rand_num_dis)
    while isodd(num_dis) 
        num_dis = rand(rdd.rand_num_dis) 
    end
    
    if rdd.truncated
        while num_dis == zero(typeof(num_dis)) || num_dis > system_size(rdd) || isodd(num_dis)
            num_dis = rand(rdd.rand_num_dis)
        end
    end
    return collect_dislocations( rdd.burgers_vector_dist, num_dis )
end

end # module RandomStrainDistributions
