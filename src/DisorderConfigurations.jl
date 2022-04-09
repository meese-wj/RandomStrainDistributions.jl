module DisorderConfigurations

using StaticArrays
using ..PhysicalVectors
using ..ShearFunctions

export DisorderConfiguration, generate_disorder, generate_disorder!

"""
Interface type for different disorder configurations
"""
abstract type DisorderConfiguration end

"""
    generate_disorder!( disorder_field, dis_config::DisorderConfiguration )

Interface function for the `DisorderConfiguration` type with an implication of an in-place modification of the `disorder_field`.
"""
generate_disorder!( disorder_field, dis_config::DisorderConfiguration ) = error("No implementation has been defined for $(typeof(dis_config)).")

"""
    generate_disorder( dis_config::DisorderConfiguration )

Interface function that creates new field to store the result of `generate_disorder!` in.
"""
generate_disorder( dis_config::DisorderConfiguration ) = error("No implementation has been defined for $(typeof(dis_config)).")

mutable struct ShearFromDislocations <: DisorderConfiguration
    vector_diff::Function
    axes::Tuple{Int, Int}
    dislocations::Array{Vector2D, 2}
    strain_fields::SVector{SMatrix}
end

"""
    ShearFromDislocations(; Lx, Ly = Lx)

Convenient keyword constructor for the `ShearFromDislocations` type. 
"""
function ShearFromDislocations(; Lx, Ly = Lx, diff )
    axes = (Lx, Ly)
    dislocations = Array{Vector2D, 2}(undef)
    b1g_field = @SMatrix zeros(Float64, ( Lx, Ly ))
    b2g_field = @SMatrix zeros(Float64, ( Lx, Ly ))
    return ShearFromDislocations( diff, axes, dislocations, @SVector [ b1g_field, b2g_field ] )
end

set_dislocations!( sfd::ShearFromDislocations, dislocation_property_array ) = ( sfd.dislocations = dislocation_property_array )


function compute_strains!( sfd::ShearFromDislocations )
    temp_eval_r = Vector2D(0., 0.)
    for (xdx, ydx) ∈ collect( Iterators.product( 1:sfd.Lx, 1:sfd.Ly ) )
        for dis_idx ∈ 1:size(sfd.dislocations[1])
            Vector2D!( temp_eval_r, (xdx, ydx) )
            (b1g, b2g) = bxg_shears!( temp_eval_r, sfd.dislocations[1, dis_idx]; source_r = sfd.dislocations[2, dis_idx], diff = sfd.vector_diff )
            sfd.strain_fields[1][xdx, ydx] += b1g
            sfd.strain_fields[2][xdx, ydx] += b2g
        end
    end
    return nothing
end

function generate_disorder!( disorder_field, sfd::ShearFromDislocations )
    compute_strains!( sfd )
    copyto!( disorder_field, sfd.strain_fields )
    return nothing 
end

function generate_disorder( sfd::ShearFromDislocations )
    disorder_fields = @SVector [ @SMatrix zeros(Float64, sfd.axes), @SMatrix zeros(Float64, sfd.axes) ]
    generate_disorder!( disorder_fields, sfd )
    return disorder_fields
end
 
end