module DisorderConfigurations

using StaticArrays
using ..PhysicalVectors
using ..ShearFunctions
using ..RandomStrainDistributions: BurgersVector, DislocationOrigin

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

mutable struct ShearFromDislocations{T <: Number} <: DisorderConfiguration
    vector_diff::Function
    axes::Tuple{Int, Int}
    dislocations::Array{Vector2D{T}, 2}
    strain_fields::Array{T, 3}
end

"""
    ShearFromDislocations(; Lx, Ly = Lx)

Convenient keyword constructor for the `ShearFromDislocations` type. 
"""
function ShearFromDislocations{T}(; Lx, Ly = Lx, diff ) where T
    axes = (Lx, Ly)
    dislocations = Array{Vector2D{T}, 2}(undef, (0, 0))
    return ShearFromDislocations( diff, axes, dislocations, zeros(T, (axes..., 2)) )
end

ShearFromDislocations(; kwargs...) = ShearFromDislocations{Float64}(; kwargs...)

set_dislocations!( sfd::ShearFromDislocations, dislocation_property_array ) = ( sfd.dislocations = dislocation_property_array )

function compute_strains!( sfd::ShearFromDislocations )
    temp_eval_r = Vector2D(0., 0.)
    for (xdx, ydx) ∈ collect( Iterators.product( 1:sfd.axes[1], 1:sfd.axes[2] ) )
        for dis_idx ∈ 1:size(sfd.dislocations)[2]
            Vector2D!( temp_eval_r, (xdx, ydx) )
            (b1g, b2g) = bxg_shears!( temp_eval_r, sfd.dislocations[Int(BurgersVector), dis_idx]; 
                                      source_r = sfd.dislocations[Int(DislocationOrigin), dis_idx], diff = sfd.vector_diff )
            sfd.strain_fields[xdx, ydx, 1] += b1g
            sfd.strain_fields[xdx, ydx, 2] += b2g
        end
    end
    return nothing
end

function generate_disorder!( disorder_field, sfd::ShearFromDislocations )
    compute_strains!( sfd )
    copyto!( disorder_field, sfd.strain_fields )
    return nothing 
end

function generate_disorder!( disorder_field, Lx::Int, Ly::Int, dislocations, diff = (x, y) -> subtract_PBC!(x, y, Lx, Ly) )
    sfd = ShearFromDislocations(; Lx, Ly, diff)
    set_dislocations!(sfd, dislocations)
    generate_disorder!(disorder_field, sfd)
end

make_fields(sfd::ShearFromDislocations{T}) where T = zeros(T, (sfd.axes..., 2))

function generate_disorder( sfd::ShearFromDislocations )
    disorder_fields = make_fields(sfd)
    generate_disorder!( disorder_fields, sfd )
    return disorder_fields
end

function generate_disorder(Lx::Int, Ly::Int, dislocations, diff = (x, y) -> subtract_PBC!(x, y, Lx, Ly))
    sfd = ShearFromDislocations(; Lx, Ly, diff)
    set_dislocations!(sfd, dislocations)
    disorder_fields = make_fields(sfd)
    generate_disorder!(disorder_fields, sfd)
    return disorder_fields
end

end # module DisorderConfigurations
