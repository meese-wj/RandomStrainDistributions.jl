module DisorderConfigurations

using StaticArrays
using ..PhysicalVectors
using ..CrystalDefects
using ..PBCImageFields
using ..ShearFunctions
using ..RandomStrainDistributions: burgersvector, dislocationorigin

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
    const include_Δ::Bool
    const coupling_ratio::T
    const tolerance::T
    axes::Tuple{Int, Int}
    dislocations::Vector{Dislocation2D{T}}
    strain_fields::Array{T, 3}
end

_field_columns(b::Bool) = b ? 3 : 2

"""
    ShearFromDislocations(; Lx, Ly = Lx)

Convenient keyword constructor for the `ShearFromDislocations` type. 
"""
function ShearFromDislocations{T}( Lx, Ly, include_Δ::Bool = false, tolerance::T = sqrt(eps()), coupling_ratio::T = oneunit(T) ) where T
    axes = (Lx, Ly)
    dislocations = Vector{Dislocation2D{T}}[]
    return ShearFromDislocations( include_Δ, coupling_ratio, tolerance, axes, dislocations, zeros(T, (axes..., _field_columns(include_Δ))) )
end

ShearFromDislocations(args...) = ShearFromDislocations{Float64}(args...)

set_dislocations!( sfd::ShearFromDislocations, dislocation_property_array ) = ( sfd.dislocations = dislocation_property_array )

function compute_strains!( sfd::ShearFromDislocations )
    include_Δ = sfd.include_Δ
    for (xdx, ydx) ∈ collect( Iterators.product( UnitRange(1, sfd.axes[1]), UnitRange(1, sfd.axes[2]) ) )
        for disloc ∈ sfd.dislocations
            temp_r = Vector2D( Float64(xdx), Float64(ydx) )
            b1g = PBCField_value( b1g_shear, temp_r, disloc, sfd.axes..., sfd.tolerance )
            b2g = PBCField_value( b2g_shear, temp_r, disloc, sfd.axes..., sfd.tolerance )
            sfd.strain_fields[xdx, ydx, 1] += b1g
            sfd.strain_fields[xdx, ydx, 2] += b2g
        end

        if include_Δ
            sfd.strain_fields[xdx, ydx, 3] = Δsplitting( sfd.strain_fields[xdx, ydx, 1], sfd.strain_fields[xdx, ydx, 2], sfd.coupling_ratio )
        end
    end
    return nothing
end

function generate_disorder!( disorder_field, sfd::ShearFromDislocations )
    compute_strains!( sfd )
    copyto!( disorder_field, sfd.strain_fields )
    disorder_field = sfd.strain_fields
    return nothing 
end

function generate_disorder!( disorder_field, Lx::Int, Ly::Int, dislocations; include_Δ = false, tolerance = sqrt(eps()), coupling_ratio = 1.0 )
    sfd = ShearFromDislocations( Lx, Ly, include_Δ, tolerance, coupling_ratio )
    set_dislocations!(sfd, dislocations)
    generate_disorder!(disorder_field, sfd)
end

make_fields(sfd::ShearFromDislocations{T}) where T = zeros(T, (sfd.axes..., _field_columns(sfd.include_Δ)) )

function generate_disorder( sfd::ShearFromDislocations )
    disorder_fields = make_fields(sfd)
    generate_disorder!( disorder_fields, sfd )
    return disorder_fields
end

function generate_disorder(Lx::Int, Ly::Int, dislocations; include_Δ = false, tolerance = sqrt(eps()), coupling_ratio = 1.0 )
    sfd = ShearFromDislocations( Lx, Ly, include_Δ, tolerance, coupling_ratio )
    set_dislocations!(sfd, dislocations)
    disorder_fields = make_fields(sfd)
    generate_disorder!(disorder_fields, sfd)
    return disorder_fields
end

end # module DisorderConfigurations
