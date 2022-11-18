module CrystalDefects

using ..PhysicalVectors
import Base: eltype

export CrystalDefect, Dislocation, Dislocation2D, BurgersVector, DislocationOrigin, shift_origin, eltype

"""
    abstract type CrystalDefect end

Parent type for all crystal defects.
"""
abstract type CrystalDefect end
"""
    abstract type Dislocation <: CrystalDefect end

Parent type for all dislocation [`CrystalDefect`](@ref)s.

# Required methods

* `BurgersVector`: return a copy of the Dislocation Burgers vector
* `DislocationOrigin`: return a copy of the Dislocation origin.

"""
abstract type Dislocation <: CrystalDefect end
BurgersVector(dis::Dislocation) = throw(MethodError(BurgersVector, (dis,)))
DislocationOrigin(dis::Dislocation) = throw(MethodError(DislocationOrigin, (dis,)))
shift_origin(dis::T, shift) where T = T( BurgersVector(dis), DislocationOrigin(dis) + shift )
Base.eltype(dis::Dislocation) = Base.eltype(BurgersVector(dis))

struct Dislocation2D{T <: AbstractFloat} <: Dislocation
    Bvector::Vector2D{T}
    Rvector::Vector2D{T}
end

BurgersVector( dis::Dislocation2D ) = dis.Bvector
DislocationOrigin( dis::Dislocation2D ) = dis.Rvector

Dislocation2D{T}(::UndefInitializer) where T = Dislocation2D{T}( Vector2D{T}(undef), Vector2D{T}(undef) )

    
end # module CrystalDefects