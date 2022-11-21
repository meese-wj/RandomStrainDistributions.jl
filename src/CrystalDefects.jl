module CrystalDefects

using ..PhysicalVectors
import Base: eltype

export CrystalDefect, Dislocation, Dislocation2D, burgersvector, dislocationorigin, shift_origin, eltype

"""
    abstract type CrystalDefect end

Parent type for all crystal defects.
"""
abstract type CrystalDefect end
"""
    abstract type Dislocation <: CrystalDefect end

Parent type for all dislocation [`CrystalDefect`](@ref)s.

# Required methods

* `burgersvector`: return a copy of the Dislocation Burgers vector
* `dislocationorigin`: return a copy of the Dislocation origin.

"""
abstract type Dislocation <: CrystalDefect end
burgersvector(dis::Dislocation) = throw(MethodError(burgersvector, (dis,)))
dislocationorigin(dis::Dislocation) = throw(MethodError(dislocationorigin, (dis,)))
shift_origin(dis::T, shift) where T = T( burgersvector(dis), dislocationorigin(dis) + shift )
Base.eltype(dis::Dislocation) = Base.eltype(burgersvector(dis))

struct Dislocation2D{T <: AbstractFloat} <: Dislocation
    Bvector::Vector2D{T}
    Rvector::Vector2D{T}
end

@inline burgersvector( dis::Dislocation2D ) = dis.Bvector
@inline dislocationorigin( dis::Dislocation2D ) = dis.Rvector

Dislocation2D{T}(::UndefInitializer) where T = Dislocation2D{T}( Vector2D{T}(undef), Vector2D{T}(undef) )

    
end # module CrystalDefects