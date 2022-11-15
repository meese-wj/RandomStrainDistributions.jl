module CrystalDefects

using ..PhysicalVectors
import Base: eltype

export CrystalDefect, Dislocation, Dislocation2D, burgers_vector, dislocation_origin, eltype

"""
    abstract type CrystalDefect end

Parent type for all crystal defects.
"""
abstract type CrystalDefect end
"""
    abstract type Dislocation <: CrystalDefect end

Parent type for all dislocation [`CrystalDefect`](@ref)s.

# Required methods

* `burgers_vector`: return a copy of the Dislocation Burgers vector
* `dislocation_origin`: return a copy of the Dislocation origin.

"""
abstract type Dislocation <: CrystalDefect end
burgers_vector(dis::Dislocation) = throw(MethodError(burgers_vector, (dis,)))
dislocation_origin(dis::Dislocation) = throw(MethodError(dislocation_origin, (dis,)))
Base.eltype(dis::Dislocation) = Base.eltype(burgers_vector(dis))

struct Dislocation2D{T <: AbstractFloat} <: Dislocation
    Bvector::Vector2D{T}
    Rvector::Vector2D{T}
end

burgers_vector( dis::Dislocation2D ) = dis.Bvector
dislocation_origin( dis::Dislocation2D ) = dis.Rvector

    
end # module CrystalDefects