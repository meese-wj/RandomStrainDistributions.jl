using DrWatson
import DrWatson: parse_savename
import Base: convert

export SimulationParameters, convert, parse_savename

abstract type SimulationParameters end

function Base.convert(::Type{T}, d::Union{Dict, NamedTuple}) where {T <: SimulationParameters}
    key_type = d isa Dict ? keytype(d) : Symbol
    keys = [ key_type(name) for name ∈ fieldnames(T) ]
    vals = [ d[key] for key ∈ keys ]
    return T(vals...)
end

DrWatson.parse_savename(::Type{T}, filename) where {T <: SimulationParameters} = convert(T, parse_savename(filename)[2]) 