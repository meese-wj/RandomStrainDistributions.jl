using DrWatson
import DrWatson: parse_savename
import Base: convert

export SimulationParameters, convert, parse_savename

abstract type SimulationParameters end

function add_defaults!(extracted_params, required_keys, defaults, key_type)
    extracted_keys = keys(extracted_params)
    if defaults isa Nothing
        isempty(setdiff(required_keys, extracted_keys)) ? nothing : throw(error("The required fields, $(required_keys), do not match the extracted fields, $(extracted_keys). The difference is:\n\t$(setdiff(required_keys, extracted_keys))."))
        return nothing
    end
    defaults = Dict( map( (key, val) -> (key_type(key), val), keys(defaults), values(defaults) ) )
    needed_keys = setdiff(required_keys, extracted_keys)
    isempty(setdiff(needed_keys, keys(defaults))) ? nothing : throw(error("The fields needed, $(needed_keys), do not match the supplied defaults, $(keys(defaults)). The difference is:\n\t$(setdiff(needed_keys, keys(defaults)))."))
    for (key, val) ∈ defaults
        extracted_params[key] = val
    end
    return
end

function Base.convert(::Type{T}, d::Union{Dict, NamedTuple}, defaults = nothing) where {T <: SimulationParameters}
    key_type = d isa Dict ? keytype(d) : Symbol
    field_keys = [ key_type(name) for name ∈ fieldnames(T) ]
    found_params = copy(d)
    add_defaults!(found_params, field_keys, defaults, key_type)
    vals = [ found_params[key] for key ∈ field_keys ]
    return T(vals...)
end

DrWatson.parse_savename(::Type{T}, filename; defaults = nothing) where {T <: SimulationParameters} = convert(T, parse_savename(filename)[2], defaults) 