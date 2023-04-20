using InteractiveUtils # for subtypes()

function _recursive_subtypes!(output, mytype)
    if isconcretetype(mytype)
        return output
    end
    st = subtypes(mytype)
    append!(output, st)
    for type ∈ st
        _recursive_subtypes!(output, type)
    end
end

function recursive_subtypes(mytype)
    output = []
    _recursive_subtypes!(output, mytype)
    return output
end