function jobname_parser(jobname, parse_variable, parse_type; delimiter = "_", connector = "=")
    (contains(jobname, parse_variable)) ? nothing : throw(ArgumentError("Job name $jobname does not contain the variable $parse_variable."))
    var = parse_variable * connector
    var_range = findfirst(var, jobname)
    var_idx_start = var_range[end] + one(Int)
    jobname[var_idx_start:end]
    value_range = findfirst(delimiter, jobname[var_idx_start:end])
    if value_range === nothing
        var_idx_end = length(jobname)
    else
        var_idx_end = var_idx_start + value_range[end] - 2 * one(Int)
    end
    return parse(parse_type, jobname[var_idx_start:var_idx_end])
end

