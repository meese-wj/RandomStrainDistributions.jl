
using FFTW
using Statistics
using JLD2

function fourierspaceCF(realspace_array; shift = true)
    arr = realspace_array .- mean(realspace_array)
    # Compute conjugate(arr) * arr 
    ft = abs2.(fft(arr))
    ft = shift ? fftshift(ft) : ft
    return ft
end

# ∑ₓ arr(x) * arr(x + y)
function fourierspaceCCF(realspace_array1, realspace_array2; shift = true)
    arr1 = realspace_array1 .- mean(realspace_array1)
    arr2 = realspace_array2 .- mean(realspace_array2)
    # Compute conjugate(arr) .* arr 
    ft = conj!( fft(arr1) ) .* fft(arr2)
    ft = shift ? fftshift(ft) : ft
    return ft
end

function realspaceCF(realspace_array; shift = true)
    fsCF = fourierspaceCF(realspace_array; shift)
    fsCF = shift ? ifftshift(fsCF) : fsCF
    # Always shift this one because there's no need to go back
    rsCF = ifft(fsCF) ./ length(fsCF) |> fftshift
    return rsCF
end

# ∑ₓ arr(x) * arr(x + y)
function realspaceCCF(realspace_array1, realspace_array2; shift = true)
    fsCCF = fourierspaceCCF(realspace_array1, realspace_array2; shift)
    fsCCF = shift ? ifftshift(fsCCF) : fsCCF
    # Always shift this one because there's no need to go back
    rsCCF = ifft(fsCCF) ./ length(fsCCF) |> fftshift
    return rsCCF
end

function extract_single_Correlations(filename, param_type::Type{<: SimulationParameters} = DistributionParameters)
    results, params, colnames, colvals = read_results_parameters(filename, param_type)

    Lsize = params.Lx
    strains = results["strains"]
    nstrainfields = size(strains, 3)
    fsCF = zeros(ComplexF64, Lsize, Lsize, nstrainfields)
    rsCF = zeros(ComplexF64, Lsize, Lsize, nstrainfields)

    for trial_dx ∈ UnitRange(1, size(strains, 4)), strain_dx ∈ UnitRange(1, nstrainfields)
        @view(fsCF[:, :, strain_dx]) .+= fourierspaceCF( @view(strains[:, :, strain_dx, trial_dx]); shift = true )
        @view(rsCF[:, :, strain_dx]) .+= realspaceCF( @view(strains[:, :, strain_dx, trial_dx]); shift = true )
    end

    @view(fsCF[:, :, :]) ./= size(strains, 4)
    @view(rsCF[:, :, :]) ./= size(strains, 4)
    
    return DataFrame( [[fsCF[:, :, 1]], [fsCF[:, :, 2]], [fsCF[:, :, 3]], 
                       [rsCF[:, :, 1]], [rsCF[:, :, 2]], [rsCF[:, :, 3]]],
                      [:FourierSpaceCFB1g, :FourierSpaceCFB2g, :FourierSpaceCFΔ, 
                       :RealSpaceCFB1g, :RealSpaceCFB2g, :RealSpaceCFΔ] )
end

# Only find the cross correlations between B1g and B2g for now
function extract_cross_Correlations(filename, param_type::Type{<: SimulationParameters} = DistributionParameters)
    results, params, colnames, colvals = read_results_parameters(filename, param_type)

    Lsize = params.Lx
    strains = results["strains"]
    fsCCF = zeros(ComplexF64, Lsize, Lsize, 2)
    rsCCF = zeros(ComplexF64, Lsize, Lsize, 2)

    for trial_dx ∈ UnitRange(1, size(strains, 4))
        @view(fsCCF[:, :, 1]) .+= fourierspaceCCF( @view(strains[:, :, 1, trial_dx]), @view(strains[:, :, 2, trial_dx]); shift = true )
        @view(fsCCF[:, :, 2]) .+= fourierspaceCCF( @view(strains[:, :, 2, trial_dx]), @view(strains[:, :, 1, trial_dx]); shift = true )
        
        @view(rsCCF[:, :, 1]) .+= realspaceCCF( @view(strains[:, :, 1, trial_dx]), @view(strains[:, :, 2, trial_dx]); shift = true )
        @view(rsCCF[:, :, 2]) .+= realspaceCCF( @view(strains[:, :, 2, trial_dx]), @view(strains[:, :, 1, trial_dx]); shift = true )
    end

    @view(fsCCF[:, :, :]) ./= size(strains, 4)
    @view(rsCCF[:, :, :]) ./= size(strains, 4)
    
    return DataFrame( [[fsCCF[:, :, 1]], [fsCCF[:, :, 2]], 
                       [rsCCF[:, :, 1]], [rsCCF[:, :, 2]]],
                      [:FourierSpaceCCFB1gB2g, :FourierSpaceCCFB2gB1g, 
                       :RealSpaceCCFB1gB2g, :RealSpaceCCFB2gB1g] )
end
