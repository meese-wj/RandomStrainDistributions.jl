
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

function realspaceCF(realspace_array; shift = true)
    fsCF = fourierspaceCF(realspace_array; shift)
    fsCF = shift ? ifftshift(fsCF) : fsCF
    # Always shift this one because there's no need to go back
    rsCF = ifft(fsCF) ./ length(fsCF) |> fftshift
    return rsCF
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