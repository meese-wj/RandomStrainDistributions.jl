

struct DistributionParameters{T <: AbstractFloat} <: SimulationParameters
    Lx::Int
    Ly::Int
    ndislocations::Int
    cratio::T
    rtol::T
    nsamples::Int
end

function DistributionParameters(; Lx = 128, Ly = Lx, ndislocations = 2, cratio = 1.0, rtol = 0.1, nsamples = 1)
    return DistributionParameters{typeof(cratio)}(Lx, Ly, ndislocations, cratio, rtol, nsamples)
end

nsize(params::DistributionParameters) = params.Lx * params.Ly
concentration(params::DistributionParameters) = params.ndislocations / nsize(params)