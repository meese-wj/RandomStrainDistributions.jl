using StatsBase
using Statistics
using HypothesisTests
using Distributions
using QuadGK
using NumericalIntegration
using SpecialFunctions
using LsqFit

# data_quantile(data, nbins) = quantile(data, LinRange(1/nbins, 1, nbins))
data_quantile(data, nbins) = ( (minx, maxx) = extrema(data); LinRange(minx, maxx, nbins) )
data_range(data) = -((reverse ∘ extrema)(data)...)
bin_width(data, nbins) = data_range(data) / nbins
data_dx(data, nbins) = 1 / bin_width(data, nbins)
function bin_centers(hist) 
    xdata = hist.edges[1]
    δx = bin_width(xdata, length(xdata) - 1)
    xdata = collect(xdata[1:end-1])
    return xdata .+ 0.5 * δx
end
function normalized_weights(data, nbins, density = true)
    vals = fill( 1/length(data), length(data) )
    if density
        return data_dx(data, nbins) * vals
    end
    return vals
end

function histogram_fit(data, nbins; normalize = false, density = true)
    quants = data_quantile(data, nbins)
    hist = fit(Histogram, data, quants)
    if normalize
        weight_s = normalized_weights(data, nbins, density)
        hist = fit(Histogram, data, weights(weight_s), quants)
    end
    return hist
end

χsqGoFTest(histogram_data::Histogram, expected_probabilities) = ChisqTest(histogram_data.weights, expected_probabilities)

function χsqGoFTest(data, nbins, fitted_dist::Distribution)
    # Make the histogram
    hist = histogram_fit(data, nbins; normalize = false)
    # Calculate expected_probabilities
    exp_prob = cdf.(fitted_dist, data_quantile(data, nbins)) |> diff # diff for the probability per bin in quants 
    exp_prob /= sum(exp_prob)
    return χsqGoFTest(hist, exp_prob)
end

function χsqGoFTest(data, nbins, fitted_dist::Tuple{<: Function, <: AbstractVector})
    # Make the histogram
    hist = histogram_fit(data, nbins; normalize = false)
    # Calculate the expected_probabilities
    centers = bin_centers(hist)
    exp_prob = bin_width(hist.edges[1], length(hist.edges[1]) - 1) * distribution_function(fitted_dist, centers)
    exp_prob /= sum(exp_prob)
    return χsqGoFTest(hist, exp_prob)
end

cVTrange(Δdata, maxTratio = 10, minTratio = 1e-4, NT = 256) = ((Δmin, Δmax) = extrema(Δdata); LinRange(minTratio * Δmin, maxTratio * Δmax, NT))
schottkycV(Δ, T) = @. (( Δ / T ) * sech( Δ / T ))^2

distribution_function(fitted_distribution::Distribution, arg) = pdf(fitted_distribution, arg)
distribution_function(fitted_distribution::Tuple{<: Function, <:AbstractVector}, arg) = (fitted_distribution[1])(arg, fitted_distribution[2])

function fitcV(Δdata, fitted_distribution, Δmin = 0, Δmax = Inf )
    Trange = cVTrange(Δdata)
    integral = similar(Trange)
    if !(fitted_distribution isa Distribution)
        # This is because the Bessel function blows up in NaNs quickly
        Δmax = maximum(fitted_distribution[2]) * 10 
    end
    for (idx, Temp) ∈ enumerate(Trange)
        integral[idx] = try quadgk( Δ -> distribution_function(fitted_distribution, Δ) * schottkycV(Δ, Temp), Δmin, Δmax )[1]
        catch err
            err isa DomainError ? ( @show fitted_distribution) : nothing
            Inf
        end
    end
    return (Trange, integral)
end


function histcV(Δdata, nbins)
    hist = histogram_fit(Δdata, nbins; normalize = true, density = false)
    Trange = cVTrange(Δdata)
    cVvals = zeros(length(Trange))
    meanbins = diff(hist.edges[1]) / 2 .+ hist.edges[1][1:end-1]
    for (idx, Temp) ∈ enumerate(Trange), bin ∈ eachindex(hist.weights)
        cVvals[idx] += hist.weights[bin] * schottkycV(meanbins[bin], Temp)
    end
    return (Trange, cVvals)
end

mse( xdata, y1data, y2data ) = integrate( xdata, (@. (y1data - y2data)^2 ), Trapezoidal() ) / data_range(xdata)
rmse( xdata, y1data, y2data ) = mse(xdata, y1data, y2data) |> sqrt
mssum(xdata, y1data, y2data) = integrate( xdata, (@. 0.5 * (y1data^2 + y2data^2) ), Trapezoidal() ) / data_range(xdata)
similarity_coeff(xdata, y1data, y2data) = 1 - mse(xdata, y1data, y2data) / mssum(xdata, y1data, y1data)

# Two Uncorrelated Gaussians
gaussian_part(Δ, σ₁, σ₂) = @. exp( -0.25 * Δ^2 * (σ₁^2 + σ₂^2) / (σ₁ * σ₂)^2 ) 
bessel_arg(Δ, σ₁, σ₂) = @. -0.25 * Δ^2 * (σ₁^2 - σ₂^2) / (σ₁ * σ₂)^2
bessel_part(Δ, σ₁, σ₂) = @. besselix(0, bessel_arg(Δ, σ₁, σ₂)) * exp( (abs ∘ bessel_arg)(Δ, σ₁, σ₂) )
pseudo_rayleigh(Δ, σ₁, σ₂) = @. Δ / (σ₁ * σ₂) * gaussian_part(Δ, σ₁, σ₂) * bessel_part(Δ, σ₁, σ₂)
pseudo_rayleigh(Δ, sigmas) = pseudo_rayleigh(Δ, sigmas...)
pseudo_rayleigh_distribution(σ₁, σ₂) = (pseudo_rayleigh, [σ₁, σ₂])
pseudo_rayleigh_distribution(sigmas) = pseudo_rayleigh_distribution(sigmas...)

function fit_pseudo_rayleigh(xdata, ydata, sigma0, lb, ub)
    return curve_fit(pseudo_rayleigh, xdata, ydata, sigma0; lower = lb, upper = ub)
end

# Assumes the Δhist is density-normalized
function fit_pseudo_rayleigh(Δhist::Histogram, sigma0 = fill(5.0, 2), lb = fill(1e-2, 2), ub = fill(Inf, 2))
    # First get the xdata
    xdata = bin_centers(Δhist)
    # Now get the y values
    ydata = Δhist.weights
    mode = xdata[ argmax(Δhist.weights) ]
    return fit_pseudo_rayleigh(xdata, ydata, mode * sigma0, mode * lb, mode * ub)
end