using StatsBase
using Statistics
using HypothesisTests
using Distributions

data_quantile(data, nbins) = quantile(data, LinRange(1/nbins, 1, nbins))

function histogram_fit(data, nbins)
    quants = data_quantile(data, nbins)
    hist = fit(Histogram, data, quants)
    return hist
end

χsqGoFTest(histogram_data::Histogram, expected_probabilities) = ChisqTest(histogram_data.weights, expected_probabilities)

function χsqGoFTest(data, nbins, fitted_dist::Distribution)
    # Make the histogram
    hist = histogram_fit(data, nbins)
    # Calculate expected_probabilities
    exp_prob = cdf.(fitted_dist, data_quantile(data, nbins)) |> diff # diff for the probability per bin in quants 
    exp_prob /= sum(exp_prob)
    return χsqGoFTest(hist, exp_prob)
end

