#< FITNESS FUNCTION HELPER FUNCTIONS ##
"""Get summed differences of the peak values from the FFT of the solution"""
function getDif(peakvals)
    differences = diff(peakvals)
    sum(abs.(differences))#/length(peakvals)
end

"""Get average standard deviation of peaks values from the FFT of the solution"""
function getSTD(fft_peakindxs::Vector{Int}, fftData, window::Int =1) 
    arrLen = length(fftData)

    sum_std = sum(std(@view fftData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs; init=0.0) #* sum rolling window of standard deviations

    return sum_std #/ length(fft_peakindxs) #* divide by number of peaks to get average std
end 

function getSTD_vector(fft_peakindxs::Vector{Int}, fft_arrayData, window::Int =1) 
    arrLen = length(fft_arrayData)

    [std(@view fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs]
end 
#> END OF FITNESS FUNCTION HELPER FUNCTIONS ##



#< FFT HELPER FUNCTIONS ##
"""
    getFrequencies(timeseries)
Return the real-valued FFT of a 1D ODE solution, will be half the length of the timeseries
"""
function getFrequencies(timeseries::AbstractVector{Float64})
    # rfft_result = rfft(timeseries)
    rfft_result = RFFT_PLAN * timeseries
    norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    abs.(rfft_result) ./ norm_val
end

"""
    getFrequencies!(fft_array, timeseries)
Computes the real-valued FFT and returns it in-place to the preallocated fft_array, which is half the length of `timeseries`.
"""
function getFrequencies!(fft_array, timeseries) 
    rfft_result = rfft(timeseries)
    norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    fft_array .= abs.(rfft_result) ./ norm_val
end

# function getFrequencies!(timeseries, jump::Int = 1) 
#     rfft_result = rfft(@view timeseries[1:jump:end])
#     norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
#     timeseries[1:length(rfft_result)] .= abs.(rfft_result) ./ norm_val
# end

"""Normalize FFT array in-place to have mean 0 and amplitude 1"""
function normalize_time_series!(fftarray)
    mu = mean(fftarray)
    amplitude = maximum(fftarray) - minimum(fftarray)
    fftarray .= (fftarray .- mu) ./ amplitude
end
#> END OF FFT HELPER FUNCTIONS ##



#<< PERIOD AND AMPLITUDE FUNCTIONS ##
# function getPerAmp(Amem_sol, solt::T) where T

#     indx_max, vals_max, indx_min, vals_min = findextrema(Amem_sol; min_height=0.1)
#     if length(indx_max) < 2 || length(indx_min) < 2
#         #* if there are not enough peaks in the time domain, return just the fitness
#         return 0.0, 0.0
#     end

#     return getPerAmp(solt, indx_max, vals_max, indx_min, vals_min)
# end

"""Calculates the period and amplitude of each individual in the population."""
function getPerAmp(solt, indx_max::Vector{Int}, vals_max::AbstractVector{Float64}, vals_min::AbstractVector{Float64}) 

    #* Calculate amplitudes and periods
    # pers = (solt[indx_max[i+1]] - solt[indx_max[i]] for i in 1:(length(indx_max)-1))
    # amps = (vals_max[i] - vals_min[i] for i in 1:min(length(indx_max), length(indx_min)))

    # return mean(pers), mean(amps) 
    return compute_period(solt, indx_max), compute_amplitude(vals_max, vals_min)
end

function compute_period(solt, indx_max::Vector{Int})
    mean(solt[indx_max[i+1]] - solt[indx_max[i]] for i in 1:(length(indx_max)-1))
end

function compute_amplitude(vals_max::AbstractVector{Float64}, vals_min::AbstractVector{Float64})
    mean(abs, (vals_max[i] - vals_min[i] for i in 1:min(length(vals_max), length(vals_min))))
end
#> END OF PERIOD AND AMPLITUDE FUNCTIONS ##






















