#< FITNESS FUNCTION HELPER FUNCTIONS ##
"""Get average difference of the first and last peak values from the FFT of the solution"""
function getDif(peakvals::Vector{Float64})
    # (peakvals[begin] - peakvals[end])/length(peakvals)
    sum(abs.(diff(peakvals)))/length(peakvals)
end

"""Get summed average standard deviation of peaks values from the FFT of the solution"""
function getSTD(fft_peakindxs::Vector{Int}, fft_arrayData, window::Int =1) 
    arrLen = length(fft_arrayData)

    sum_std = sum(std(@view fft_arrayData[max(1, ind - window):min(arrLen, ind + window)]) for ind in fft_peakindxs; init=0.0) #* sum rolling window of standard deviations

    return sum_std #/ length(fft_peakindxs) #* divide by number of peaks to get average std
end 
#> END OF FITNESS FUNCTION HELPER FUNCTIONS ##



#< FFT HELPER FUNCTIONS ##
"""
    getFrequencies(timeseries)
Return the real-valued FFT of a timeseries, will be half the length of the timeseries
"""
function getFrequencies(timeseries::Vector{Float64})#, rfft_plan)
    # sampled_timeseries = @view timeseries[1:jump:end]
    rfft_result = rfft(timeseries)
    # rfft_result = rfft_plan * timeseries
    norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    abs.(rfft_result) ./ norm_val
end

"""
    getFrequencies!(fft_array, timeseries)
Computes the real-valued FFT and returns it in-place to the preallocated fft_array, which is half the length of `timeseries`.
"""
function getFrequencies!(fft_array, timeseries::Vector{Float64}, jump::Int = 1) 
    rfft_result = rfft(@view timeseries[1:jump:end])
    norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    fft_array[1:length(rfft_result)] .= abs.(rfft_result) ./ norm_val
end

function getFrequencies!(timeseries::Vector{Float64}, jump::Int = 1) 
    rfft_result = rfft(@view timeseries[1:jump:end])
    norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    timeseries[1:length(rfft_result)] .= abs.(rfft_result) ./ norm_val
end

"""Normalize FFT array in-place to have mean 0 and amplitude 1"""
function normalize_time_series!(fftarray)
    mu = mean(fftarray)
    amplitude = maximum(fftarray) - minimum(fftarray)
    fftarray .= (fftarray .- mu) ./ amplitude
end
#> END OF FFT HELPER FUNCTIONS ##



#<< PERIOD AND AMPLITUDE FUNCTIONS ##
function getPerAmp(Amem_sol, solt::T) where T

    indx_max, vals_max, indx_min, vals_min = findextrema(Amem_sol; min_height=0.1)
    # println("Number of FFT peaks: $(length(indx_max))")
    if length(indx_max) < 2 || length(indx_min) < 2
        #* if there are not enough peaks in the time domain, return just the fitness
        # println("Not enough peaks in the time domain")
        return 0.0, 0.0
    end

    return getPerAmp(solt, indx_max, vals_max, indx_min, vals_min)
end

"""Calculates the period and amplitude of each individual in the population."""
function getPerAmp(solt, indx_max::Vector{Int}, vals_max::Vector{T}, indx_min::Vector{Int}, vals_min::Vector{T}) where T<:Real

    #* Calculate amplitudes and periods
    pers = (solt[indx_max[i+1]] - solt[indx_max[i]] for i in 1:(length(indx_max)-1))
    amps = (vals_max[i] - vals_min[i] for i in 1:min(length(indx_max), length(indx_min)))

    return mean(pers), mean(amps) 
end
#> END OF PERIOD AND AMPLITUDE FUNCTIONS ##




#< OSCILLATION DETECTION HEURISTICS ##
"""
    is_steadystate(solu::Vector{Float64}, solt::Vector{Float64})

Checks if the last tenth of the solution array is steady state
"""
function is_steadystate(solu::Vector{Float64}, solt::Vector{Float64})
    tstart = cld(length(solt),10) 

    #* Check if last tenth of the solution array is steady state
    testwindow = solu[end-tstart:end]
    if std(testwindow; mean=mean(testwindow)) < 0.01  
        return true
    else
        return false
    end 
end

function get_std_last10th(solu::Vector{Float64}, solt::Vector{Float64})
    tstart = cld(length(solt),10) 

    #* Test window of last 10% of solution
    testwindow = solu[end-tstart:end]
    
    return std(testwindow; mean=mean(testwindow))
end

"""
    is_oscillatory(solu::Vector{Float64}, solt::Vector{Float64}, max_idxs::Vector{Int}, min_idxs::Vector{Int})

Checks if the solution is oscillatory by checking if there are more than 1 maxima and minima in the solution array and whether the solution is steady state
"""
function is_oscillatory(solu::Vector{Float64}, solt::Vector{Float64}, max_idxs::Vector{Int}, min_idxs::Vector{Int})
    if !is_steadystate(solu, solt) && length(max_idxs) > 1 && length(min_idxs) > 1
        return true
    else
        return false
    end
end
#> END OF OSCILLATION DETECTION HEURISTICS ##





















