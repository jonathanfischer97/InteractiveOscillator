function get_fitness(fftData)
    #* get the indexes of the peaks in the fft
    fft_peakindexes, fft_peakvals = findmaxpeaks(fftData) 

    #* if there is no signal in the frequency domain, return 0.0s
    if isempty(fft_peakvals)
        return [0.0, 0.0]
    end

    normalized_peakvals = fft_peakvals ./ maximum(fft_peakvals)
    sum_of_peakvals = sum(normalized_peakvals)

    #* get the summed standard deviation of the peaks in frequency domain
    standard_deviation = getSTD(fft_peakindexes, fftData) / sum_of_peakvals

    #* get the summed difference between the first and last peaks in frequency domain
    # sum_diff = OscTools.getDif(fft_peakvals) 
    sum_diff = maximum(fft_peakvals) / sum_of_peakvals

    #* fitness is the sum of the standard deviation and the difference between the first and last peaks
    return [sum_diff, standard_deviation] #.* 1e2
end