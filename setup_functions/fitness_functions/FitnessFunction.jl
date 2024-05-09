function calculate_fitness(sol::ODESolution)

    phenotype = zeros(3)

    calculate_fitness!(phenotype, sol)
    return phenotype
end



function calculate_fitness!(phenotype, sol::ODESolution)

    Amem = get_Amem(sol)

    #* Get the rfft of the solution and normalize it
    fftData = getFrequencies(Amem)

    #* get the indexes of the peaks in the fft
    fft_peakindexes = findmaxpeaks(fftData) 

    #* if there is no signal in the frequency domain, return 0.0s
    if length(fft_peakindexes) < 2 
        return phenotype
    end

    #* get the average standard deviation of the window around each peak in frequency domain
    standard_deviation = getSTD(fft_peakindexes, fftData)/length(fft_peakindexes)

    #* get the summed differences between all the peaks in frequency domain
    fft_peakvals = @view fftData[fft_peakindexes]
    sum_diff = getDif(fft_peakvals) 

    #* fitness is the sum of the standard deviation and the difference between the first and last peaks
    phenotype[1] = standard_deviation + sum_diff


    start_idx = cld(length(Amem), 10)
    indx_max, indx_min = findextrema(@view Amem[start_idx:end]; min_height=0.1) 
    if length(indx_max) < 2 || length(indx_min) < 2
        #* if there are not enough peaks in the time domain, return just the fitness
        return phenotype
    end

    vals_max = @view Amem[indx_max]
    vals_min = @view Amem[indx_min]

    #* get the period and amplitude of the solution
    phenotype[2:3] .= getPerAmp(sol.t, indx_max, vals_max, vals_min)
    return nothing
end

function calculate_fitness!(phenotype, saved_values::SavedValues)

    Amem = saved_values.saveval

    #* Get the rfft of the solution and normalize it
    fftData = getFrequencies(Amem)

    #* get the indexes of the peaks in the fft
    fft_peakindexes = findmaxpeaks(fftData) 

    #* if there is no signal in the frequency domain, return 0.0s
    if length(fft_peakindexes) < 2 
        return phenotype
    end

    #* get the average standard deviation of the window around each peak in frequency domain
    standard_deviation = getSTD(fft_peakindexes, fftData)/length(fft_peakindexes)

    #* get the summed differences between all the peaks in frequency domain
    fft_peakvals = @view fftData[fft_peakindexes]
    sum_diff = getDif(fft_peakvals) 

    #* fitness is the sum of the standard deviation and the difference between the first and last peaks
    phenotype[1] = standard_deviation + sum_diff


    start_idx = cld(length(Amem), 10)
    indx_max, indx_min = findextrema(@view Amem[start_idx:end]; min_height=0.1) 
    if length(indx_max) < 2 || length(indx_min) < 2
        #* if there are not enough peaks in the time domain, return just the fitness
        return phenotype
    end

    vals_max = @view Amem[indx_max]
    vals_min = @view Amem[indx_min]

    #* get the period and amplitude of the solution
    phenotype[2:3] .= getPerAmp(saved_values.t, indx_max, vals_max, vals_min)
    return nothing
end




function calculate_fitness!(phenotype, sol::ODESolution, Amem, fftData)


    Amem .= get_Amem(sol)

    #* Get the rfft of the solution and normalize it
    getFrequencies!(fftData, Amem)

    #* get the indexes of the peaks in the fft
    fft_peakindexes = findmaxpeaks(fftData) 

    #* if there is no signal in the frequency domain, return 0.0s
    if length(fft_peakindexes) < 2 
        return phenotype
    end

    #* get the average standard deviation of the window around each peak in frequency domain
    standard_deviation = getSTD(fft_peakindexes, fftData) 

    #* get the summed differences between all the peaks in frequency domain
    fft_peakvals = @view fftData[fft_peakindexes]
    sum_diff = getDif(fft_peakvals) 

    #* fitness is the sum of the standard deviation and the difference between the first and last peaks
    phenotype[1] = standard_deviation + sum_diff


    start_idx = cld(length(Amem), 10)
    indx_max, indx_min = findextrema(@view Amem[start_idx:end]; min_height=0.1) 
    if length(indx_max) < 2 || length(indx_min) < 2
        #* if there are not enough peaks in the time domain, return just the fitness
        return phenotype
    end

    vals_max = @view Amem[indx_max]
    vals_min = @view Amem[indx_min]

    #* get the period and amplitude of the solution
    phenotype[2:3] .= getPerAmp(sol.t, indx_max, vals_max, vals_min)
    return phenotype
end





























