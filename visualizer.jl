#< IMPORTS >#
begin 
    using DataFrames
    using CSV
    using Catalyst
    using ModelingToolkit
    using OrdinaryDiffEq
    using FFTW
    using StatsBase
end


include("full_model.jl")
include("ode_solving_functions.jl")
include("custom_peakfinder.jl")
include("fitness_function_helpers.jl")
#> IMPORTS <#


# LOAD DATA
df = CSV.read("smalldata.csv", DataFrame)



using GLMakie
theme = merge(theme_ggplot2(), theme_latexfonts())
set_theme!(theme)

function make_ode_solver(prob::ODEProblem)

    osys = prob.f.sys

    #* Get vector of all tunable parameters and states
    tunable_variables = get_tunable_variables(osys)

    num_tunable_parameters = length(tunable_parameters(osys; default = true))

    saved_idxs = get_Amem_indices(osys)

    function solver(input)::ODESolution

        #* Create a mapping the tunable variables to their values
        tunable_var_value_map = [var => val for (var, val) in zip(tunable_variables, input)]

        #* Split the input into parameters and initial conditions
        params = @view tunable_var_value_map[1:num_tunable_parameters]
        u0 = @view tunable_var_value_map[num_tunable_parameters+1:end]        

        #* Remake the problem with the new parameters
        newprob = remake(prob; p = params, u0 = u0)

        #* Solve the problem and return the solution
        solve(newprob, Rosenbrock23(); saveat = 0.1, save_idxs = saved_idxs, verbose=false, maxiters=1e6)
    end
    return solver
end



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

#- Plots interactive Makie timeseries plot with sliders for each parameter
function plot_interactive_parameters_timeseries(df::AbstractDataFrame)
    #* Make ODEProblem from ReactionSystem
    odeprob = make_ODEProblem((;));

    #* Make solver function from ODEProblem
    ode_solver = make_ode_solver(odeprob);


    fig = Figure(; size = (1200, 900));
    time_ax = Axis(fig[1, 1:3], title = "ODE Solution",
                xlabel = "Time (s)", xlabelsize = 18, 
                ylabel = "AP2 Membrane Localization", ylabelsize = 12, 
                limits = ((0.0, 2000.0), (0.0, 1.0)))

    fft_ax = Axis(fig[1, 4], title = "FFT",
                xlabel = "Period (s)", xlabelsize = 18, xtickformat = values -> ["$(round(Int, cld(1, value)))" for value in values], xscale = log10,
                ylabel = "Power", ylabelsize = 18,
                limits = ((0.0001, 1.0), nothing))


    #- Menu to select different regime groupings of the DataFrame
    subset_functions = @views [
        df[(df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :],
        df[(df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :],
        df[(df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :],
        df[(df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :]
    ]
    option_labels = ["K > P, Kmᴸᴷ > Kmᴸᴾ", "K > P, Kmᴸᴷ < Kmᴸᴾ", "K < P, Kmᴸᴷ > Kmᴸᴾ", "K < P, Kmᴸᴷ < Kmᴸᴾ"]
    fig[3, 3:4] = menu_grid = GridLayout()
    menu_label = Label(menu_grid[1, 1:2], "Regime", fontsize = 25, color = :black, valign = :top)
    menu = Menu(menu_grid[2, 1:2], options = zip(option_labels, subset_functions), default = "K > P, Kmᴸᴷ > Kmᴸᴾ", valign = :center)

    #* Technically just a subsetting function, because making an Observable DataFrame behaves weirdly
    # df_subsetter = Observable{Any}(subset_functions[1])
    data = Observable{Any}(subset_functions[1])

    on(menu.selection) do selection
        data[] = selection
    end
    notify(menu.selection)



    
    #- Row slider
    #* Make row slider to cycle through the dataframe
    row_slider = Slider(fig[5, 1:4], range = 1:1:1000, startvalue = 1, horizontal = true, color_inactive = :pink, color_active = :red, color_active_dimmed = :pink)

    row_slider_label = lift(row_slider.value) do rownumber
        # n_data = nrow(subsetter(df))
        return "Row $(rownumber)/1000"
    end 

    #* Label for row slider
    Label(fig[6, 1:4], row_slider_label, fontsize = 20, color = :black)

    #* Make observable for the row number, returns vector of input parameters
    row = lift(row_slider.value) do rownumber
        # data = subsetter(df)
        if isempty(df)
            return zeros(17)
        else
            return collect(Float64, df[rownumber, Between(:kfᴸᴬ, :A)])
        end
    end


    #- Reset button that resets the row to what the slider is at
    reset_button = Button(fig[4, 1], label = "Reset", labelcolor = :red)

    on(reset_button.clicks) do n 
        notify(row_slider.value)
    end


    slider_labels = names(df[!, Between(:kfᴸᴬ, :A)])
    parameter_names = slider_labels[1:13]
    species_names = slider_labels[14:end]
    lb, ub = get_tunable_bounds(odeprob.f.sys)
    param_lb, param_ub = lb[1:13], ub[1:13]

    #- Parameter value sliders
    # fig[3, 1:2] = parameter_slider_grid = GridLayout()

    # Label(parameter_slider_grid[1, 1:2], "Parameters", fontsize = 25)
    Label(fig[2, 1:2], "Parameters", fontsize = 25)
    #* Make slider grid with slider for each parameter in ind 
    parameter_sliders = SliderGrid(fig[3,1:2],
                                ((label = label, range = range(param_lb[i], param_ub[i], 100000), startvalue = df[1, label]) for (i, label) in enumerate(parameter_names))...)

    # Adjust DF range 
    parameter_sliders.sliders[13].range = range(0.0, 10000.0, 100000)


    #- Have parameter sliders listen to row
    on(row) do row
        set_close_to!.(parameter_sliders.sliders, row[1:13])
    end

    #* Makes the vector of observables (DIFFERENT FROM OBSERVABLE VECTOR)
    parameter_slider_observables = [s.value for s in parameter_sliders.sliders]

    #- Initial conditions sliders
    # fig[3, 3:4] = species_slider_grid = GridLayout()

    ic_lb, ic_ub = lb[14:end], ub[14:end]
    # Label(species_slider_grid[1, 1:2], "Initial Conditions", fontsize = 25)
    Label(fig[2, 3:4], "Initial Conditions", fontsize = 25)
    species_sliders = SliderGrid(fig[3,3:4],
                                ((label = label, range = range(ic_lb[i], ic_ub[i], 100000), startvalue = df[1, label]) for (i, label) in enumerate(species_names))...; valign = :top)

    on(row) do row
        set_close_to!.(species_sliders.sliders, row[14:17])
    end

    species_slider_observables = [s.value for s in species_sliders.sliders]


    all_slider_observables = vcat(parameter_slider_observables, species_slider_observables)

    #* Make observable vector that is lifted from the slider_observables
    observable_vector = lift(all_slider_observables...) do slvalues...
        [slvalues...]
    end

    function compute_Km_vals(observables)
        Kmᴸᴷ = round((observables[4] + observables[5])/observables[3]; digits = 2)
        Kmᴸᴾ = round((observables[7] + observables[8])/observables[6]; digits = 2)
        return (Kmᴸᴷ, Kmᴸᴾ)
    end

    Km_vals = lift(observable_vector) do observables
        return compute_Km_vals(observables)
    end

    # Function to determine the regime based on parameter values
    function get_regime(input_vector)
        K_P_label = input_vector[15] > input_vector[16] ? "K > P" : "K < P"
        Kmᴸᴷ, Kmᴸᴾ = Km_vals.val
        Kmᴸᴷ_label = Kmᴸᴷ > Kmᴸᴾ ? "Kmᴸᴷ > Kmᴸᴾ" : "Kmᴸᴷ < Kmᴸᴾ"
        return (K_P_label, Kmᴸᴷ_label)
    end

    # Define colors for each regime combination
    regime_color_dict = Dict(
        ("K > P", "Kmᴸᴷ > Kmᴸᴾ") => :blue,
        ("K > P", "Kmᴸᴷ < Kmᴸᴾ") => :red,
        ("K < P", "Kmᴸᴷ > Kmᴸᴾ") => :green,
        ("K < P", "Kmᴸᴷ < Kmᴸᴾ") => :purple
    )

    regime = lift(observable_vector) do observables
        get_regime(observables)
    end

    regime_color = lift(regime) do regime
        return regime_color_dict[regime]
    end

    # regime_label = lift(regime) do regime
    #     return join(regime, "\n")
    # end

    Km_vals_label = lift(Km_vals) do Km_vals
        return join(["Kmᴸᴷ: $(Km_vals[1])", "Kmᴸᴾ: $(Km_vals[2])"], "\n")
    end


    # Label(fig[4, 3], regime_label, fontsize = 25, color = regime_color, valign = :top)
    Label(fig[4, 4], Km_vals_label, fontsize = 25, color = regime_color, valign = :top)




    #* Update function that takes the observables vector and resolves the ODE
    function update_function(params)
        sol = ode_solver(params)
        return compute_Amem(sol)
    end

    testsol = ode_solver(collect(df[1, Between(:kfᴸᴬ, :A)]))
    frequencies = [i*10/length(testsol) for i in 1:(length(testsol)/2)+1]
    # periods = 1 ./ frequencies


    #* Make observable for AP2 Membrane Localization from calling the update function on the observable vector
    Amem = lift(observable_vector) do observables
        update_function(observables)
    end


    function getFrequencies(timeseries::Vector{Float64})
        rfft_result = rfft(timeseries)
        norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
        abs.(rfft_result) ./ norm_val
    end

    fft_Amem = lift(Amem) do Amem
        getFrequencies(Amem)
    end

    diff_and_std_label = lift(fft_Amem) do fft_Amem
        difs, stds = round.(get_fitness(fft_Amem);digits = 5)
        fitness = round(difs + stds; digits = 5)
        return "Diff: $(difs)\nSTD: $(stds)\nFitness: $fitness"
    end
    Label(fig[4, 3], diff_and_std_label, fontsize = 25, color = :black, valign = :top)

    ln = lines!(time_ax, testsol.t, Amem, color = regime_color, linewidth = 3)
    lines!(fft_ax, frequencies, fft_Amem, color = regime_color, linewidth = 3)
    fig
end

plot_interactive_parameters_timeseries(df)



