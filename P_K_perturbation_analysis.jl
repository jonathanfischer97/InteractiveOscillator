#< IMPORTS >#
using DrWatson
@quickactivate "GeometricallyTunableOscillator"
# include(srcdir("OscTools", "OscTools.jl"))
include("/home/jfisch27/Desktop/ThesisStuff/GeometricallyTunableOscillator/src/OscTools/OscTools.jl")

using .OscTools

begin 
    using DataFrames
    using CSV
    using StatsBase
    using CategoricalArrays
    using DataFramesMeta
    using ModelingToolkit
    using OrdinaryDiffEq
end

include("/home/jfisch27/Desktop/ThesisStuff/GeometricallyTunableOscillator/src/Kd_analysis_functions.jl")

using ColorSchemes
#> IMPORTS <#


#* Load and process the data
function load_data()
    path1 = datadir("full_df.csv")
    df1 = load_and_process_genotype_data(path1)

    path2 = datadir("KmLK_lt_KmLP_results.csv")
    df2 = CSV.read(path2, DataFrame)[!, Between(:per, :A)] |> process_genotype_data

    # path2 = datadir("ROCKFISH_DATA", "3Fixed", "PopSize_15000", "L_K_P", "combined_df.csv")
    # df2 = load_and_process_genotype_data(path2)
    # path3 = datadir("ROCKFISH_DATA", "3Fixed", "PopSize_15000", "L_K_A", "combined_df.csv")
    # df3 = load_and_process_genotype_data(path3)
    
    return vcat(df1, df2)
end

df = load_data()

df.A_L_ratio = df.A ./ df.L
df.Km_ratio = df.Kmᴸᴷ ./ df.Kmᴸᴾ

df.DF .= unwrap.(df.DF)





using GLMakie
theme = merge(theme_ggplot2(), theme_latexfonts())
set_theme!(theme)

function make_ode_solver(prob::ODEProblem)

    osys = prob.f.sys

    #* Get vector of all tunable parameters and states
    tunable_variables = get_tunable_variables(osys)

    num_tunable_parameters = length(tunable_parameters(osys; default = true))

    saved_idxs = OscTools.get_Amem_indices(osys)

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

function fitness_function(fftData)

    sliced_fftData = @view fftData[1:cld(length(fftData), 2)]

    #* get the indexes of the peaks in the fft
    fft_peakindexes, fft_peakvals = OscTools.findmaxpeaks(sliced_fftData) 

    #* if there is no signal in the frequency domain, return 0.0s
    if isempty(fft_peakvals)
        return [0.0, 0.0]
    end

    #* get the summed standard deviation of the peaks in frequency domain
    standard_deviation = OscTools.getSTD(fft_peakindexes, sliced_fftData) 

    #* get the summed difference between the first and last peaks in frequency domain
    sum_diff = OscTools.getDif(fft_peakvals) 
    # sum_diff = maximum(fft_peakvals)

    #* fitness is the sum of the standard deviation and the difference between the first and last peaks
    return [sum_diff, standard_deviation] .* 1e4
end


#- Plots interactive Makie timeseries plot with sliders for each parameter
function plot_interactive_parameters_timeseries(df::AbstractDataFrame)
    #* Make ODEProblem from ReactionSystem
    odeprob = make_ODEProblem((;));

    #* Make solver function from ODEProblem
    ode_solver = make_ode_solver(odeprob);


    fig = Figure(; size = (1400, 1000));
    time_ax = Axis(fig[1, 1:2], 
                xlabel = "Time (s)", xlabelsize = 20, 
                ylabel = "AP2 Membrane Localization", ylabelsize = 20, 
                limits = ((0, 2000), (0, 1)))

    fft_ax = Axis(fig[1, 3:4], 
                xlabel = "Frequency (Hz)", xlabelsize = 20, 
                ylabel = "Power", ylabelsize = 20,
                limits = ((2, 500), nothing))


    #- Menu to select different regime groupings of the DataFrame
    subset_functions = [
        (df) -> df,
        (df) -> df[(df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :],
        (df) -> df[(df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :],
        (df) -> df[(df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :],
        (df) -> df[(df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :],
        (df) -> df[(df.A .> df.L), :],
        (df) -> df[(df.A .< df.L), :]
    ]
    option_labels = ["All", "K > P, Kmᴸᴷ > Kmᴸᴾ", "K > P, Kmᴸᴷ < Kmᴸᴾ", "K < P, Kmᴸᴷ > Kmᴸᴾ", "K < P, Kmᴸᴷ < Kmᴸᴾ", "A > L", "A < L"]
    fig[3, 3:4] = menu_grid = GridLayout()
    menu_label = Label(menu_grid[1, 1:2], "Regime", fontsize = 25, color = :black, valign = :top)
    menu = Menu(menu_grid[2, 1:2], options = zip(option_labels, subset_functions), default = "All", valign = :center)

    #* Technically just a subsetting function, because making an Observable DataFrame behaves weirdly
    df_subsetter = Observable{Any}(subset_functions[1])

    on(menu.selection) do selection
        df_subsetter[] = selection
    end
    notify(menu.selection)



    
    #- Row slider
    #* Make row slider to cycle through the dataframe
    row_slider = Slider(fig[5, 1:4], range = 1:1:500, startvalue = 1, horizontal = true, color_inactive = :pink, color_active = :red, color_active_dimmed = :pink)

    row_slider_label = lift(row_slider.value, df_subsetter) do rownumber, subsetter
        n_data = nrow(subsetter(df))
        return "Row $(rownumber)/$(n_data)"
    end 

    #* Label for row slider
    Label(fig[6, 1:4], row_slider_label, fontsize = 20, color = :black)

    #* Make observable for the row number, returns vector of input parameters
    row = lift(row_slider.value, df_subsetter) do rownumber, subsetter
        data = subsetter(df)
        if isempty(data)
            return zeros(17)
        else
            return collect(Float64, data[rownumber, Between(:kfᴸᴬ, :A)])
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
    param_ub = ub[1:13]

    #- Parameter value sliders
    Label(fig[2, 1:2], "Parameters", fontsize = 25)
    #* Make slider grid with slider for each parameter in ind 
    parameter_sliders = SliderGrid(fig[3,1:2],
                                ((label = label, range = range(minimum(df[!,label]), param_ub[i], 100000), startvalue = df[1, label]) for (i, label) in enumerate(parameter_names))...)

    # Adjust DF range 
    parameter_sliders.sliders[13].range = range(0.0, 10000.0, 100000)


    #- Have parameter sliders listen to row
    on(row) do row
        set_close_to!.(parameter_sliders.sliders, row[1:13])
    end

    #* Makes the vector of observables (DIFFERENT FROM OBSERVABLE VECTOR)
    parameter_slider_observables = [s.value for s in parameter_sliders.sliders]

    #- Initial conditions sliders
    ic_ub = ub[14:end]
    Label(fig[2, 3:4], "Initial Conditions", fontsize = 25)
    species_sliders = SliderGrid(fig[3,3:4],
                                ((label = label, range = range(0.0, ic_ub[i], 100000), startvalue = df[1, label]) for (i, label) in enumerate(species_names))...; valign = :top)

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

    regime_label = lift(regime) do regime
        return join(regime, "\n")
    end

    Km_vals_label = lift(Km_vals) do Km_vals
        return join(["Kmᴸᴷ: $(Km_vals[1])", "Kmᴸᴾ: $(Km_vals[2])"], "\n")
    end


    Label(fig[4, 3], regime_label, fontsize = 25, color = regime_color, valign = :top)
    Label(fig[4, 4], Km_vals_label, fontsize = 25, color = regime_color, valign = :top)




    #* Update function that takes the observables vector and resolves the ODE
    function update_function(params)
        sol = ode_solver(params)
        return compute_Amem(sol)
    end


    #* Make observable for AP2 Membrane Localization from calling the update function on the observable vector
    Amem = lift(observable_vector) do observables
        update_function(observables)
    end

    fft_Amem = lift(Amem) do Amem
        return OscTools.getFrequencies(Amem; jump = 1)
    end

    diff_and_std_label = lift(fft_Amem) do fft_Amem
        difs, stds = round.(fitness_function(fft_Amem);digits = 2)
        fitness = round(difs + stds; digits = 2)
        return "Diff: $(difs)\nSTD: $(stds)\nFitness: $fitness"
    end
    Label(fig[4, 1], diff_and_std_label, fontsize = 25, color = :black, valign = :top)

    ln = lines!(time_ax, 0.0:0.1:2000.0, Amem, color = regime_color, linewidth = 3)
    lines!(fft_ax, fft_Amem, color = regime_color, linewidth = 3)
    fig
end

plot_interactive_parameters_timeseries(df[!, Not(:Period, :Amplitude)])






ind = [1.18902, 19.37129, 0.12065, 0.00216, 0.01814, 0.02333, 0.06029, 9.75757, 0.12002, 0.00624, 0.04642, 0.02592, 10000.0, 8.04727, 2.05306, 0.04264, 1.51217]

solver = make_ode_solver(make_ODEProblem((;)))

sol = solver(ind) |> compute_Amem

fftData = OscTools.getFrequencies(sol;jump = 1)

fitness = fitness_function(fftData[2:end])

fft_peakindexes, fft_peakvals = OscTools.findmaxpeaks(fftData) 

sum(abs.(diff(fft_peakvals)))/length(fft_peakvals)


OscTools.getDif(fft_peakvals)
OscTools.getSTD(fft_peakindexes, fftData)