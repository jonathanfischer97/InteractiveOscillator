#< IMPORTS >#
begin 
    using DataFrames
    using CSV
    using Catalyst
    using ModelingToolkit
    using OrdinaryDiffEq
    # using FFTW
    using StatsBase
    using GLMakie
    using SymbolicIndexingInterface
    GLMakie.activate!(; title = "Interactive Oscillator Visualizer", float = true, focus_on_show = true)
    theme = merge(theme_ggplot2(), theme_latexfonts())
    set_theme!(theme)
end

#* Set the script directory to the current directory
script_dir = @__DIR__
cd(script_dir)

# Assert we are in the correct directory
@assert basename(pwd()) == "InteractiveOscillator" 


# Load dependencies
begin
    include("setup_functions/full_model.jl")
    include("setup_functions/ode_solving_functions.jl")
    # include("setup_functions/fitness_functions/helpers/custom_peakfinder.jl")
    # include("setup_functions/fitness_functions/helpers/fitness_function_helpers.jl")
    # include("setup_functions/fitness_functions/FitnessFunction.jl")
    include("setup_functions/load_data.jl")
end
#> IMPORTS <#



#- Plots interactive Makie timeseries plot with sliders for each parameter
function plot_interactive_parameters_timeseries(df::AbstractDataFrame)

    println("Plotting...")
    fig = Figure(; size = (1200, 900));
    time_ax = Axis(fig[1, 1:5], title = "ODE Solution",
                xlabel = "Time (s)", xlabelsize = 18, 
                ylabel = "AP2 Membrane Localization", ylabelsize = 12, 
                limits = (ODEPROB.tspan, (0.0, 1.0)))

    # fft_ax = Axis(fig[1, 4:5], title = "FFT",
    #             xlabel = "Period (s)", xlabelsize = 18, xtickformat = values -> ["$(round(Int, cld(1, value)))" for value in values], xscale = log10,
    #             ylabel = "Power", ylabelsize = 18,
    #             limits = ((0.0001, 1.0), nothing))

    #< GRIDS ##
    fig[2, 1:3] = parameter_slider_grid = GridLayout()
    fig[2, 4:5] = species_slider_grid = GridLayout()
    species_slider_grid[3, 1:2] = menu_grid = GridLayout()
    species_slider_grid[4, 1:2] = misc_grid =  GridLayout()


    #- Menu to select different regime groupings of the DataFrame
    # subset_views = [
    #     view(df, (df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :),
    #     view(df, (df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :),
    #     view(df, (df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ), :),
    #     view(df, (df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ), :)
    # ]

    # function get_subset_indices(df)
    #     subsets = [
    #         (df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
    #         (df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ),
    #         (df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
    #         (df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ)
    #     ]
    #     return [findall(subset) for subset in subsets]
    # end
    
    # subset_indices = get_subset_indices(df)

    # function get_matrix_subsets(indices, df)
    #     df_matrix = permutedims(Matrix(df[!, Between(:kfᴸᴬ, :A)]))
    #     return [df_matrix[:, idx] for idx in indices]
    # end
    # subsetted_matrices = get_matrix_subsets(subset_indices,df)


    function get_matrix_subsets(df)
        subsets = [
            (df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
            (df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ),
            (df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
            (df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ)
        ]
        return [permutedims(Matrix(df[subset, Between(:kfᴸᴬ, :APLp)])) for subset in subsets]
    end
    subsetted_matrices = get_matrix_subsets(df)



    option_labels = ["K > P, Kmᴸᴷ > Kmᴸᴾ", "K > P, Kmᴸᴷ < Kmᴸᴾ", "K < P, Kmᴸᴷ > Kmᴸᴾ", "K < P, Kmᴸᴷ < Kmᴸᴾ"]
    menu_label = Label(menu_grid[1, 1:2], "Regime Selector", fontsize = 25, color = :black, valign = :top, tellheight = true)
    menu = Menu(menu_grid[2, 1:2], options = zip(option_labels, subsetted_matrices), default = "K > P, Kmᴸᴷ > Kmᴸᴾ", valign = :top, tellheight = false)

    #* Make observable for the data that is lifted from the menu selection
    data_observable = Observable{Any}(subsetted_matrices[1])

    # data = lift(menu.selection) do selection
    #     return selection
    # end
    
    #- Row slider
    #* Make row slider to cycle through the dataframe
    row_slider = Slider(fig[3, 2:end], range = 1:1:1000, startvalue = 1, horizontal = true, color_inactive = :pink, color_active = :red, color_active_dimmed = :pink, valign = :center, tellheight = true)

    row_slider_label = lift(row_slider.value) do rownumber
        return "Row $(rownumber)/1000"
    end 

    #* Label for row slider
    Label(fig[3, 1], row_slider_label, fontsize = 20, color = :black, valign = :top, tellheight = true)

    #* Make observable for the row number, returns vector of input parameters
    # row = lift(row_slider.value) do rownumber
    #     # if nrow(data) == 0
    #     #     return zeros(17)
    #     # else
    #     #     return collect(Float64, data[rownumber, Between(:kfᴸᴬ, :A)])
    #     # end
    #     subsetted_matrix = data[]
    #     return subsetted_matrix[:, rownumber]
    # end


    #- Reset button that resets the row to what the slider is at
    reset_button = Button(misc_grid[1, 2], label = "Reset", labelcolor = :red)

    on(reset_button.clicks) do n 
        notify(row_slider.value)
    end


    on(menu.selection) do selection
        data_observable[] = selection
        notify(row_slider.value)
    end
    notify(menu.selection)


    function logrange(start, stop, steps)
        startval = start == 0.0 ? 1.0 : log10(start)
        return 10 .^ range(startval, log10(stop), length=steps)
    end

    slider_labels = names(df[!, Between(:kfᴸᴬ, :APLp)])
    parameter_names = @view slider_labels[1:13]
    species_names = @view slider_labels[14:end]
    # lb, ub = get_tunable_bounds(odeprob.f.sys)
    # param_lb, param_ub = @views lb[1:13], ub[1:13]

    #- Parameter value sliders
    Label(parameter_slider_grid[1, 1:3], "Parameters", fontsize = 25)
    #* Make slider grid with slider for each parameter in ind 
    parameter_sliders = SliderGrid(parameter_slider_grid[2, 1:3],
                                ((label = label, range = logrange(INPUT_BOUNDS[Symbol(label)]..., 100000), startvalue = df[1, Symbol(label)]) for (i, label) in enumerate(parameter_names))...)

    # Adjust DF range 
    parameter_sliders.sliders[13].range = range(0.0, 10000.0, 100000)


    #- Have parameter sliders listen to row
    # on(row) do row
    #     new_parameters = view(row, 1:13)
    #     set_close_to!.(parameter_sliders.sliders, new_parameters)
    # end

    on(row_slider.value) do rownumber
        subsetted_matrix = data_observable[]
        new_parameters = view(subsetted_matrix, 1:13, rownumber)

        # new_parameters = view(row, 1:13)
        set_close_to!.(parameter_sliders.sliders, new_parameters)
    end

    #* Makes the vector of observables (DIFFERENT FROM OBSERVABLE VECTOR)
    parameter_slider_observables = [s.value for s in parameter_sliders.sliders]

    #- Initial conditions sliders
    # ic_lb, ic_ub = @views lb[14:end], ub[14:end]
    Label(species_slider_grid[1, 1:2], "Initial Conditions", fontsize = 25)
    species_sliders = SliderGrid(species_slider_grid[2, 1:2],
                                ((label = label, range = logrange(INPUT_BOUNDS[Symbol(label)]..., 100000), startvalue = df[1, Symbol(label)]) for (i, label) in enumerate(species_names))...; valign = :top, tellheight = false)

    # on(row) do row
    #     new_species = view(row, 14:17)
    #     set_close_to!.(species_sliders.sliders, new_species)
    # end

    on(row_slider.value) do rownumber
        subsetted_matrix = data_observable[]
        new_species = view(subsetted_matrix, 14:17, rownumber)

        # new_parameters = view(row, 1:13)
        set_close_to!.(species_sliders.sliders, new_species)
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
    Label(misc_grid[1, 1], Km_vals_label, fontsize = 25, color = regime_color, valign = :center, tellwidth = false)




    # #* Update function that takes the observables vector and resolves the ODE
    # function update_ode_solution(params)
    #     solve_odes(params)
    #     # return compute_Amem(sol)
    # end

    # testsol = solve_odes(collect(df[1, Between(:kfᴸᴬ, :APLp)]))
    # frequencies = [i*10/length(testsol) for i in 1:(length(testsol)/2)+1]
    # periods = 1 ./ frequencies


    #* Make observable for AP2 Membrane Localization from calling the update function on the observable vector
    Amem = lift(observable_vector) do observables
        get_Amem(solve_odes(observables))
    end

    # ode_sol = lift(observable_vector) do observables
    #     ode_solver(observables)
    # end

    # function getFrequencies(timeseries::Vector{Float64})
    #     rfft_result = RFFT_PLAN * timeseries
    #     norm_val = length(timeseries)/ 2 #* normalize by length of timeseries
    #     abs.(rfft_result) ./ norm_val
    # end

    # fft_Amem = lift(Amem) do Amem
    #     getFrequencies(Amem)
    # end

    # diff_and_std_label = lift(fft_Amem) do fft_Amem
    #     difs, stds = round.(get_fitness(fft_Amem);digits = 5)
    #     fitness = round(difs + stds; digits = 5)
    #     return "Diff: $(difs)\nSTD: $(stds)\nFitness: $fitness"
    # end
    # Label(fig[1, 4:5], diff_and_std_label, fontsize = 20, color = :black, valign = :top, tellwidth = false, tellheight=false)

    ln = lines!(time_ax, range(0.0, ODEPROB.tspan[2], length(Amem.val)), Amem, color = regime_color, linewidth = 3)
    # lines!(fft_ax, frequencies, fft_Amem, color = regime_color, linewidth = 3)
    fig
end

data = load_data()
data.P

names(data)
fig = plot_interactive_parameters_timeseries(load_data())
display(fig)

println("Press enter to exit")
readline()

ODEPROB.tspan

