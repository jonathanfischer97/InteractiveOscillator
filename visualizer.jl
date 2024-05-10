#< IMPORTS >#
begin 
    using DataFrames
    using CSV
    using Catalyst
    using ModelingToolkit
    using OrdinaryDiffEq
    using StatsBase
    using GLMakie
    using SymbolicIndexingInterface
    # using DimensionalData
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
    include("setup_functions/data_handling.jl")
end
#> IMPORTS <#


function logrange(start, stop, steps)
    startval = start == 0.0 ? 1.0 : log10(start)
    logrange_vec = 10 .^ range(startval, log10(stop), length=steps)
    if start == 0.0
        logrange_vec[1] = 0.0
    end
    return logrange_vec
end




#- Plots interactive Makie timeseries plot with sliders for each parameter
function plot_interactive_parameters_timeseries(df::AbstractDataFrame)

    println("Plotting...")
    fig = Figure(; size = (1200, 900));
    time_ax = Axis(fig[1, 1:5], title = "ODE Solution",
                xlabel = "Time (s)", xlabelsize = 18, 
                ylabel = "Percentage of Total AP2 %", ylabelsize = 14, yscale = Makie.pseudolog10,
                limits = (ODEPROB.tspan, (0.0, 100.0)))


    #< GRIDS ##
    fig[2, 1:3] = parameter_slider_grid = GridLayout()
    fig[2, 4:5] = species_slider_grid = GridLayout()
    species_slider_grid[3, 1:2] = menu_grid = GridLayout()
    species_slider_grid[4, 1:2] = misc_grid =  GridLayout()


    #- Menu to select different regime groupings of the DataFrame

    # data_dimarray = convert_to_dimarray(df)

    function get_matrix_subsets(df)
        subsets = [
            (df.K .> df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
            (df.K .> df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ),
            (df.K .< df.P) .& (df.Kmᴸᴷ .> df.Kmᴸᴾ),
            (df.K .< df.P) .& (df.Kmᴸᴷ .< df.Kmᴸᴾ),
            df.A .> df.L,
            df.A .< df.L,
        ]
        return [permutedims(Matrix(@view df[subset, DataFrames.Between(:kfᴸᴬ, :Lp)])) for subset in subsets]
    end
    subsetted_matrices = get_matrix_subsets(df)



    option_labels = ["K > P, Kmᴸᴷ > Kmᴸᴾ", "K > P, Kmᴸᴷ < Kmᴸᴾ", "K < P, Kmᴸᴷ > Kmᴸᴾ", "K < P, Kmᴸᴷ < Kmᴸᴾ", "A > L", "A < L"]
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


    slider_labels = names(df[!, DataFrames.Between(:kfᴸᴬ, :Lp)])
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
    # parameter_sliders.sliders[13].range = range(0.0, 10000.0, 100000)


    on(row_slider.value) do rownumber
        subsetted_matrix = data_observable[]
        new_parameters = view(subsetted_matrix, 1:13, rownumber)

        # new_parameters = view(row, 1:13)
        set_close_to!.(parameter_sliders.sliders, new_parameters)
    end

    #* Makes the vector of observables (DIFFERENT FROM OBSERVABLE VECTOR)
    parameter_slider_observables = [s.value for s in parameter_sliders.sliders]

    #- Initial conditions sliders
    Label(species_slider_grid[1, 1:2], "Initial Conditions", fontsize = 25)
    species_sliders = SliderGrid(species_slider_grid[2, 1:2],
                                ((label = label, range = logrange(INPUT_BOUNDS[Symbol(label)]..., 100000), startvalue = df[1, Symbol(label)]) for (i, label) in enumerate(species_names))...; valign = :top, tellheight = false)


    on(row_slider.value) do rownumber
        subsetted_matrix = data_observable[]
        new_species = view(subsetted_matrix, 14:18, rownumber)

        set_close_to!.(species_sliders.sliders, new_species)
    end

    species_slider_observables = [s.value for s in species_sliders.sliders]

    all_slider_observables = vcat(parameter_slider_observables, species_slider_observables)


    tmp_vector = copy(FULL_INDIVIDUAL)

    #* Make observable vector that is lifted from the slider_observables
    observable_vector = lift(all_slider_observables...) do slvalues...
        tmp_vector[1:18] .= [slvalues...]
        tmp_vector
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


    #* Make observable for AP2 Membrane Localization from calling the update function on the observable vector


    sol = lift(observable_vector) do observables
        solve_odes(observables)
    end

    Amem = lift(sol) do sol
        return get_Amem(sol) .* 100
    end

    # L = lift(sol) do sol
    #     return get_L(sol)
    # end

    A = lift(sol) do sol
        return (get_A(sol) ./ sol[4,1]) .* 100
    end

    LpA = lift(sol) do sol
        return (get_LpA(sol) ./ sol[4,1]) .* 100
    end

    LpAK = lift(sol) do sol
        return (get_LpAK(sol) ./ sol[4,1]) .* 100
    end

    LpAP = lift(sol) do sol
        return (get_LpAP(sol) ./ sol[4,1]) .* 100
    end

    LpAKL = lift(sol) do sol
        return (get_LpAKL(sol) ./ sol[4,1]) .* 100
    end

    LpAPLp = lift(sol) do sol
        return get_LpAPLp(sol) ./ sol[4,1] .* 100 
    end

    tspan = 0.0:0.1:ODEPROB.tspan[2]
    lines!(time_ax, tspan, Amem, color = :black, linewidth = 4, linestyle = :dash, label = "Amem")
    # lines!(time_ax, tspan, L, color = :blue, linewidth = 3)
    lines!(time_ax, tspan, A, color = :brown, linewidth = 3, label = "A", linestyle = :dash)
    lines!(time_ax, tspan, LpA, color = :red, linewidth = 3, label = "LpA")
    lines!(time_ax, tspan, LpAK, color = :blue, linewidth = 3, label = "LpAK")
    lines!(time_ax, tspan, LpAP, color = :orange, linewidth = 3, label = "LpAP")
    lines!(time_ax, tspan, LpAKL, color = :green, linewidth = 3, label = "LpAKL")
    lines!(time_ax, tspan, LpAPLp, color = :purple, linewidth = 3, label = "LpAPLp")
    axislegend(time_ax, position = :rt, labelsize = 20)

    # lines!(time_ax, sol, color = regime_color, linewidth = 3)

    fig
end

data = load_df()
data = load_df("/home/jfisch27/Desktop/ThesisStuff/GeometricallyTunableOscillator/data/alldata.csv"; islocal = false)
# data[data.DF .!= 50.0 .&& data.DF .!= 100.0 .&& data.DF .!= 500.0 .&& data.DF .!= 1000.0, :]
fig = plot_interactive_parameters_timeseries(data)

display(fig)

println("Press enter to exit")
readline()





p = [0.639, 0.376, 0.626, 87.8, 614.0, 43.0, 4.02, 6.63, 67.4, 0.501, 0.278, 0.001, 27.0]

u0 = [:L => 15.4, :K => 0.44, :P => 1.98, :A => 1.48]

newprob = remake(ODEPROB, p = p, u0 = u0)

sol = solve(newprob, Rosenbrock23(), saveat = 0.1)




data = CSV.read("/home/jfisch27/Desktop/ThesisStuff/GeometricallyTunableOscillator/data/alldata.csv", DataFrame)

scatter(data.Kmᴸᴷ./data.Kmᴸᴾ, data.K./data.P, legend = false, axis = (;xscale = log10, yscale = log10))


scatter(data.DF, data.Kdᴸᴬ, legend = false, axis = (;xscale = log10, yscale = identity))

colormap = Vector(data.Amplitude)
scatter(data.K ./ data.P, data.Kmᴸᴷ./data.Kmᴸᴾ, color = colormap, legend = false, axis = (;xscale = log10, yscale = log10, xlabel = "K/P", ylabel = "Kmᴸᴷ/Kmᴸᴾ"))

