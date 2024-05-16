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
    using DynamicalSystems

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






diffeq = (alg = Rosenbrock23(autodiff= false),  adaptive = false, dt = 0.01, abstol = 1e-10)
dsprob = CoupledODEs(ODEPROB, diffeq)
osys_full = OSYS


# propertynames(osys_full)
# fieldnames(typeof(osys_full))




# using Symbolics
# Define which parameters will be interactive during the simulation
parameter_sliders = Dict(
    # can use integer indexing
    # 1 => 0:10:1000,
    # the global scope symbol
    :DF => 100:10:10000,
)

# Define what variables will be visualized as timeseries
observables = [
    # 1,         # can use integer indexing
    # z,         # MTK state variable (called "unknown")
    osys_full.Amem, # MTK observed variable
    # osys_full.LpAPLp 
    # :y,        # `Symbol` instance with same name as symbolic variable
    # power,     # arbitrary function of the state
    # x^2 - y^2, # arbitrary symbolic expression of symbolic variables
]

# Define what variables will be visualized as state space trajectory
# same as above, any indexing works, but ensure to make the vector `Any`
# so that integers are not converted to symbolic variables
idxs = Any[osys_full.Amem, osys_full.LpAPLp, osys_full.LpA]

u0s = [
    # we can specify dictionaries, each mapping the variable to its value
    # un-specified variables get the value they currently have in `ds`
    Dict(),
    Dict(:L => DEFAULT_TUPLE.Lp),
    # Dict(:L => DEFAULT_TUPLE.A, :A => DEFAULT_TUPLE.L)
    Dict(:K => DEFAULT_TUPLE.P, :P => DEFAULT_TUPLE.K)
]

update_theme!(fontsize = 14)
tail = 2000

fig, dsobs = interactive_trajectory_timeseries(dsprob, observables, u0s;
    parameter_sliders, Δt = 0.01, tail, idxs,
    figure = (size = (1100, 650),)
)

step!(dsobs, 2tail)

display(fig)

fig

