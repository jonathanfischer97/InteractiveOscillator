function solve_odes(individual, odeprob::ODEProblem = ODEPROB)
    newprob = remake_odeprob(individual, odeprob)
    # saved_values = SavedValues(Float64, Float64)
    # cb = SavingCallback(save_func, saved_values, saveat=0.1)
    sol = solve(newprob, Rosenbrock23(), saveat =0.1)#, callback=cb)

    return sol
end

function remake_odeprob(individual::Union{Vector, SubArray}, odeprob::ODEProblem)
    u0 = @view individual[14:end]
    p = @view individual[1:13]
    return remake(odeprob, u0=u0, p=p)
end



# For analysis of datasets 
function get_u0_params(individual_row::DataFrameRow)
    u0 = NamedTuple(individual_row[Between(:L, :A)]) |> pairs |> collect 
    params = NamedTuple(individual_row[Between(:kfᴸᴬ, :DF)]) |> pairs |> collect
    return u0, params
end

function solve_odes(individual_row::DataFrameRow, odeprob::ODEProblem = ODEPROB)
    u0, p = get_u0_params(individual_row)
    newprob = remake(odeprob, u0 = u0, p = p)
    return solve(newprob, Rodas4P(chunk_size = length(odeprob.u0)), saveat = 0.1, save_on = true)
end
