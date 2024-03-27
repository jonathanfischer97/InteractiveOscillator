function compute_Amem(sol::ODESolution)::Vector{Float64}
    initialAP2 = sol.prob.u0[4] #* initial AP2 concentration
    #* sum all AP2 species on the membrane and normalize by initial AP2 
    map(sum, sol.u) ./ initialAP2 
end


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

function make_ode_solver()
    odeprob = make_ODEProblem((;))
    make_ode_solver(odeprob)
end




















