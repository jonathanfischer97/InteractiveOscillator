"""
    function make_fullrn(; fixed_inputs = (; DF = 1000.0))

Construct reaction network for the full oscillator model. Give named tuple to set parameters.
"""
function make_fullrn(fixed_inputs = (; DF = 1000.0))

    kfrange = (1e-3, 1e2)
    krrange = (1e-3, 1e3)
    kcatrange = (1e-3, 1e3)
    DFrange = (1e-1, 1e4)

    Lrange = (1e-1, 1e2)
    Krange = (1e-2, 1e2)
    Prange = (1e-2, 1e2)
    Arange = (1e-1, 1e2)

    fullrn = @reaction_network fullrn begin
        # unit = u"μM^-1*s^-1"; 
        #- Parameters to be optimized, both rate constants and DF (dimensionality factor)
        @parameters kfᴸᴬ = 12.9 [bounds = $kfrange] krᴸᴬ = 6.1 [bounds = $krrange] kfᴸᴷ = 0.009 [bounds = $kfrange] krᴸᴷ = 2.3 [bounds = $krrange] kcatᴸᴷ = 832.7 [bounds = $kcatrange] kfᴸᴾ = 57.3 [bounds = $kfrange] krᴸᴾ = 0.04 [bounds = $krrange] kcatᴸᴾ = 42.2 [bounds = $kcatrange] kfᴬᴷ = 1.3 [bounds = $kfrange] krᴬᴷ = 0.006 [bounds = $krrange] kfᴬᴾ = 0.006 [bounds = $kfrange] krᴬᴾ = 0.9 [bounds = $krrange] DF = 1000.0 [description = "Dimensionality factor"; bounds = $DFrange]


        #- Species
        @species L(t) = 3.0 [description = "PIP"; bounds = $Lrange] K(t) = 0.5 [description = "PIP5K"; bounds = $Krange] P(t) = 0.3 [description = "Synaptojanin"; bounds = $Prange] A(t) = 2.0 [description = "AP2"; bounds = $Arange] Lp(t) = 0.0 [description = "PIP2"; tunable = false] LpA(t) = 0.0 [description = "PIP2-AP2"; tunable = false] LK(t) = 0.0 [description = "PIP-Kinase"; tunable = false] LpP(t) = 0.0 [description = "PIP2-Phosphatase"; tunable = false] LpAK(t) = 0.0 [description = "PIP2-AP2-Kinase"; tunable = false] LpAP(t) = 0.0 [description = "PIP2-AP2-Phosphatase"; tunable = false] LpAKL(t) = 0.0 [description = "PIP2-AP2-Kinase-PIP"; tunable = false] LpAPLp(t) = 0.0 [description = "PIP2-AP2-Phosphatase-PIP2"; tunable = false] AK(t) = 0.0 [description = "AP2-Kinase"; tunable = false] AP(t) = 0.0 [description = "AP2-Phosphatase"; tunable = false] AKL(t) = 0.0 [description = "AP2-Kinase-PIP"; tunable = false] APLp(t) = 0.0 [description = "AP2-Phosphatase-PIP2"; tunable = false] 
        
        

        #* ALIASES: L = PIP, Lp = PIP2, K = Kinase, P = Phosphatase, A = AP2 
        #* reactions between the same binding interfaces will have the same rate constant no matter the dimensionality or complex
        #* DF is the dimensionality factor that scales the 2D rate constants, and is proportional to the membrane surface to volume ratio
        (kfᴸᴷ,krᴸᴷ), L + K <--> LK # L binding to kinase
        kcatᴸᴷ, LK --> Lp + K # L phosphorylation by kinase into Lp
        (kfᴸᴬ,krᴸᴬ), Lp + A <--> LpA # Lp binding to AP2 adaptor 
        (kfᴬᴷ,krᴬᴷ), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
        (kfᴸᴷ*DF,krᴸᴷ), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by DF
        kcatᴸᴷ, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
        (kfᴸᴾ,krᴸᴾ), Lp + P <--> LpP # Lp binding to phosphatase 
        kcatᴸᴾ, LpP --> L + P # L dephosphorylation by phosphatase
        (kfᴬᴾ,krᴬᴾ), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 
        (kfᴸᴾ*DF,krᴸᴾ), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by DF
        kcatᴸᴾ, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality

        #* peripheral reactions, just need to be included in the model for conservation and balance
        (kfᴸᴬ,krᴸᴬ), Lp + AK <--> LpAK
        (kfᴸᴬ*DF,krᴸᴬ), Lp + AKL <--> LpAKL
        (kfᴸᴬ,krᴸᴬ), Lp + AP <--> LpAP
        (kfᴸᴬ*DF,krᴸᴬ), Lp + APLp <--> LpAPLp
        (kfᴬᴷ,krᴬᴷ), A + K <--> AK
        (kfᴬᴾ,krᴬᴾ), A + P <--> AP
        (kfᴬᴷ,krᴬᴷ), A + LK <--> AKL
        (kfᴬᴾ,krᴬᴾ), A + LpP <--> APLp
        (kfᴬᴷ*DF,krᴬᴷ), LpA + LK <--> LpAKL
        (kfᴬᴾ*DF,krᴬᴾ), LpA + LpP <--> LpAPLp
        (kfᴸᴷ,krᴸᴷ), AK + L <--> AKL #binding of kinase to lipid
        kcatᴸᴷ, AKL --> Lp + AK #phosphorylation of lipid
        (kfᴸᴾ,krᴸᴾ), AP + Lp <--> APLp #binding of phosphatase to lipid
        kcatᴸᴾ, APLp --> L + AP #dephosphorylation of lipid
    end  

    setdefaults!(fullrn, collect(pairs(fixed_inputs)))

    #- Get dictionary mapping of symbols to variables
    symtovar_dict = fullrn.var_to_name

    #- Set the tunable metadata to false for all fixed variables
    for fixed_sym in keys(fixed_inputs)
        fixed_var = symtovar_dict[fixed_sym] #* get the corresponding Symbolics variable
        fixed_var = setmetadata(fixed_var, VariableTunable, false) #*set the tunable metadata to false
    end

    return fullrn
end



"""
    function make_ODEProblem(fixed_inputs = (; DF = 1000.0); tend::Float64=2000.)

Construct `ODEProblem` from `ReactionSystem` model, give named tuple to set parameters.
"""
function make_ODEProblem(fixed_inputs = (; DF = 1000.0); tend::Float64=2187.0)

    rn = make_fullrn(fixed_inputs);

    tspan = (0.0, tend)

    osys = convert(ODESystem, rn)

    ODEProblem(osys, Float64[], tspan, jac = true)
end


indexof(sym, syms)::Int = findfirst(isequal(sym), syms)

get_variables(osys) = vcat(parameters(osys), states(osys))


function get_tunable_variables(osys) 
    #* Get vector of all parameters and states
    variables = get_variables(osys)

    #* Get vector of all tunable parameters and states
    filter!(x -> istunable(x, true), variables)

    return variables
end


get_tunable_bounds(osys) = getbounds(get_tunable_variables(osys))


function get_Amem_indices(sys::ODESystem)::Vector{Int}
    #* Get indices of membrane-bound AP2 species to tell the solver to save them
    @unpack LpA, LpAK, LpAP, LpAKL, LpAPLp, AKL, APLp = sys
    Amem_species = (LpA, LpAK, LpAP, LpAKL, LpAPLp, AKL, APLp)
    [indexof(var, states(sys)) for var in Amem_species]
end


