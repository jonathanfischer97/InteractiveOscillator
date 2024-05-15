const KF_RANGE = (1e-3, 1e2)
const KR_RANGE = (1e-3, 1e3)
const KCAT_RANGE = (1e-3, 1e3)
const DF_RANGE = (1e-1, 1e4)

const L_RANGE = (1e-1, 1e2)
const K_RANGE = (1e-2, 1e2)
const P_RANGE = (1e-2, 1e2)
const A_RANGE = (1e-1, 1e2)
const ZERO_SPECIES_RANGE = (0.0, 1e2)


const INPUT_BOUNDS = Dict([:kfᴸᴬ => KF_RANGE, :krᴸᴬ => KR_RANGE, :kfᴸᴷ => KF_RANGE, :krᴸᴷ => KR_RANGE, :kcatᴸᴷ => KCAT_RANGE, :kfᴸᴾ => KF_RANGE, :krᴸᴾ => KR_RANGE, :kcatᴸᴾ => KCAT_RANGE, :kfᴬᴷ => KF_RANGE, :krᴬᴷ => KR_RANGE, :kfᴬᴾ => KF_RANGE, :krᴬᴾ => KR_RANGE, :DF => DF_RANGE, :L => L_RANGE, :K => K_RANGE, :P => P_RANGE, :A => A_RANGE, :Lp => ZERO_SPECIES_RANGE, :LpA => ZERO_SPECIES_RANGE, :LK => ZERO_SPECIES_RANGE, :LpP => ZERO_SPECIES_RANGE, :LpAK => ZERO_SPECIES_RANGE, :LpAP => ZERO_SPECIES_RANGE, :LpAKL => ZERO_SPECIES_RANGE, :LpAPLp => ZERO_SPECIES_RANGE, :AK => ZERO_SPECIES_RANGE, :AP => ZERO_SPECIES_RANGE, :AKL => ZERO_SPECIES_RANGE, :APLp => ZERO_SPECIES_RANGE])


const INPUT_DEFAULTS = Dict([:kfᴸᴬ => 1.8438562425888723, :krᴸᴬ => 14.208807856126104, :kfᴸᴷ => 0.0015636561694145224, :krᴸᴷ => 79.75320912981275, :kcatᴸᴷ => 49.35726678915171, :kfᴸᴾ => 24.785882197512002, :krᴸᴾ => 122.13694437630977, :kcatᴸᴾ => 39.85674324994366, :kfᴬᴷ => 8.915650787390948, :krᴬᴷ => 0.03319110625610763, :kfᴬᴾ => 0.024471394984824354, :krᴬᴾ => 0.03517891784666366, :DF => 2582.7847577988523, :L => 43.15801010657622, :K => 21.026988187748856, :P => 5.484455252993349, :A => 9.83564656410484])


# const DEFAULT_INDIVIDUAL = [1.8438562425888723, 14.208807856126104, 0.0015636561694145224, 79.75320912981275, 49.35726678915171, 24.785882197512002, 122.13694437630977, 39.85674324994366, 8.915650787390948, 0.03319110625610763, 0.024471394984824354, 0.03517891784666366, 2582.7847577988523, 16.675796876637364, -16.646566665834015, 21.026988187748856, 5.4552250421899995, 43.15801010657622, 21.026988187748856, 5.484455252993349, 9.83564656410484, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

const FULL_INDIVIDUAL = [1.8438562425888723, 14.208807856126104, 0.0015636561694145224, 79.75320912981275, 49.35726678915171, 24.785882197512002, 122.13694437630977, 39.85674324994366, 8.915650787390948, 0.03319110625610763, 0.024471394984824354, 0.03517891784666366, 2582.7847577988523, 43.15801010657622, 21.026988187748856, 5.484455252993349, 9.83564656410484, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


const PARAMETER_SYMBOLS = [:kfᴸᴬ, :krᴸᴬ, :kfᴸᴷ, :krᴸᴷ, :kcatᴸᴷ, :kfᴸᴾ, :krᴸᴾ, :kcatᴸᴾ, :kfᴬᴷ, :krᴬᴷ, :kfᴬᴾ, :krᴬᴾ, :DF]
# const CONSERVED_CONSTANTS_SYMBOLS = [:Γ1, :Γ2, :Γ3, :Γ4]
# const PARAMETERS_AND_CONSTANTS_SYMBOLS = vcat(PARAMETER_SYMBOLS, CONSERVED_CONSTANTS_SYMBOLS)
const NONZERO_SPECIES_SYMBOLS = [:L, :K, :P, :A]
const ZERO_SPECIES_SYMBOLS = [:Lp, :LpA, :LK, :LpP, :LpAK, :LpAP, :LpAKL, :LpAPLp, :AK, :AP, :AKL, :APLp]
# const SOLVABLE_SPECIES_SYMBOLS = [:L, :K, :P, :A, :Lp, :LpA, :LK, :LpAK, :LpAP, :LpAKL, :LpAPLp, :AK]
const ALL_SPECIES_SYMBOLS = vcat(NONZERO_SPECIES_SYMBOLS, ZERO_SPECIES_SYMBOLS)

# const ODE_SYMBOLS = vcat(PARAMETERS_AND_CONSTANTS_SYMBOLS, SOLVABLE_SPECIES_SYMBOLS)
const ODE_SYMBOLS = vcat(PARAMETER_SYMBOLS, NONZERO_SPECIES_SYMBOLS)

# const IGNORED_SYMBOLS = setdiff(vcat(PARAMETERS_AND_CONSTANTS_SYMBOLS, ALL_SPECIES_SYMBOLS), ODE_SYMBOLS)

"""
    function make_fullrn()

Construct reaction network for the full oscillator model.
"""
function make_fullrn()

    fullrn = @reaction_network fullrn begin
        # unit = u"μM^-1*s^-1"; 
        #- Parameters to be optimized, both rate constants and DF (dimensionality factor)
        @parameters kfᴸᴬ = $(INPUT_DEFAULTS[:kfᴸᴬ]) [bounds = $(INPUT_BOUNDS[:kfᴸᴬ])] krᴸᴬ = $(INPUT_DEFAULTS[:krᴸᴬ]) [bounds = $(INPUT_BOUNDS[:krᴸᴬ])] kfᴸᴷ = $(INPUT_DEFAULTS[:kfᴸᴷ]) [bounds = $(INPUT_BOUNDS[:kfᴸᴷ])] krᴸᴷ = $(INPUT_DEFAULTS[:krᴸᴷ]) [bounds = $(INPUT_BOUNDS[:krᴸᴷ])] kcatᴸᴷ = $(INPUT_DEFAULTS[:kcatᴸᴷ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴷ])] kfᴸᴾ = $(INPUT_DEFAULTS[:kfᴸᴾ]) [bounds = $(INPUT_BOUNDS[:kfᴸᴾ])] krᴸᴾ = $(INPUT_DEFAULTS[:krᴸᴾ]) [bounds = $(INPUT_BOUNDS[:krᴸᴾ])] kcatᴸᴾ = $(INPUT_DEFAULTS[:kcatᴸᴾ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴾ])] kfᴬᴷ = $(INPUT_DEFAULTS[:kfᴬᴷ]) [bounds = $(INPUT_BOUNDS[:kfᴬᴷ])] krᴬᴷ = $(INPUT_DEFAULTS[:krᴬᴷ]) [bounds = $(INPUT_BOUNDS[:krᴬᴷ])] kfᴬᴾ = $(INPUT_DEFAULTS[:kfᴬᴾ]) [bounds = $(INPUT_BOUNDS[:kfᴬᴾ])] krᴬᴾ = $(INPUT_DEFAULTS[:krᴬᴾ]) [bounds = $(INPUT_BOUNDS[:krᴬᴾ])] DF = $(INPUT_DEFAULTS[:DF]) [description = "Dimensionality factor"; bounds = $(INPUT_BOUNDS[:DF])] 


        #- Species
        @species L(t) = $(INPUT_DEFAULTS[:L]) [description = "PIP"; bounds = $(INPUT_BOUNDS[:L])] K(t) = $(INPUT_DEFAULTS[:K]) [description = "PIP5K"; bounds = $(INPUT_BOUNDS[:K])] P(t) = $(INPUT_DEFAULTS[:P]) [description = "Synaptojanin"; bounds = $(INPUT_BOUNDS[:P])] A(t) = $(INPUT_DEFAULTS[:A]) [description = "AP2"; bounds = $(INPUT_BOUNDS[:A])] Lp(t) = 0.0 [description = "PIP2"; tunable = false] LpA(t) = 0.0 [description = "PIP2-AP2"; tunable = false] LK(t) = 0.0 [description = "PIP-Kinase"; tunable = false] LpP(t) = 0.0 [description = "PIP2-Phosphatase"; tunable = false] LpAK(t) = 0.0 [description = "PIP2-AP2-Kinase"; tunable = false] LpAP(t) = 0.0 [description = "PIP2-AP2-Phosphatase"; tunable = false] LpAKL(t) = 0.0 [description = "PIP2-AP2-Kinase-PIP"; tunable = false] LpAPLp(t) = 0.0 [description = "PIP2-AP2-Phosphatase-PIP2"; tunable = false] AK(t) = 0.0 [description = "AP2-Kinase"; tunable = false] AP(t) = 0.0 [description = "AP2-Phosphatase"; tunable = false] AKL(t) = 0.0 [description = "AP2-Kinase-PIP"; tunable = false] APLp(t) = 0.0 [description = "AP2-Phosphatase-PIP2"; tunable = false] 


        #* Phosphorylation Reactions
        (kfᴸᴷ, krᴸᴷ), L + K <--> LK                   # L binding to kinase
        kcatᴸᴷ, LK --> Lp + K                        # L phosphorylation by kinase into Lp

        (kfᴸᴷ * DF, krᴸᴷ), LpAK + L <--> LpAKL       # Membrane-bound K binds to L
        kcatᴸᴷ, LpAKL --> Lp + LpAK                  # L phosphorylation by membrane-bound kinase

        (kfᴸᴷ, krᴸᴷ), AK + L <--> AKL                # Peripheral: AK binds to L
        kcatᴸᴷ, AKL --> Lp + AK                      # L phosphorylation by peripheral AK

        #* Dephosphorylation Reactions
        (kfᴸᴾ, krᴸᴾ), Lp + P <--> LpP                # Lp binding to phosphatase
        kcatᴸᴾ, LpP --> L + P                        # L dephosphorylation by phosphatase

        (kfᴸᴾ * DF, krᴸᴾ), Lp + LpAP <--> LpAPLp     # Membrane-bound phosphatase binds to Lp
        kcatᴸᴾ, LpAPLp --> L + LpAP                  # L dephosphorylation by membrane-bound phosphatase

        (kfᴸᴾ, krᴸᴾ), AP + Lp <--> APLp              # Peripheral: AP binds to Lp
        kcatᴸᴾ, APLp --> L + AP                      # L dephosphorylation by peripheral AP

        #* Binding and Ternary Complex Formation
        (kfᴸᴬ, krᴸᴬ), Lp + A <--> LpA                # Lp binding to AP2

        (kfᴬᴷ, krᴬᴷ), LpA + K <--> LpAK              # LpA binding to kinase
        (kfᴬᴾ, krᴬᴾ), LpA + P <--> LpAP              # LpA binding to phosphatase

        (kfᴸᴷ * DF, krᴸᴷ), LpAK + L <--> LpAKL       # Membrane-bound K binds to L

        (kfᴸᴾ * DF, krᴸᴾ), Lp + LpAP <--> LpAPLp     # Membrane-bound P binds to Lp

        (kfᴬᴷ * DF, krᴬᴷ), LpA + LK <--> LpAKL       # Peripheral: LpA binds to LK
        (kfᴬᴾ * DF, krᴬᴾ), LpA + LpP <--> LpAPLp     # Peripheral: LpA binds to LpP

        #- Peripheral Ternary Complexes
        (kfᴸᴬ, krᴸᴬ), Lp + AK <--> LpAK
        (kfᴸᴬ * DF, krᴸᴬ), Lp + AKL <--> LpAKL
        (kfᴸᴬ, krᴸᴬ), Lp + AP <--> LpAP
        (kfᴸᴬ * DF, krᴸᴬ), Lp + APLp <--> LpAPLp
        (kfᴬᴷ, krᴬᴷ), A + K <--> AK
        (kfᴬᴾ, krᴬᴾ), A + P <--> AP
        (kfᴬᴷ, krᴬᴷ), A + LK <--> AKL
        (kfᴬᴾ, krᴬᴾ), A + LpP <--> APLp
    end  

    return fullrn
end

# function make_Amem_sys(osys::ODESystem)
#     @unpack L, K, P, A, Lp, LpA, LK, LpAK, LpAP, LpAKL, LpAPLp, AK, Γ = osys
#     t = independent_variable(osys) 
#     @variables Amem(t) = 0.0 [tunable = false]

#     Amem_eq = Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)

#     Amem_eq = Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + (LpAPLp + LpAKL + LpA + Lp - K - AK + L - P + Γ[2]) + (-LpAKL - LpAK - LK - K - AK + Γ[3]) + (-LpAP - 2.0LpAPLp - LpAKL - 2.0LpA - Lp - A + LK + 2.0K + AK - L + P + Γ[4]))

#     @named Amem_sys = ODESystem(Amem_eq, t)

#     return Amem_sys
# end

function make_Amem_sys(osys::ODESystem)
    @unpack L, K, P, A, Lp, LpA, LK, LpP, LpAK, LpAP, LpAKL, LpAPLp, AK, AP, APLp, AKL = osys
    t = independent_variable(osys) 
    @variables Amem(t) = 0.0 [tunable = false]

    Amem_eq = Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp)

    @named Amem_sys = ODESystem(Amem_eq, t)

    return Amem_sys
end

"""
    function make_ODEProblem(fixed_inputs = (; DF = 1000.0); tend::Float64=2000.)

Construct `ODEProblem` from `ReactionSystem` model, give named tuple to set parameters.
"""
function make_ODEProblem(; tend::Float64=1968.2) 

    rn = make_fullrn();

    osys = convert(ODESystem, rn, remove_conserved = false)

    Amem_sys = make_Amem_sys(osys)

    newsys = extend(Amem_sys, osys)
    newsys = structural_simplify(newsys, simplify = true, may_be_zero = true)

    ODEProblem{true, SciMLBase.FullSpecialize}(newsys, Float64[], (0.0, tend); jac = true, verbose=false, maxiters=1e6, abstol = 1e-8)
end


const ODEPROB = remake(make_ODEProblem(), u0 = view(FULL_INDIVIDUAL, 14:29), p = view(FULL_INDIVIDUAL, 1:13))


get_osys(odeprob::ODEProblem) = odeprob.f.sys


const OSYS = get_osys(ODEPROB)


const get_Amem = getu(ODEPROB, :Amem)

const get_L = getu(ODEPROB, :L)
const get_K = getu(ODEPROB, :K)
const get_P = getu(ODEPROB, :P)
const get_A = getu(ODEPROB, :A)
const get_LpA = getu(ODEPROB, :LpA)
const get_LpAK = getu(ODEPROB, :LpAK)
const get_LpAP = getu(ODEPROB, :LpAP)
const get_LpAKL = getu(ODEPROB, :LpAKL)
const get_LpAPLp = getu(ODEPROB, :LpAPLp)
const get_AK = getu(ODEPROB, :AK)
const get_AP = getu(ODEPROB, :AP)
const get_AKL = getu(ODEPROB, :AKL)
const get_APLp = getu(ODEPROB, :APLp)



# const RFFT_PLAN = plan_rfft(zeros(19683); flags = FFTW.PATIENT)


# save_func(u, t, integrator) = (u[6] + u[9] + u[10] + u[11] + u[12]) / (u[4] + u[6] + u[9] + u[10] + u[11] + u[12] + u[13] + u[14] + u[15] + u[16])

# save_func(u, t, integrator) = (u[8] + u[9] + u[10] + u[11] + u[6]) / (integrator.p[15] + integrator.p[16] + integrator.p[17])

# (LpAP(t) + LpAPLp(t) + LpAKL(t) + LpAK(t) + LpA(t)) / (Γ[2] + Γ[3] + Γ[4])






