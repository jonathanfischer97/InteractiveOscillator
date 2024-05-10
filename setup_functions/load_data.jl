# LOAD DATA
function load_data(filename::String = "smalldata.csv")
    loadpath = joinpath(pwd(), filename)
    try
        data = CSV.read(loadpath, DataFrame)
        select!(data, Cols(Between(:kfᴸᴬ, :A) ,Cols([:Kmᴸᴷ ,:Kmᴸᴾ])))

        println("Data loaded!")
        # add_zero_species!(data)
        return data
    catch e
        error("Failed to load data, check that you have downloaded it into this repository first: ", e)
    end
end

# # Add columns for zero species to dataframe
# function add_zero_species!(data)
#     zerocols = ZERO_SPECIES_SYMBOLS .=> zeros(length(ZERO_SPECIES_SYMBOLS))
#     insertcols!(data, :Kmᴸᴷ, zerocols...)
#     # for sym in ZERO_SPECIES_SYMBOLS
#     #     data[!, sym] .= 0.0
#     # end
#     return nothing
# end