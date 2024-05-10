# DimensionalData.@dim Genes "Genes"
# DimensionalData.@dim Individuals "Individuals"


# cateogorical_lookup = DimensionalData.Categorical(ODE_SYMBOLS; order=DimensionalData.ForwardOrdered())

# const genes_dim = Genes(cateogorical_lookup)


# LOAD DATAFRAME
function load_df(filename::String = "smalldata.csv"; islocal = true)
    loadpath = islocal ? joinpath(pwd(), filename) : filename
    try
        data = CSV.read(loadpath, DataFrame)
        select!(data, Cols(DataFrames.Between(:kfᴸᴬ, :A) ,Cols([:Kmᴸᴷ ,:Kmᴸᴾ])))
        insertcols!(data, :Kmᴸᴷ, :Lp => 0.0)

        println("Data loaded!")
        return data
    catch e
        error("Failed to load data, check that you have downloaded it into this repository first: ", e)
    end
end



# function convert_to_dimarray(df::AbstractDataFrame)
#     data_matrix = Matrix(df[!, DataFrames.Between(:kfᴸᴬ, :A)]) |> permutedims

#     DimArray(data_matrix, (genes_dim, Individuals))
# end

