# LOAD DATA
function load_data()
    loadpath = joinpath(pwd(), "smalldata.csv")
    try
        data = CSV.read(loadpath, DataFrame)
        println("Data loaded!")
        return data
    catch e
        error("Failed to load data, check that you have downloaded it into this repository first: ", e)
    end
end