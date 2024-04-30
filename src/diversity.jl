using CSV, DataFrames, NCDatasets
using SparseArrays, LinearAlgebra
    

function load_spot_data()

    path = "NPZBD_1D/data/spot_diversity"
    spot_diversity = DataFrame(CSV.File("$(path)/spot_shannon_means.csv"))

    return spot_diversity

end


function get_proportions()

end

function get_richness()

end

function get_evenness()

end


function get_shannon_diversity()

end














spot_diversity = load_spot_data()