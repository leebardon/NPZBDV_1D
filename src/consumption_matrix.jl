# using Statistics, StatsBase
# using SparseArrays, JLD

# g = gmax * prey/(prey + k)
# g = gmax * (P + B)/ (P + B) + k
# dPdt = -g*Z(t) * GrM(k, 1:np)* P/prey
# dBdt = -g*Z(t) * GrM(k, np+1:end) * B/prey 


function build_consumption_matrix(resources, consumers)

    if resources == 1
        CM = convert(Array{Bool}, ones(resources, consumers) .== 1)
    else 
        CM1 = generalists_matrix(resources, consumers)
        CM = convert(Array{Bool}, CM1 .== 1)
    end
    
    return CM

end

function generalists_matrix(resources, consumers)

    # Method 1: assign num substrates first. 
    # nupr = n substrates consumed by each consumer 
    xhigh, xlow = log10(resources), log10(1)
    nupr = 10 .^ ((rand(consumers).-0.5).*(xhigh-xlow).+mean([xhigh xlow]))
    nupr = round.(Int,nupr) #round to integers:
    wv = StatsBase.ProbabilityWeights([1/resources:1/resources:1;])
    food = [sample(1:resources,wv,nupr[j], replace = false) for j=1:consumers]

    CM = zeros(resources,consumers)
    for j = 1:consumers
        CM[food[j],j] .= 1
    end

    return CM

end



#NOTE below is code for specialistsa vs generalists consumtion matrix. Results inconclusive, needs further work
# function build_consumption_matrix(resources, consumers)

#     if resources == 1
#         CM = convert(Array{Bool}, ones(resources, consumers) .== 1)
#     else 
#         CMone = specialists_matrix(resources, consumers)
#         CMtwo = generalists_matrix(resources, consumers)
#         CM = convert(Array{Bool}, hcat(CMone,CMtwo) .== 1)
#     end
    
#     return CM

# end
 

# #TODO - these only work if mod(nb/2) == 0
# function specialists_matrix(resources, consumers)
#     """TODO need to square the matrix before halfing, or we end up with lots of specialists that consume nothing

#     e.g. 
#     1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  1  0  0
#     0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1
#     0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0
#     0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  1  0  0  1  0  0  0  0  0  0  0  0  1  0  0
#     0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  0  1  0  1  0  0  1  0  1
#     0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  0  1  0  1  1  1  1  1  0  0  1  1  0  1  0  1
#     0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  1  1  1  1  0  0  1  0  1  0  0  0  1  1  1
#     0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  0  1  1  1  1  0  0  0  0  0  1  0  0
#     0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  1  0  1  0  1  0  1  1  1  0  1  0  0  1  1  1  1  0  1
#     0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  1  1  1  1  0  0  1  0  0  0  1  1  1  0  1

#     """
 
#     CMone = zeros(resources, consumers÷2)
#     for j = 1:resources
#         CMone[j,j] = 1
#     end

#     return CMone

# end


# function generalists_matrix(resources, consumers)

#     # Method 1: assign num substrates first. 
#     # nupr = n substrates consumed by each consumer 
#     xhigh, xlow = log10(resources), log10(1)
#     nupr = 10 .^ ((rand(consumers÷2).-0.5).*(xhigh-xlow).+mean([xhigh xlow]))
#     nupr = round.(Int,nupr) #round to integers:
#     wv = StatsBase.ProbabilityWeights([1/resources:1/resources:1;])
#     food = [sample(1:resources,wv,nupr[j], replace = false) for j=1:consumers÷2]

#     CMtwo = zeros(resources,consumers÷2)
#     for j = 1:consumers÷2
#         CMtwo[food[j],j] .= 1
#     end

#     return CMtwo

# end








