using NCDatasets
using Plots, ColorSchemes, LaTeXStrings
using DataFrames, CSV
using SparseArrays, LinearAlgebra, Statistics

include("utils/utils.jl")
include("utils/save_utils.jl")
include("plotting/rstar_plots.jl")


function rstar_analysis(fsaven, season_num, lysis, graze, bloom)

    ds = NCDataset(fsaven)
    season_num == 1 ? season = "Winter" : season = "Summer"

    if bloom==true
        N, P, Z, B, D, V, O = get_bloom_means(["n", "p", "z", "b", "d", "v", "o"], ds)
    else
        if ds["pulse"][:] == 1
            N, P, Z, B, D, V, O = get_endpoints(["n", "p", "z", "b", "d", "v", "o"], ds)
        else
            N, P, Z, B, D, V, O = mean_over_time(["n", "p", "z", "b", "d", "v", "o"], ds, season_num)
        end
    end

    P = set_extinct_to_zero(P)
    B = set_extinct_to_zero(B)
    Z = set_extinct_to_zero(Z)
    # V = set_extinct_to_zero(V)

    rstar_b, rstar_p, rstar_z = get_rstar(B, P, Z, V, lysis, graze, ds)
    # plot_rstar(rstar_b, rstar_p, rstar_z, fsaven)
    plot_rstar_dar(rstar_b, rstar_p, rstar_z, fsaven)

    return rstar_b, rstar_p, rstar_z

end 


function get_rstar(B, P, Z, V, lysis, graze, ds)

    nb, nz, np, nv = get_size([B, Z, P, V])

    # B loss from lysis (if explicitly included) added to mort
    mort_b = mortality(B, V, ds, nb, lysis, graze, "B")
    mort_p = mortality(P, V, ds, np, lysis, graze, "P")

    # calc b_grazing only if Z and  calc b_lysis only if V
    if graze == 1 
        grz_b = b_grazing(B, Z, np, nb, nz, ds) 
        grz_p = p_grazing(P, Z, ds, np, nz)
    else    
        grz_b, grz_p = 0, 0
    end
    

    loss_b = loss(B, mort_b, grz_b, graze, nb)
    rstar_b = RstarB(loss_b, ds)

    loss_p = loss(P, mort_p, grz_p, graze, np)
    rstar_p = RstarP(loss_p, ds, np)

    #TODO fix rstar z calcs (check yield in particular, should not be y_ij)
    # loss_z = mortality(Z, ds, nz, "Z")
    # rstar_z = RstarZ(loss_z, ds, nz)
    rstar_z = NaN

    return rstar_b, rstar_p, rstar_z

end


#-----------------------------------------------------------------------------------
#                                     MORTALITY
#-----------------------------------------------------------------------------------

function mortality(Biomass, V, ds, n, lysis, graze, group)

    ngrid = length(Biomass[:,1])
    mort = zeros(Float64, ngrid, n) 

    if group == "B"
        if lysis == 1
            for i in range(1, n)
                lysis_Bi = get_lysis_Bi(ds, Biomass[:, i], V[:, i])
                mort[:, i] += (ds["m_lb"][i] .* Biomass[:,i]) .+ lysis_Bi
            end
        else
            for i in range(1, n)
                mort[:, i] += (ds["m_lb"][i] .+ ds["m_qb"][i] .* Biomass[:,i])
            end
        end
    end

    if group == "P"
        for i in range(1, n)
                # mort[:, i] += (m_lp[i] .+ m_qp[i] .* Biomass[:,i])
            mort[:, i] += (ds["m_lp"][i] .+ ds["m_qp"][i] .* Biomass[:,i])
        end
    end

    if group == "Z"
        if graze == 1
            for i in range(1, n)
                mort[:, i] += (ds["m_lz"][i] .+ ds["m_qz"][i] .* Biomass[:,i])
            end
        end
    end

    return mort

end


function get_lysis_Bi(ds, Bi, Vi)
       
    lysis_Bi = ds["vly"] .* Bi .* Vi
    lysis_Bi = replace!(lysis_Bi, NaN => 0.0)
        
    return lysis_Bi

end

#-----------------------------------------------------------------------------------
#                                     GRAZING
#-----------------------------------------------------------------------------------

function b_grazing(B, Z, np, nb, nz, ds)

    GrM = ds["GrM"][:,:]
    g_max = ds["g_max"][:]
    K_g = ds["K_g"][:]
    ngrid = length(B[:,1])

    grazing = Any[]
    grazing_zi = zeros(Float64, ngrid, nb) 
    for z in range(1, nz)
        if sum(GrM[z,np+1:end]) > 0 
            prey = GrM[z,np+1:end]' .*B[:,1:end]
            g = g_max[z] .* prey ./ (prey .+ K_g[z])
            grazing_zi += (g .* Z[:,z] .* GrM[z,np+1:end]') ./ prey
            grazing_zi = replace!(grazing_zi, NaN => 0.0)
            push!(grazing, grazing_zi)
        end   
    end

    return sum(grazing)

end


function p_grazing(P, Z, ds, np, nz)

    GrM = ds["GrM"][:,:]
    grazing = Any[]
    g_max = ds["g_max"][:]
    K_g = ds["K_g"][:]

    grazing_zi = zeros(Float64, length(P[:,1]), np) 
    for z in range(1, nz)
        if sum(GrM[z,1:np]) > 0 
            prey = GrM[z,1:np]' .* P[:,1:end]
            g = g_max[z] .* prey ./ (prey .+ K_g[z])
            grazing_zi += (g .* Z[:,z] .* GrM[z,1:np]') ./ prey
            grazing_zi = replace!(grazing_zi, NaN => 0.0)
            push!(grazing, grazing_zi)
        end   
    end

    return sum(grazing)

end



#-----------------------------------------------------------------------------------
#                                     TOTAL LOSS
#-----------------------------------------------------------------------------------
# MORTALITY + GRAZING

function loss(Biomass, mortality, grazing, graze, n)

    ngrid = length(Biomass[:,1])
    loss = zeros(Float64, ngrid, n)

    for i in range(1, n)
        loss[:,i] = mortality[:,i] 
    end

    if graze == 1 
        for i in range(1, n)
            loss[:,i] += grazing[:,i]
        end
    end

    return loss
end



#-----------------------------------------------------------------------------------
#                                     RSTAR CALCS
#-----------------------------------------------------------------------------------

function RstarB(loss, ds)

    umax_ij = ds["umax_ij"][:,:]
    Km_ij = ds["Km_ij"][:,:]
    yield = ds["y_ij"][:,:]
    temp_mod = ds["temp_fun"][:]
    II, JJ = get_nonzero_axes(ds["CM"][:,:])

    RS = Any[]
    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* umax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
    end

    # When R* plots are logscale, replace 0 with NaN or plots won't work
    RS_out = check_for_negatives(RS)
    RS_out = replace.(RS, 0.0 => NaN)


    return RS_out


end


function RstarP(loss, ds, np)

    vmax_ij = ds["vmax_ij"][:,:]
    Kp_ij = ds["Kp_ij"][:,:]
    temp_mod = ds["temp_fun"][:]
    II, JJ = get_nonzero_axes(ds["CMp"][:,:])

    RS = Any[]
    for j = axes(II, 1)
        rs_j = Kp_ij[II[j],JJ[j]] .* loss[:, j] ./ (vmax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j])
        RS_j = replace!(rs_j, NaN => 0.0)
        push!(RS, RS_j)
    end

    RS_out = check_for_negatives(RS)

    return RS_out

end


function RstarZ(loss, ds, nz)

    g_max = ds["g_max"][:]
    K_g = ds["K_g"][:]
    yield = ds["y_ij"][1]
    temp_mod = get_temp_mod(ds)

    RS = Any[]
    for z in range(1, nz)
        push!(RS, K_g[z] .* loss[:, z] ./ (yield[z] * g_max[z] .* temp_mod .- loss[:, i]))
    end

    RS_out = check_for_negatives(RS)

    return RS_out

end

# lysis=1
# graze=2
# season=1
# bloom=false
# fsaven="results/outfiles/240501_14:06_Wi2yNP_6P3Z13B8D13V.nc"
# rstar_analysis(fsaven, season, lysis, graze, bloom)

