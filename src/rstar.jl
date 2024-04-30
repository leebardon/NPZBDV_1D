# using NCDatasets
# using Plots, ColorSchemes, LaTeXStrings
# using DataFrames, CSV
# using SparseArrays, LinearAlgebra, Statistics

include("utils/utils.jl")
include("utils/save_utils.jl")
include("plotting/rstar_plots.jl")


function rstar_analysis(fsaven, season_num, bloom)

    ds = NCDataset(fsaven)
    season_num == 1 ? season = "Winter" : season = "Summer"

    if bloom==true
        N, P, Z, B, D, O = get_bloom_means(["n", "p", "z", "b", "d", "o"], ds)
    else
        if ds["pulse"][:] == 1
            N, P, Z, B, D, O = get_endpoints(["n", "p", "z", "b", "d", "o"], ds)
        else
            N, P, Z, B, D, O = mean_over_time(["n", "p", "z", "b", "d", "o"], ds, season_num)
        end
    end

    # P = set_extinct_to_zero(P)
    # B = set_extinct_to_zero(B)
    # Z = set_extinct_to_zero(Z)

    rstar_b, rstar_p, rstar_z = get_rstar(B, P, Z, ds)
    # plot_rstar(rstar_b, rstar_p, rstar_z, fsaven)
    plot_rstar_dar(rstar_b, rstar_p, rstar_z, fsaven)

    return rstar_b, rstar_p, rstar_z

end 


function get_rstar(B, P, Z, ds)

    nb, nz, np = get_size([B, Z, P])
    
    mort_b = mortality(B, ds, nb, "B")
    grz_b = b_grazing(B, Z, np, nb, nz, ds)
    loss_b = loss(B, mort_b, grz_b, nb)
    rstar_b = RstarB(loss_b, ds)

    mort_p = mortality(P, ds, np, "P")
    grz_p = p_grazing(P, Z, ds, np, nz)
    loss_p = loss(P, mort_p, grz_p, np)
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

function mortality(Biomass, ds, n, group)

    ngrid = length(Biomass[:,1])
    mort = zeros(Float64, ngrid, n) 

    if group == "B"
        for i in range(1, n)
            mort[:, i] += (ds["m_lb"][i] .+ ds["m_qb"][i] .* Biomass[:,i])
        end
    elseif group == "P"
        for i in range(1, n)
            m_lp = ones(n) * 0 
            m_qp = ones(n) * 0.1
            mort[:, i] += (m_lp[i] .+ m_qp[i] .* Biomass[:,i])
            # mort[:, i] += (ds["m_lp"][i] .+ ds["m_qp"][i] .* Biomass[:,i])
        end
    elseif group == "Z"
        for i in range(1, n)
            mort[:, i] += (ds["m_lz"][i] .+ ds["m_qz"][i] .* Biomass[:,i])
        end
    end

    return mort

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

function loss(Biomass, mortality, grazing, n)

    ngrid = length(Biomass[:,1])
    loss = zeros(Float64, ngrid, n)

    for i in range(1, n)
        loss[:,i] = mortality[:,i] .+ grazing[:,i]
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
    temp_mod = get_temp_mod(ds)
    II, JJ = get_nonzero_axes(ds["CM"][:,:])

    RS = Any[]
    for j = axes(II, 1)
        push!(RS, Km_ij[II[j],JJ[j]] .* loss[:, j] ./ (yield[II[j],JJ[j]] .* umax_ij[II[j],JJ[j]] .* temp_mod .- loss[:, j]))
    end

    RS_out = check_for_negatives(RS)

    # RS_out = replace!(RS, NaN => 0.0)

    return RS_out


end


function RstarP(loss, ds, np)

    vmax_ij = ds["vmax_ij"][:,:]
    Kp_ij = ds["Kp_ij"][:,:]
    temp_mod = get_temp_mod(ds)
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




# fsaven = "results/outfiles/Wi50y_231202_15:33_6P3Z13B8D.nc" # meso steady 50yrs
# fsaven = "results/outfiles/Wi100y_231202_16:19_6P3Z13B8D.nc" # meso steady 100yrs
# fsaven = "results/outfiles/Wi50y_231202_16:48_6P3Z13B8D.nc" # meso pulse 50yrs
# fsaven = "results/outfiles/Wi100y_231202_17:10_6P3Z13B8D.nc" # meso pulose 100yrs

# fsaven="results/outfiles/Wi50y_231202_23:38_6P3Z13B8D.nc"
# fsaven = "results/outfiles/Wi50y_231203_00:10_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi50y_231203_00:42_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi50y_231203_01:17_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi50y_231203_01:42_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi100y_231203_10:39_6P3Z13B8D.nc"
# fsaven="results/outfiles/Wi100y_231203_11:03_6P3Z13B8D.nc"

# fsaven="results/outfiles/Wi100y_231203_10:39_6P3Z13B8D.nc" #winter steady
# fsaven="results/outfiles/Wi100y_231203_11:03_6P3Z13B8D.nc" #winter pulse
# fsaven="results/outfiles/Su100y_231203_14:48_6P3Z13B8D.nc" #summer steady
# fsaven="results/outfiles/Su100y_231203_19:58_6P3Z13B8D.nc"  #summer pulse

# fsaven="results/outfiles/Wi50y_231212_16:28_6P3Z13B8D.nc"
# rstar_analysis(fsaven, 1)

# fsaven="results/outfiles/Su50y_231212_12:30_6P3Z13B8D.nc"
# rstar_analysis(fsaven, 2)
