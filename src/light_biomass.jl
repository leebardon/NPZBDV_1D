# #Based on: https://bg.copernicus.org/articles/12/4447/2015/

using DataFrames, NCDatasets
using SparseArrays, LinearAlgebra
using Plots, ColorSchemes, LaTeXStrings

using DelimitedFiles


include("utils/utils.jl")



function phyto_uptake(ds, N, P)
    dz = 10
    H = 890
    zc = [dz/2 : dz : H - dz/2;]    # centered depth 

    euz = 25                        # euphotic zone lengthscale
    I_in = 700                      # Incident radiation at surface W/m2

    CMp = ds["CMp"][:,:]
    temp_fun = ds["temp_fun"][:]
    vmax_ij = ds["vmax_ij"][:,:]
    Kp_ij = ds["Kp_ij"][:,:]
    K_I = 10.0         
    e_o = 150/16   

    dNdt = zeros(89, 1)
    dPdt = zeros(89, 6)
    dOdt = zeros(89, 1)
    II, JJ = get_nonzero_axes(CMp)
  
    light, light_euz_only = light_attenuation(P, 89, zc, euz, I_in)

    plot_light(light_euz_only, light, zc)

    for j = axes(II, 1)
        uptake = P[:,JJ[j]] .* temp_fun .* vmax_ij[II[j],JJ[j]] .* min.(N ./ (N .+ Kp_ij[II[j],JJ[j]]), light ./ (light .+ K_I))
        dNdt += -uptake
        dOdt += uptake * e_o
        dPdt[:,JJ[j]] += uptake 
    end

    return dPdt, dNdt, dOdt

end 


function light_attenuation(P, dz, zc, euz, I_in)

    a_chlD = 0.04                # Absorption coeff incorporating Chl-a and CDOM (m2/mg Chl) (Zakem et al 2015)
    chl2c_max = 0.2              # max chlorophyll to carbon ratio (mg Chl/mmol C)
    chl2c_min = 0.02
    #NOTE chl2c is variable in emily's model, ranging from 0.02 at top to 0.2
    kw = 0.04                    # attenuation coeff of water (m2/mg Chl)

    light = zeros(dz)
    # chl_tot = sum(P, dims=2) .* chl2c_max .* 6.6   # mmolN/m3 * mgChl/mmolC (redfield ratio -> 6.6 C for every N)
    chl_tot = sum(P, dims=2) .* chl2c_min .* 6.6 

    for d in range(1, dz)
		light[d] = I_in * exp(-(zc[d])*(kw + sum(chl_tot[1:d] * a_chlD)))
    end

    light_euz_only = I_in .* exp.( -zc ./ euz)

    return light, light_euz_only

end



function plot_light(old_light, light, zc)

    ls=7
    ab=0.4
    lfs = 8

    d=20 #plot top 200m only
    
    p1 = plot(old_light[1:d], -zc[1:d], lw=ls, lc="black", grid=false, xrotation=45, label=" Euz only", 
    alpha=ab, labelfontsize=lfs, ylabel="Depth (m)", xlabel=L"W/m^2")
    plot!(light[1:d], -zc[1:d], lw=ls, lc="darkgreen", label=" B'mass scaled", alpha=ab, labelfontsize=lfs)

    f1 = plot(p1, 
    layout = [1],
    size=(250,400),
    fg_legend = :transparent,
    title = "Light Attenuation"
    )

    # p2 = plot(light[1:d], -zc[1:d], lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" B'mass", 
    # alpha=ab, labelfontsize=lfs, ylabel="Depth (m)", xlabel=L"W/m^2")

    # f1 = plot(p1, p2, 
    # layout = [1 1],
    # size=(500,350),
    # fg_legend = :transparent,
    # title = ["Euz only" "B'mass scaled"]
    # )

    savefig(f1,"light_compared_chl2cmax.png")
end

# ds = NCDataset("results/outfiles/endpoints/Wi50y_231202_16:19_6P3Z13B8D_ep.nc")
# N = ds["n"][:,:]
# P = ds["p"][:,:]

ds = NCDataset("results/outfiles/Wi50y_231214_20:10_6P3Z13B8D.nc")
N = ds["n"][:,:,end]
P = ds["p"][:,:,end]

dPdt, dNdt, dOdt = phyto_uptake(ds, N, P)