# using NCDatasets
# using Plots, ColorSchemes, LaTeXStrings
# using DataFrames
# using SparseArrays, LinearAlgebra, Statistics

include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/save_utils.jl")

#------------------------------------------------------------------------------
#                            TEST FOR EQUILIBRIUM
#------------------------------------------------------------------------------

function equilibrium_test(fsaven, season_num)

    ds = NCDataset(fsaven)
    parent_folder = "results/plots/equilib/"
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    season_num == 1 ? season = "Winter" : season = "Summer"

    f1 = equilib_depth_plots(ds, season)
    f2 = equilib_heatmaps(ds, season)

    savefig(f1, "$(dir)/$(filename)_tesy.png")

end


function equilib_heatmaps(ds, season)

end


function equilib_depth_plots(ds, season)

    final_2_years = final2(ds, ["p", "z", "b", "d"])
    final_yr1 = mean_penultimate_year(final_2_years)
    final_yr2 = mean_final_year(final_2_years)
    groups = ["P", "Z", "B", "D"]

    p1 = plot_equilib_depth(final_yr1[1], final_yr2[1], groups[1])
    p2 = plot_equilib_depth(final_yr1[2], final_yr2[2], groups[2])
    p3 = plot_equilib_depth(final_yr1[3], final_yr2[3], groups[3])
    p4 = plot_equilib_depth(final_yr1[4], final_yr2[4], groups[4])

    f1 = plot(p1, p2, p3, p4,
    layout = [1 1 ; 1 1],
    size=(1000,700),
    plot_title = "Equilib. States ($(season))",
    plot_titlefontsize = 24
    )

    return f1

end


function final2(ds, vars)

    final_2yrs = Vector{Any}()

    for v in vars
        append!(final_2yrs, [ds[v][:, :, end-14640:end]])
    end

    return final_2yrs

end


function mean_penultimate_year(ds)

    final2_yr1 = Vector{Any}()

    for i in 1:length(["p", "z", "b", "d"])
        append!(final2_yr1, [mean(ds[i][:, :, 1:7320], dims=3)])
    end

    return [final2_yr1[1], final2_yr1[2], final2_yr1[3], final2_yr1[4]]

end


function mean_final_year(ds)

    final2_yr2 = Vector{Any}()

    for j in 1:length(["p", "z", "b", "d"])
        append!(final2_yr2, [mean(ds[j][:, :, end-7320:end], dims=3)])
    end

    return [final2_yr2[1], final2_yr2[2], final2_yr2[3], final2_yr2[4]]

end


function get_pct_change(yr1, yr2)

    year1 = set_extinct_to_zero(yr1[:,:,:])
    year2 = set_extinct_to_zero(yr2[:,:,:])

    pct_chg = ((year2 .- year1)./ year1) * 100
    pct_chg .= ifelse.(isnan.(pct_chg), 0.0, pct_chg)

    return year1, year2, pct_chg

end


function plot_equilib_depth(yr1, yr2, group)

    H = 890
    zc = get_zc(H)
    n = get_size([yr1])
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs=9
    al=0.8
    xl=(-2.0, 2.0)

    if group == "P"
        cols = pcols
    elseif group == "B"
        cols = bcols
    elseif group == "Z"
        cols = zcols
    elseif group == "D"
        cols = dcols
    else
        cols = ncols
    end

    year1, year2, pct_chg = get_pct_change(yr1, yr2)

    p1 = plot(year1[:,1,end], -zc, grid=false, lw=ls, lc=cols[1], title="Penultimate Year", label="", xrotation=45,  
    titlefontsize=tfs, xlabel=L" mm/m^3", alpha=al, ylabel="Depth (m)")
    for i in 2:n[1]
        plot!(year1[:,i,end], -zc, lw=ls, lc=cols[i], label="", labelfontsize=lfs, alpha=al)
    end

    p2 = plot(year2[:,1,end], -zc, grid=false, lw=ls, lc=cols[1], title="Final Year", label="", xrotation=45, 
    titlefontsize=tfs, xlabel=L" mm/m^3", alpha=al, yformatter=Returns(""))
    for i in 2:n[1]
        plot!(year2[:,i,end], -zc, lw=ls, lc=cols[i], alpha=al, label="")
    end

    p3 = plot(pct_chg[:, 1], -zc, grid=false, lw=ls, lc=cols[1], labelfontsize=lfs, titlefontsize=tfs, xlabel=L" \% ", 
    title="% Change", alpha=al, yformatter=Returns(""), label=" $(group)1", legend=lg, xrotation=45)
    for i in 2:n[1]
        plot!(pct_chg[:, i], -zc, lw=ls, lc=cols[i], label=" $(group)$(i)", labelfontsize=lfs, alpha=al)
    end

    f = plot(p1, p2, p3,
    layout = [1 1 1],
    border=:box,
    fg_legend = :transparent,
    size=(600,450),
    )

    return f


end


# fsaven="results/outfiles/Wi100y_230928_23:10_8P6Z13B5D.nc"
# equilibrium_test(fsaven, 1)
# fsaven2="results/outfiles/endpoints/Su100y_230923_19:36_8P6Z13B5D_ep.nc"
# equilibrium_test(fsaven2, 2)

