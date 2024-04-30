
using NCDatasets
using Plots, Colors, LaTeXStrings, Measures
using DataFrames, Statistics

include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/save_utils.jl")
include("/home/lee/Dropbox/Development/NPZBDV_1D/src/rstar.jl")


function plot_state_vars(fsaven, season_num, bloom=false)

    H = 890
    zc = get_zc(H)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    ds.attrib["Season"] == "winter" ? season_num == 1 : season_num == 2
    season_num == 1 ? season = "Mesotrophic" : season = "Oligotrophic"

    rstar_b, _, _ = rstar_analysis(fsaven, season_num, bloom)

    if bloom==true
        N, P, Z, B, D, V, O = get_bloom_means(["n", "p", "z", "b", "d", "v", "o"], ds)
    else
        if ds["pulse"][:] == 1
            N, P, Z, B, D, V, O = get_endpoints(["n", "p", "z", "b", "d", "v", "o"], ds)
        else
            N, P, Z, B, D, V, O = mean_over_time(["n", "p", "z", "b", "d", "v", "o"], ds, season_num)
        end
    end
    
    f1, dir1 = plot_stacked_dar(N, P, Z, B, D, V, O, zc, filename)
    println("saving to: $(dir1)/$(filename)_dar.png")
    savefig(f1,"$(dir1)/$(filename)_dar.png")

    f2, dir2 = plot_individual_dar(P, Z, B, D, V, N, zc, filename)
    println("saving to: $(dir2)/$(filename)_dar.png")
    savefig(f2,"$(dir2)/$(filename)_dar.png")

    f3, dir3 = plot_combined_OSM(P, Z, B, D, V, N, zc, filename, rstar_b)
    println("saving to: $(dir3)/$(filename).png")
    savefig(f3,"$(dir3)/$(filename).png")

end


function plot_stacked_dar(N, P, Z, B, D, V, O, zc, filename)

    parent_folder = "results/plots/bmass_stacked/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls=7
    ab=1.0
    lfs = 8
    
    p1 = plot(sum(P, dims=2), -zc, lw=ls, lc="darkgreen", grid=false, xrotation=45, label=" Total P", 
    ylimits=yl, alpha=ab, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol ~N/m^3")
    plot!(sum(B, dims=2), -zc, lw=ls, lc="skyblue", label=" Total B", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(Z, dims=2), -zc, lw=ls, lc="red", label=" Total Z", alpha=ab, labelfontsize=lfs, legend=lg) 
    plot!(sum(V, dims=2), -zc, lw=ls, lc="purple", label=" Total V", alpha=ab, labelfontsize=lfs, legend=lg) 

    p2 = plot(sum(D, dims=2), -zc, lw=ls, lc="orange", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")

    p3 = plot(sum(N, dims=2), -zc, lw=ls, lc="grey", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")

    p4 = plot(O[:,1], -zc, lw=ls, lc="pink", ls=:dot, grid=false, xrotation=45, label="", ylimits=yl, 
    alpha=ab, labelfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")

    f1 = plot(p1, p2, p3, p4,
        layout = [1 1 1 1],
        size=(600,350),
        fg_legend = :transparent,
        title = ["Biomass" "OM" "N" "O"]
    )

    return f1, dir
end

function plot_individual_dar(P, Z, B, D, V, N, zc, filename)

    parent_folder = "results/plots/bmass_individual/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    ab=0.6
    lg=:bottomright
    tfs = 9
    ls=7
    ab=0.7
    lfs = 8
    ls2=5

    dcols = ["teal", "lightblue1", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "red2"]
    bcols = ["teal", "lightblue1", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "coral", "grey", "lime", "orchid", "pink2"]
    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "black", "brown", "wheat2", "mediumpurple3", "darkseagreen" ]
    tls = ["POM", "POM Consumers", "DOM Lab", "DOM S.Lab", "DOM Mid", "DOM S.Rec", "DOM Rec", "Phyto"]
    ncols = ["blue2"]

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    p1 = plot(D[:,1], -zc, lc=dcols[1], lw=ls2, linestyle=:dot, grid=false,  label=" Most Labile", xrotation=45, ylimits=yl, title=tls[1], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel="")
    plot!(D[:,2], -zc, lc=dcols[2], lw=ls2, linestyle=:dot, label=" ", labelfontsize=lfs, legend=lg)
    plot!(D[:,3], -zc, lc=dcols[3], lw=ls2, linestyle=:dot, label=" Least Labile", labelfontsize=lfs, legend=lg)

    p2 = plot(B[:,1], -zc, lc=bcols[1], lw=ls, grid=false,  label="", xrotation=45, ylimits=yl, title=tls[2], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,2], -zc, lc=bcols[2], lw=ls, label="", labelfontsize=lfs , legend=lg)
    plot!(B[:,3], -zc, lc=bcols[3], lw=ls, label="", labelfontsize=lfs , legend=lg)

    p3 = plot(D[:,4], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D4", xrotation=45, ylimits=yl, title=tls[3], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,4], -zc, lc=bcols[4], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,9], -zc, lc=bcols[9], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    p4 = plot(D[:,5], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D5", xrotation=45, ylimits=yl, title=tls[4], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel="")
    plot!(B[:,5], -zc, lc=bcols[5], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,10], -zc, lc=bcols[10], lw=ls, label=" Oli.", labelfontsize=lfs , legend=lg)

    p5 = plot(D[:,6], -zc, lc=dcols[9], lw=ls2, linestyle=:dot, grid=false, label=" D6", xrotation=45, ylimits=yl, title=tls[5], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, ylabel="Depth (m)", xlabel=L" mmol ~N/m^3")
    plot!(B[:,6], -zc, lc=bcols[6], lw=ls, label=" Cop", labelfontsize=lfs , legend=lg)
    plot!(B[:,11], -zc, lc=bcols[11], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    #NOTE changed D7 and D8 to remove resource line to see biomasss better
    p6 = plot(B[:,7], -zc, lc=bcols[7], lw=ls, grid=false, label=" Cop", xrotation=45, ylimits=yl, title=tls[6], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")
    plot!(B[:,12], -zc, lc=bcols[12], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    p7 = plot(B[:,8], -zc, lc=bcols[8], lw=ls, grid=false, label=" Cop", xrotation=45, ylimits=yl, title=tls[7], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, yformatter=Returns(""), xlabel=L" mmol ~N/m^3")
    plot!(B[:,13], -zc, lc=bcols[13], lw=ls, label=" Oli", labelfontsize=lfs , legend=lg)

    zcp = zc[1:40]
    p8 = plot(P[1:40,1], -zcp, lc=ncols[1], grid=false, lw=ls, label=" Oli", xrotation=45, title=tls[8], 
    titlefontsize=tfs, labelfontsize=lfs, legend=lg, xlabel=L" mmol ~N/m^3", alpha=ab)
    plot!(P[1:40,2], -zcp, lc=pcols[2], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)    
    plot!(P[1:40,3], -zcp, lc=pcols[3], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,4], -zcp, lc=pcols[4], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,5], -zcp, lc=pcols[5], lw=ls, label=" ",labelfontsize=lfs , legend=lg, alpha=ab)
    plot!(P[1:40,6], -zcp, lc=pcols[6], lw=ls, label=" Cop",labelfontsize=lfs , legend=lg, alpha=ab)


    f = plot(p1, p2, p3, p4, p5, p6, p7, p8,
            layout = [1 1 1 1 ; 1 1 1 1],
            fg_legend = :transparent,
            size=(700,600),
            # plot_title = "$season $type",
        )

    
    return f, dir

end

function plot_combined(P, Z, B, D, V, N, zc, filename, rstar)

    parent_folder = "results/plots/combined/"
    dir = check_subfolder_exists(filename, parent_folder)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls = 9
    ls2 = 4
    ls3 = 7
    ab = 0.7
    lfs = 7
    ls2 = 4
    xtfs = 8

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    dcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "honeydew3"]
    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    rscols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["POM", "POM Consumers", "DOM Lab", "DOM S.Lab", "DOM Mid", "DOM S.Rec", "DOM Rec", "Phyto"]

    fig1 = Array{Plots.Plot, 1}(undef, 12);
    fig1[1] =   plot(D[1:89, 4], -zc, lw=ls3, lc=dcols[9], label=" DOM1", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"log(mmol ~N/m^3)", 
                xrotation=45, title=tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[4][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[4])
                plot!(rstar[9][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[9])

    fig1[2] =   plot(D[1:89, 5], -zc, lw=ls3, lc=dcols[9], label=" DOM2", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", 
                xrotation=45, title=tls[4], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[5][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[5])
                plot!(rstar[10][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[10])

    fig1[3] =   plot(D[1:89, 6], -zc, lw=ls3, lc=dcols[9], label=" DOM3", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", 
                xrotation=45, title=tls[5], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[6][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[6])
                plot!(rstar[11][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[11])

    fig1[4] =   plot(D[1:89, 7], -zc, lw=ls3, lc=dcols[9], label=" DOM4", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", xrotation=45, title=tls[6], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[7][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[7])
                plot!(rstar[12][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[12])

    fig1[5] =   plot(D[1:89, 8], -zc, lw=ls3, lc=dcols[9], label=" DOM5", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", xrotation=45, title=tls[7], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(rstar[8][1:89], -zc, lw=ls2, label=" R*Cop", linestyle=:dot, lc=rscols[8])
                plot!(rstar[13][1:89], -zc, lw=ls2, label=" R*Oli", linestyle=:dot, lc=rscols[13])

    fig1[6] =   plot(D[1:89, 1], -zc, lw=ls, lc=dcols[1], label=" Lab", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab)
                plot!(D[1:89, 2], -zc, lw=ls, lc=dcols[2], label=" S.Lab", legendfontsize=lfs, alpha=ab)
                plot!(D[1:89, 3], -zc, lw=ls, lc=dcols[3], label=" Rec", legendfontsize=lfs, alpha=ab)
                plot!(rstar[1][1:89], -zc, lw=ls2, label=" R*B1", linestyle=:dot, lc=rscols[1])
                plot!(rstar[2][1:89], -zc, lw=ls2, label=" R*B2", linestyle=:dot, lc=rscols[2])
                plot!(rstar[3][1:89], -zc, lw=ls2, label=" R*B3", linestyle=:dot, lc=rscols[3])


    fig1[7] =   plot(B[1:89, 4], -zc, lw=ls, lc=bcols[4], label=" Cop", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L" mmol ~N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 9], -zc, lw=ls, lc=bcols[9], label=" Oli", legendfontsize=lfs)

    fig1[8] =   plot(B[1:89, 5], -zc, lw=ls, lc=bcols[5], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol ~N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 10], -zc, lw=ls, lc=bcols[10], label=" Oli", legendfontsize=lfs)

    fig1[9] =   plot(B[1:89, 6], -zc, lw=ls, lc=bcols[6], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol ~N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 11], -zc, lw=ls, lc=bcols[11], label=" Oli", legendfontsize=lfs)

    fig1[10] =   plot(B[1:89, 7], -zc, lw=ls, lc=bcols[7], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol ~N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 12], -zc, lw=ls, lc=bcols[12], label=" Oli", legendfontsize=lfs)

    fig1[11] =   plot(B[1:89, 8], -zc, lw=ls, lc=bcols[8], label=" Cop", legendfontsize=lfs, yformatter=Returns(""), xlabel=L" mmol ~N/m^3", xrotation=45, 
                grid=false, border=:box, legend=lg, xtickfontsize=xtfs)               
                plot!(B[1:89, 13], -zc, lw=ls, lc=bcols[13], label=" Oli", legendfontsize=lfs)

    fig1[12] =  plot(B[1:89, 1], -zc, lw=ls, lc=bcols[1], label=" B1", legendfontsize=lfs,
                ylabel="", xlabel=L" mmol ~N/m^3", xrotation=45, titlefontsize=tfs, grid=false, border=:box, legend=lg,
                yformatter=Returns(""), xtickfontsize=xtfs)               
                plot!(B[1:89, 2], -zc, lw=ls, lc=bcols[2], label=" B2", legendfontsize=lfs)
                plot!(B[1:89, 3], -zc, lw=ls, lc=bcols[3], label=" B3", legendfontsize=lfs)


    f1 = plot(fig1..., 
    fg_legend = :transparent,
    layout = (2,6),
    size=(900,700),
    )

    return f1, dir

end

function plot_combined_OSM(P, Z, B, D, V, N, zc, filename, rstar)

    parent_folder = "results/plots/combined/"
    dir = check_subfolder_exists(filename, parent_folder)

    P, Z, B = set_extinct_to_zero(P), set_extinct_to_zero(Z), set_extinct_to_zero(B)

    yl=(-890.0, 0)
    bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg = get_plot_vars()
    tfs = 9
    ls = 9
    ls2 = 4
    ls3 = 7
    ab = 0.7
    lfs = 7
    ls2 = 4
    xtfs = 8
    lg=:topleft

    dcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "honeydew3"]
    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    rscols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["Labile DOM", "Semi-Labile DOM", "Refractory DOM"]

    all_lab_D = sum(D[:, 4:6], dims=2)
    all_lab_cop = rstar[4][:] .+ rstar[5][:] .+ rstar[6][:]
    all_lab_oli = rstar[9][:] .+ rstar[10][:] .+ rstar[11][:]

    fig1 = Array{Plots.Plot, 1}(undef, 3);
    fig1[1] =   plot(all_lab_D, -zc, lw=ls3, lc=dcols[9], label=" L.DOM", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"log(mmol ~N/m^3)", 
                xrotation=45, title=tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, xscale=:log10, alpha=ab, 
                bottom_margin=5mm)
                plot!(all_lab_cop, -zc, lw=ls2, label=" R* Copio", linestyle=:dot, lc=rscols[4])
                plot!(all_lab_oli, -zc, lw=ls2, label=" R* Oligo", linestyle=:dot, lc=rscols[9])

    fig1[2] =   plot(D[1:89, 7], -zc, lw=ls3, lc=dcols[9], label=" SL.DOM", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", xrotation=45, title=tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab, bottom_margin=5mm)
                plot!(rstar[7][1:89], -zc, lw=ls2, label=" R* Copio", linestyle=:dot, lc=rscols[7])
                plot!(rstar[12][1:89], -zc, lw=ls2, label=" R* Oligo", linestyle=:dot, lc=rscols[12])

    fig1[3] =   plot(D[1:89, 8], -zc, lw=ls3, lc=dcols[9], label=" R.DOM", legendfontsize=lfs,
                yformatter=Returns(""), xlabel=L"log(mmol ~N/m^3)", xrotation=45, title=tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, 
                xtickfontsize=xtfs, xscale=:log10, alpha=ab, bottom_margin=5mm)
                plot!(rstar[8][1:89], -zc, lw=ls2, label=" R* Copio", linestyle=:dot, lc=rscols[8])
                plot!(rstar[13][1:89], -zc, lw=ls2, label=" R* Oligo", linestyle=:dot, lc=rscols[13])


    f1 = plot(fig1..., 
    fg_legend = :transparent,
    layout = (1,3),
    size=(600,300),
    )

    return f1, dir

end

bloom=false
# fsaven = "results/outfiles/240430_01:08_Su30yNP_6P3Z13B8D13V.nc"
fsaven = "results/outfiles/240430_10:54_Su30yNP_6P3Z13B8D13V.nc"

plot_state_vars(fsaven, 2, bloom)
