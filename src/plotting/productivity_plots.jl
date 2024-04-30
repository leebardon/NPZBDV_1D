# using Plots, LatexStrings

function plot_seasonal_BP(BP_model, PP_model, BP_Leu, BP_Thy, PP, ds, fsaven, season_num)

    H = 890
    zc = get_zc(H)

    filename = replace(fsaven, ".nc" => "", "results/outfiles/" => "")
    parent_folder = "results/plots/productivity/"
    dir = check_subfolder_exists(filename, parent_folder)

    season_num == 1 ? season = "Meso." : season = "Oligo."
    leu_spot = dropmissing(BP_Leu, disallowmissing=true)
    thy_spot = dropmissing(BP_Thy, disallowmissing=true)

    tfs = 9 
    ls = 5
    ls2 = 7
    ab = 0.5
    lfs = 7
    ls3 = 3
    xtfs = 8
    lg = :bottomright

    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "mediumpurple3"]
    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["\nPOM Consumers (model)", "\nDOM Consumers (model)", "\nTotal BP", "\nPhytoplankton (model)","\nTotal PP"]

    p1 = plot(BP_model[:, 1], -zc, lw=ls2, lc=bcols[1], label=" B1", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"cells/mL/day", 
        xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        plot!(BP_model[:, 2], -zc, lw=ls2, lc=bcols[2], label=" B2", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 3], -zc, lw=ls2, lc=bcols[3], label=" B3", legendfontsize=lfs, alpha=ab)

    p2 = plot(BP_model[:, 4], -zc, lw=ls2, lc=bcols[4], label=" B4c", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"cells/mL/day" ,
        xrotation=45, title =tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 9], -zc, lw=ls2, lc=bcols[9], label=" B9o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 5], -zc, lw=ls2, lc=bcols[5], label=" B5c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 10], -zc, lw=ls2, lc=bcols[10], label=" B10o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 6], -zc, lw=ls2, lc=bcols[6], label=" B6c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 11], -zc, lw=ls2, lc=bcols[11], label=" B11o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 7], -zc, lw=ls2, lc=bcols[7], label=" B7c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 12], -zc, lw=ls2, lc=bcols[12], label=" B12o", legendfontsize=lfs, alpha=ab)
        plot!(BP_model[:, 8], -zc, lw=ls2, lc=bcols[8], label=" B8c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(BP_model[:, 13], -zc, lw=ls2, lc=bcols[13], label=" B13o", legendfontsize=lfs, alpha=ab)

    p3 = plot(sum(BP_model[:, 1:3], dims=2), -zc, lw=ls2, lc="cyan", label=" POM", legendfontsize=lfs, yformatter=Returns(""),  xlabel=L"cells/mL/day" ,
        xrotation=45, title =tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        plot!(sum(BP_model[:, 4:13], dims=2), -zc, lw=ls2, lc="purple", label=" DOM", legendfontsize=lfs, alpha=ab)
        plot!(sum(BP_model[:, :], dims=2), -zc, lw=ls, lc="red", label=" Total BP", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        scatter!(leu_spot.mean_Leu, -leu_spot.depth,  markersize=7, markercolor="yellow2", markerstrokecolor=:black, markershape=:diamond, 
        label=" SPOT")

    maxdep=20
    p4 = plot(PP_model[1:maxdep, 1], -zc[1:maxdep], lw=ls2, lc=pcols[1], label=" P1c", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"mg~C/m^{3}/day" ,
        xrotation=45, title =tls[4], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, linestyle=:dot)
        plot!(PP_model[1:maxdep, 2], -zc[1:maxdep], lw=ls2, lc=pcols[2], label=" P2c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_model[1:maxdep, 3], -zc[1:maxdep], lw=ls2, lc=pcols[3], label=" P3c", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_model[1:maxdep, 4], -zc[1:maxdep], lw=ls2, lc=pcols[4], label=" P4o", legendfontsize=lfs, alpha=ab)
        plot!(PP_model[1:maxdep, 5], -zc[1:maxdep], lw=ls2, lc=pcols[5], label=" P5o", legendfontsize=lfs, alpha=ab)
        plot!(PP_model[1:maxdep, 6], -zc[1:maxdep], lw=ls2, lc=pcols[6], label=" P6o", legendfontsize=lfs, alpha=ab)
        # plot!(sum(PP_model[1:20, :], dims=2), -zc[1:20], lw=ls2, lc="black", label=" Total PP", legendfontsize=lfs, yformatter=Returns(""),  
        # xlabel=L"mg~C/m3/day", xrotation=45, title =tls[5], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        # scatter!([PP_spot], [-5.0],  markersize=7, markercolor="green3", markerstrokecolor=:black, markershape=:diamond, 
        # label=" SPOT")
        # lens!([0, 1200], [-20, 0], inset=(1, bbox(0.03, 0.15, 0.45, 0.7, :top, :right)))

    season_num == 1 ? PP_spot = PP[1] : PP_spot = PP[3]

    p5 = plot(sum(PP_model[1:maxdep, :], dims=2), -zc[1:maxdep], lw=ls2, lc="black", label=" Total PP", legendfontsize=lfs, yformatter=Returns(""),  
        xlabel=L"mg~C/m^{3}/day", xrotation=45, title =tls[5], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab)
        scatter!([PP_spot], [-5.0],  markersize=7, markercolor="green3", markerstrokecolor=:black, markershape=:diamond, 
        label=" SPOT")
        lens!([800, 3000], [-60, 5], inset=(1, bbox(0.05, 0.43, 0.6, 0.3, :top, :right)), subplot=2, xrotation=45, 
        border=:box, bg_color_inside=:gray98)

    l = @layout [ a b c ; d e ] 
   

    f = plot(p1, p2, p3, p4, p5,
    layout = l,
    fg_legend = :transparent,
    size=(600,1000),
    plot_title = "SPOT vs. Model Productivity ($season)\n",
    )

    savefig(f,"$(dir)/$(filename).pdf")

    return f, dir, filename 

end


function bloom_productivity_plot(BP_pre, BP_blm, BP_spot_meso, BP_spot_oli, ds)

    H = 890
    zc = get_zc(H)

    tfs = 10 
    ls = 3
    ls2 = 7
    ab = 0.3
    lfs = 8
    ls3 = 3
    xtfs = 8
    lg = :bottomright

    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "mediumpurple3"]
    bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
    tls = ["\nModel BP (pre-bloom)", "\nModel BP (bloom)", "\nTotal BP", "\nModel PP","\nTotal PP"]

    lab_c = dropdims(sum(BP_pre[1:89, 4:6, :], dims=2), dims=2);
    lab_o = dropdims(sum(BP_pre[1:89, 9:11, :], dims=2), dims=2);
    lab_c_blm = dropdims(sum(BP_blm[1:89, 4:6, :], dims=2), dims=2);
    lab_o_blm = dropdims(sum(BP_blm[1:89, 9:11, :], dims=2), dims=2);

    spot_oli = dropmissing(BP_spot_oli, disallowmissing=true);
    spot_meso = dropmissing(BP_spot_meso, disallowmissing=true);


    p1 = plot(lab_c[:, :], -zc, lw=ls2, lc=bcols[4], label=" Lab C", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"cells/mL/day" ,
            xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab,
            foreground_color_legend = nothing, bottom_margin=5mm, top_margin=3mm)
            plot!(lab_o[:, :], -zc, lw=ls, lc=bcols[9], label=" Lab O", legendfontsize=lfs)
            plot!(BP_pre[:, 7], -zc, lw=ls2, lc=bcols[7], label=" S.Lab C", legendfontsize=lfs, alpha=ab)
            plot!(BP_pre[:, 12], -zc, lw=ls, lc=bcols[12], label=" S.Lab O", legendfontsize=lfs)
            plot!(BP_pre[:, 8], -zc, lw=ls2, lc=bcols[8], label=" Ref C", legendfontsize=lfs, alpha=ab)
            plot!(BP_pre[:, 13], -zc, lw=ls, lc=bcols[13], label=" Ref O", legendfontsize=lfs)
            lens!([100, 11000], [-450, -100], inset=(1, bbox(0.05, 0.18, 0.6, 0.3, :top, :right)), subplot=2, xrotation=45, 
            border=:box, bg_color_inside=:gray98, yticks=false)

    p2 = plot(lab_c_blm[:, :], -zc, lw=ls2, lc=bcols[4], label=" Lab C", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"cells/mL/day" ,
            xrotation=45, title =tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab,
            foreground_color_legend = nothing, bottom_margin=5mm, top_margin=3mm)
            plot!(lab_o_blm[:, :], -zc, lw=ls, lc=bcols[9], label=" Lab O", legendfontsize=lfs)
            plot!(BP_blm[:, 7], -zc, lw=ls2, lc=bcols[7], label=" S.Lab C", legendfontsize=lfs, alpha=ab)
            plot!(BP_blm[:, 12], -zc, lw=ls, lc=bcols[12], label=" S.Lab O", legendfontsize=lfs)
            plot!(BP_blm[:, 8], -zc, lw=ls2, lc=bcols[8], label=" Ref C", legendfontsize=lfs, alpha=ab)
            plot!(BP_blm[:, 13], -zc, lw=ls, lc=bcols[13],label=" Ref O", legendfontsize=lfs); 
            lens!([100, 11000], [-450, -100], inset=(1, bbox(0.05, 0.18, 0.6, 0.3, :top, :right)), subplot=2, xrotation=45, 
            border=:box, bg_color_inside=:gray98, yticks=false)

    p3 = plot(sum(BP_blm[:, 1:13], dims=2), -zc, lw=ls2, lc="cyan", label=" Model Meso", legendfontsize=lfs, yformatter=Returns(""),  xlabel=L"cells/mL/day" ,
            xrotation=45, title =tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, bottom_margin=5mm, top_margin=3mm)
            plot!(sum(BP_pre[:, 1:13], dims=2), -zc, lw=ls2, lc="purple", label=" Model Oli", legendfontsize=lfs, alpha=ab)
            scatter!(spot_meso.mean_Leu, -spot_meso.depth,  markersize=7, markercolor="cyan2", markerstrokecolor=:black, markershape=:diamond, 
            label=" SPOT Meso", foreground_color_legend = nothing, alpha=ab)
            scatter!(spot_oli.mean_Leu, -spot_oli.depth,  markersize=7, markercolor="purple2", markerstrokecolor=:black, markershape=:diamond, 
            label=" SPOT Oli", foreground_color_legend = nothing, alpha=ab)

    fig_BP = plot(p1, p2, p3, layout=(1,3), size=(700,400))

    savefig(fig_BP,"bloom_prod_BP.png")

    return fig_BP
end

function bloom_PP_plot(PP_pre, PP_blm, spot_meso, spot_oli, ds)

    H = 890
    zc = get_zc(H)

    tfs = 10 
    ls = 3
    ls2 = 7
    ab = 0.3
    lfs = 8
    ls3 = 3
    xtfs = 8
    lg = :bottomright

    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "mediumpurple3"]
    tls = ["\nModel PP (pre-bloom)", "\nModel PP (bloom)", "\nTotal PP"]

    maxdep=15
    p1 = plot(PP_pre[1:maxdep, 1], -zc[1:maxdep], lw=ls2, lc=pcols[1], label=" P1 C", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"mg~C/m^{3}/day" ,
        xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, linestyle=:dot, foreground_color_legend = nothing)
        plot!(PP_pre[1:maxdep, 2], -zc[1:maxdep], lw=ls2, lc=pcols[2], label=" P2 ", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_pre[1:maxdep, 3], -zc[1:maxdep], lw=ls2, lc=pcols[3], label=" P3 ", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_pre[1:maxdep, 4], -zc[1:maxdep], lw=ls2, lc=pcols[4], label=" P4 ", legendfontsize=lfs, alpha=ab)
        plot!(PP_pre[1:maxdep, 5], -zc[1:maxdep], lw=ls2, lc=pcols[5], label=" P5 ", legendfontsize=lfs, alpha=ab)
        plot!(PP_pre[1:maxdep, 6], -zc[1:maxdep], lw=ls2, lc=pcols[6], label=" P6 O", legendfontsize=lfs, alpha=ab)

    p2 = plot(PP_blm[1:maxdep, 1], -zc[1:maxdep], lw=ls2, lc=pcols[1], label=" P1 C", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"mg~C/m^{3}/day" ,
        xrotation=45, title =tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, linestyle=:dot, foreground_color_legend = nothing)
        plot!(PP_blm[1:maxdep, 2], -zc[1:maxdep], lw=ls2, lc=pcols[2], label=" P2 ", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_blm[1:maxdep, 3], -zc[1:maxdep], lw=ls2, lc=pcols[3], label=" P3 ", legendfontsize=lfs, alpha=ab, linestyle=:dot)
        plot!(PP_blm[1:maxdep, 4], -zc[1:maxdep], lw=ls2, lc=pcols[4], label=" P4 ", legendfontsize=lfs, alpha=ab)
        plot!(PP_blm[1:maxdep, 5], -zc[1:maxdep], lw=ls2, lc=pcols[5], label=" P5 ", legendfontsize=lfs, alpha=ab)
        plot!(PP_blm[1:maxdep, 6], -zc[1:maxdep], lw=ls2, lc=pcols[6], label=" P6 O", legendfontsize=lfs, alpha=ab)

    p3 = plot(sum(PP_blm[1:maxdep, :], dims=2), -zc[1:maxdep], lw=ls2, lc="green3", label=" Model Meso", legendfontsize=lfs, yformatter=Returns(""),  xlabel=L"mg~C/m^{3}/day",
        xrotation=45, title =tls[3], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab, bottom_margin=5mm, top_margin=3mm)
        plot!(sum(PP_pre[1:maxdep, :], dims=2), -zc[1:maxdep], lw=ls2, lc="black", label=" Model Oli", legendfontsize=lfs, alpha=ab)
        scatter!([spot_meso], [-5.0],  markersize=7, markercolor="green3", markerstrokecolor=:black, markershape=:diamond, 
        label=" SPOT Meso", foreground_color_legend = nothing, alpha=ab)
        scatter!([spot_oli], [-5.0],  markersize=7, markercolor="black", markerstrokecolor=:black, markershape=:diamond, 
        label=" SPOT Oli", foreground_color_legend = nothing, alpha=ab)


    fig_PP = plot(p1, p2, p3, layout=(1,3), size=(700,400))

    savefig(fig_PP,"bloom_prod_PP.png")

    return fig_PP

end