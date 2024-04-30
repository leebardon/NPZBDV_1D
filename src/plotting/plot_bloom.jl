# using NCDatasets
# using CairoMakie
# using LinearAlgebra, Statistics

include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBDV_1D/src/utils/save_utils.jl")

function plot_bloom_heatmaps(fsaven, varname, xunit, pulse_day)

    zc = get_zc(890);
    ds = NCDataset(fsaven);
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/blooms/" => "")
    p = ds["p"][:,:,:];
    b = ds["b"][:,:,:];
    d = ds["d"][:,:,:];
    z = ds["z"][:,:,:];
    n = ds["n"][:,:,:];

    # p = ds["p"][:,:,358680:end];
    # b = ds["b"][:,:,358680:end];
    # d = ds["d"][:,:,358680:end];
    # z = ds["z"][:,:,358680:end];
    # n = ds["n"][:,:,358680:end];

    # varname == "P" ? bloom_heatmaps(p, zc, filename, varname) : bloom_heatmaps_B_summed(b, zc, filename, xunit, pulse_day)

    # bloom_heatmaps_B_summed(b, zc, filename, xunit, pulse_day)
    bloom_heatmaps_B_RA(b, zc, filename, xunit, pulse_day)
    
    # bloom_heatmaps_BZ_OSM(b, z, zc, filename, xunit, pulse_day)
    bloom_heatmaps_DNZ_OSM(d, n, z, zc, filename, xunit, pulse_day)

end

function concat_outputs(fsaven1, fsaven2, state_var)

    ds1=NCDataset(fsaven1);
    ds2=NCDataset(fsaven2);
    statevar1 = ds1[state_var][:,:,:];
    statevar2 = ds2[state_var][:,:,:];

    return cat(statevar1, statevar2, dims=3)

end 

function get_hmap_data_ts(state_var, i, d, max_d)

    data = state_var[1:max_d, i, :]
    survivors = set_extinct_to_zero(data)

    ts = collect(1:size(survivors, 2))
    z = get_hmap_z_axis(d, ts, survivors)

    return z, ts

end

function get_hmap_data_days(state_var, i, d, max_d)

    data = state_var[1:max_d, i, :]
    survivors = set_extinct_to_zero(data)

    daily_data = survivors[:, 1:200:end]
    days = collect(1:size(daily_data, 2))
    z = get_hmap_z_axis(d, days, daily_data)

    return z, days

end

function bloom_heatmaps_B(state_var, zc, filename, xunit, pulse_day)

    parent_folder = "results/plots/blooms/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    tits = ["Labile POM Consumer", "POM Consumer", "Recalc. POM Consumer", "Fast-Lab Copio", "Lab Copio", "Semi-Lab Copio", "Semi-Rec Copio", "Rec Copio",
            "Fast-Lab Oligo", "Lab Oligo", "Semi-Lab Oligo", "Semi-Rec Oligo", "Rec Oligo"]

    max_d = 89
    d = -zc[1:max_d]
    num_species = size(state_var, 2)
    res1, res2 =(800, 350), (1500, 700)
    cmap1, cmap2 = :hawaii, :thermal

    f0 = Figure(size=res1);
    row_a, col_a = 1, 1
    for i in range(1, 3)
        xunit == "days" ? (z, days) = get_hmap_data_days(state_var, i, d, max_d) : (z, days) = get_hmap_data_ts(state_var, i, d, max_d)
        Axis(f0[row_a, col_a], xticklabelsize=11, yticklabelsize=11, title = tits[i], titlesize=11)
        hm = heatmap!(f0[row_a, col_a], days, d, z, colormap=cmap2, interpolate=true)
        # hm = heatmap!(f0[row_a, col_a], days, reverse(d), reverse(z), colormap=cmap2)
        Colorbar(f0[row_a, col_a+1], hm, ticklabelsize=11)
        col_a += 2
    end

    f1 = Figure(size=res2);
    row, col = 1, 1
    for i in range(4, num_species)
        xunit == "days" ? (z, days) = get_hmap_data_days(state_var, i, d, max_d) : (z, days) = get_hmap_data_ts(state_var, i, d, max_d)
        ax, hm = heatmap(f1[row, col], days, d, z, xlabelrotation=45, colormap=cmap2, interpolate=true)
        # ax, hm = heatmap(f1[row, col], days, reverse(d), reverse(z), xlabelrotation=45, colormap=cmap2)
        Colorbar(f1[row, col+1], hm)
        ax.title = tits[i]
        col += 2
        if col > 9 
            row += 1
            col = 1
        end
    end

    #---------------------------------------------------------------------------------------------------
    xunit == "days" ? xlab="Days" : xlab="Timesteps"
    ylab = "Depth (m)"

    f0[2, :] = Label(f0, xlab, fontsize=12)
    f0[:, 0] = Label(f0, ylab, rotation = pi/2, fontsize=12)
    f0[:, 7] = Label(f0, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=12)
    f0[0, :] = Label(f0, "Response of POM Consumers to Nutrient Pulse at $(pulse_day) days", fontsize=18)

    f1[0, :] = Label(f1, "Response of DOM Consumers to Nutrient Pulse at $(pulse_day) days", fontsize=35)
    f1[1:end, 0] = Label(f1, ylab, rotation = pi/2, fontsize=20)
    f1[1:end, 11] = Label(f1, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=20)
    f1[3, :] = Label(f1, xlab, fontsize=20)

    save("$(dir)/$(filename)_pom.png", f0)
    save("$(dir)/$(filename)_dom.png", f1)
    # save("_pom2.png", f0)
    # save("dom_1000days_bloom.png", f1)
end

function bloom_heatmaps_B_summed(state_var, zc, filename, xunit, pulse_day)

    parent_folder = "results/plots/blooms/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    tits = ["Labile Copio", "Semi-Labile Copio", "Refractory Copio","Labile Oligo", "Semi-Labile Oligo", "Refractory Oligo" ]
    
    max_d = 60;
    d = -zc[1:max_d];
    res1, res2 =(800, 350), (900, 700)
    cmap1, cmap2 = :solar, :thermal

    #---------------------------------------------------------------------------------------------------

    # f0 = Figure(size=res1);
    # row_a, col_a = 1, 1
    # for i in range(1, 3)
    #     xunit == "days" ? (z, days) = get_hmap_data_days(state_var, i, d, max_d) : (z, days) = get_hmap_data_ts(state_var, i, d, max_d)
    #     Axis(f0[row_a, col_a], xticklabelsize=11, yticklabelsize=11, title = tits[i], titlesize=11)
    #     hm = heatmap!(f0[row_a, col_a], days, d, z, colormap=cmap2, interpolate=true)
    #     Colorbar(f0[row_a, col_a+1], hm, ticklabelsize=11)
    #     col_a += 2
    # end

    #---------------------------------------------------------------------------------------------------
    sum_lab_cop = dropdims(sum(state_var[1:max_d, 4:6, :], dims=2), dims=2);
    sum_lab_oli = dropdims(sum(state_var[1:max_d, 9:11, :], dims=2), dims=2);

    st=1;
    # step=20;   # for steady state 'regular' runs
    step=200 # for higher-res bloom runs

    zmax, z97 = set_zmax(state_var[:,4:end,:], size(state_var[:,4:end,:],2))
    joint_limits = (0, zmax)

    f1 = Figure(size=res2);
    data1 = sum_lab_cop;
    survivors1 = set_extinct_to_zero(data1);
    daily_data1 = survivors1[:, st:step:end];
    days = collect(1:size(daily_data1, 2));
    z1 = get_hmap_z_axis(d, days, daily_data1);
    ax1, hm1 = heatmap(f1[1, 1], days, d, z1, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 2], hm1)
    ax1.title = tits[1]

    data2 = state_var[1:max_d, 7, :];
    survivors2 = set_extinct_to_zero(data2);
    daily_data2 = survivors2[:, st:step:end];
    z2 = get_hmap_z_axis(d, days, daily_data2);
    ax2, hm2 = heatmap(f1[1, 2], days, d, z2, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 4], hm2)
    ax2.title = tits[2]

    data3 = state_var[1:max_d, 8, :];
    survivors3 = set_extinct_to_zero(data3);
    daily_data3 = survivors3[:, st:step:end];
    z3 = get_hmap_z_axis(d, days, daily_data3);
    ax3, hm3 = heatmap(f1[1, 3], days, d, z3, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 6], hm3)
    ax3.title = tits[3]

    data4 = sum_lab_oli;
    survivors4 = set_extinct_to_zero(data4);
    daily_data4 = survivors4[:, st:step:end];
    z4 = get_hmap_z_axis(d, days, daily_data4);
    ax4, hm4 = heatmap(f1[2, 1], days, d, z4, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[2, 2], hm4)
    ax4.title = tits[4]

    data5 = state_var[1:max_d, 12, :];
    survivors5 = set_extinct_to_zero(data5);
    daily_data5 = survivors5[:, st:step:end];
    z5 = get_hmap_z_axis(d, days, daily_data5);
    ax5, hm5 = heatmap(f1[2, 2], days, d, z5, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[2, 4], hm5)
    ax5.title = tits[5]

    data6 = state_var[1:max_d, 13, :];
    survivors6 = set_extinct_to_zero(data6);
    daily_data6 = survivors6[:, st:step:end];
    z6 = get_hmap_z_axis(d, days, daily_data6);
    ax6, hm6 = heatmap(f1[2, 3], days, d, z6, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[2, 6], hm6)
    ax6.title = tits[6]


    Colorbar(f1[:, end+1], colormap=cmap2, colorrange=joint_limits, size=30, label=rich("mmol N/m", superscript("3")), labelsize=20)

    #---------------------------------------------------------------------------------------------------

    xunit == "days" ? xlab="Days" : xlab="Timesteps"
    ylab = "Depth (m)"

    # f0[2, :] = Label(f0, xlab, fontsize=12)
    # f0[:, 0] = Label(f0, ylab, rotation = pi/2, fontsize=12)
    # f0[:, 7] = Label(f0, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=12)
    # f0[0, :] = Label(f0, "Response of POM Consumers to Nutrient Pulse at $(pulse_day) days", fontsize=18)

    f1[0, :] = Label(f1, "Response of DOM Consumers to Simulated Bloom", fontsize=25)
    f1[1:end, 0] = Label(f1, ylab, rotation = pi/2, fontsize=20)
    # f1[1:end, 7] = Label(f1, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=20)
    f1[3, :] = Label(f1, xlab, fontsize=20)

    save("testbloom.png", f1)

end

function bloom_heatmaps_BZ_OSM(state_var, Z, zc, filename, xunit, pulse_day)

    parent_folder = "results/plots/blooms/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    tits = ["Labile Oligo", "Semi-Labile Copio", "Refractory Copio", "Picohet. Grazer" ]
    
    max_d = 60;
    d = -zc[1:max_d];
    res1, res2 =(1000, 350), (900, 700)
    cmap1, cmap2 = :hawaii, :thermal

    #---------------------------------------------------------------------------------------------------

    f0 = Figure(size=res1);
    row_a, col_a = 1, 1
    for i in range(1, 3)
        xunit == "days" ? (z, days) = get_hmap_data_days(state_var, i, d, max_d) : (z, days) = get_hmap_data_ts(state_var, i, d, max_d)
        Axis(f0[row_a, col_a], xticklabelsize=11, yticklabelsize=11, title = tits[i], titlesize=11)
        hm = heatmap!(f0[row_a, col_a], days, d, z, colormap=cmap2, interpolate=true)
        Colorbar(f0[row_a, col_a+1], hm, ticklabelsize=11)
        col_a += 2
    end

    #---------------------------------------------------------------------------------------------------

    sum_lab_cop = dropdims(sum(state_var[1:max_d, 4:6, :], dims=2), dims=2);
    sum_lab_oli = dropdims(sum(state_var[1:max_d, 9:11, :], dims=2), dims=2);


    f1 = Figure(size=res1);
    data1 = sum_lab_oli;
    survivors1 = set_extinct_to_zero(data1);
    daily_data1 = survivors1[:, 1:200:end];
    days = collect(1:size(daily_data1, 2));
    z1 = get_hmap_z_axis(d, days, daily_data1);
    ax1, hm1 = heatmap(f1[1, 1], days, d, z1, xlabelrotation=45, colormap=cmap2, interpolate=true)
    Colorbar(f1[1, 2], hm1)
    ax1.title = tits[1]

    data2 = state_var[1:max_d, 7, :];
    survivors2 = set_extinct_to_zero(data2);
    daily_data2 = survivors2[:, 1:200:end];
    z2 = get_hmap_z_axis(d, days, daily_data2);
    ax2, hm2 = heatmap(f1[1, 3], days, d, z2, xlabelrotation=45, colormap=cmap2, interpolate=true)
    Colorbar(f1[1, 4], hm2)
    ax2.title = tits[2]

    data3 = state_var[1:max_d, 8, :];
    survivors3 = set_extinct_to_zero(data3);
    daily_data3 = survivors3[:, 1:200:end];
    z3 = get_hmap_z_axis(d, days, daily_data3);
    ax3, hm3 = heatmap(f1[1, 5], days, d, z3, xlabelrotation=45, colormap=cmap2, interpolate=true)
    Colorbar(f1[1, 6], hm3)
    ax3.title = tits[3]

    data4 = Z[1:max_d, 3, :];
    daily_data4 = data4[:, 1:200:end];
    z4 = get_hmap_z_axis(d, days, daily_data4);
    ax4, hm4 = heatmap(f1[1, 7], days, d, z4, xlabelrotation=45, colormap=cmap2, interpolate=true)
    Colorbar(f1[1, 8], hm4)
    ax4.title = tits[4]

    #---------------------------------------------------------------------------------------------------

    xunit == "days" ? xlab="Days" : xlab="Timesteps"
    ylab = "Depth (m)"

    f0[2, :] = Label(f0, xlab, fontsize=12)
    f0[:, 0] = Label(f0, ylab, rotation = pi/2, fontsize=12)
    f0[:, 7] = Label(f0, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=12)
    f0[0, :] = Label(f0, "Response of POM Consumers to Nutrient Pulse at $(pulse_day) days", fontsize=18)

    f1[0, :] = Label(f1, "Response of DOM Consumers to Simulated Bloom", fontsize=18)
    f1[:, 0] = Label(f1, ylab, rotation = pi/2, fontsize=14)
    f1[:, 9] = Label(f1, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=14)
    f1[2, :] = Label(f1, xlab, fontsize=14)

    save("bz_bloom_4_reduced.png", f1)
end

function bloom_heatmaps_DNZ_OSM(D, N, Z, zc, filename, xunit, pulse_day)

    parent_folder = "results/plots/blooms/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    tits = ["Labile DOM", "Semi-Labile DOM", "Refractory DOM", "Inorganic Nutrients", "Zooplank. Conc."]
    
    max_d = 60;
    d = -zc[1:max_d];
    res1 =(1300, 400);
    cmap2 =:thermal

    #---------------------------------------------------------------------------------------------------
    start=1;
    # step=200; stop=37001 # for higher-res bloom runs, plots to 180days
    step=200; stop=73001 # plots to 365days

    f1 = Figure(size=res1);
    zmax, z97 = set_zmax(D[:,4:8,:], size(D[:,4:8,:],2));
    joint_limits = (0, zmax);

    sum_lab_dom = dropdims(sum(D[1:max_d, 4:6, :], dims=2), dims=2);
    data1 = sum_lab_dom;
    daily_data1 = data1[:, 1:step:stop];
    days = collect(1:size(daily_data1, 2));
    z1 = get_hmap_z_axis(d, days, daily_data1);
    ax1, hm1 = heatmap(f1[1, 1], days, d, z1, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 2], hm1)
    ax1.title = tits[1]

    data2 = D[1:max_d, 7, :];
    daily_data2 = data2[:, 1:step:stop];
    z2 = get_hmap_z_axis(d, days, daily_data2);
    ax2, hm2 = heatmap(f1[1, 2], days, d, z2, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 4], hm2)
    ax2.title = tits[2]

    data3 = D[1:max_d, 8, :];
    daily_data3 = data3[:, 1:step:stop];
    z3 = get_hmap_z_axis(d, days, daily_data3);
    ax3, hm3 = heatmap(f1[1, 3], days, d, z3, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 6], hm3)
    ax3.title = tits[3]

    data4 = dropdims(N[1:max_d, :, :], dims=2);
    daily_data4 = data4[:, 1:step:stop];
    z4 = get_hmap_z_axis(d, days, daily_data4);
    ax4, hm4 = heatmap(f1[1, 4], days, d, z4, xlabelrotation=45, colormap=cmap2, interpolate=true)
    ax4.title = tits[4]

    data5 = Z[1:max_d, 3, :];
    daily_data5 = data5[:, 1:step:stop];
    z5 = get_hmap_z_axis(d, days, daily_data5);
    ax5, hm5 = heatmap(f1[1, 5], days, d, z5, xlabelrotation=45, colormap=cmap2, interpolate=true)
    ax5.title = tits[5]

    Colorbar(f1[:, end+1], colormap=cmap2, colorrange=joint_limits, size=30, label=rich("mmol N/m", superscript("3")), labelsize=20)

    #---------------------------------------------------------------------------------------------------

    xunit == "days" ? xlab="Days" : xlab="Timesteps"
    ylab = "Depth (m)"

    f1[0, :] = Label(f1, "DOM & Nutrient Concentration", fontsize=18)
    f1[:, 0] = Label(f1, ylab, rotation = pi/2, fontsize=14)
    # f1[:, 5] = Label(f1, rich("mmol N/m", superscript("3")), rotation = pi/2, fontsize=14)
    f1[2, :] = Label(f1, xlab, fontsize=14)

    save("test7_1816_DNZ.png", f1)
end


# ----------------------------------------------------------------------------------------------------------

function calculate_bdom_RA(bdom)

    tsteps = size(bdom)[3]
    nb = size(bdom)[2]
    nboxes = size(bdom)[1]

    total_biomass = Array{Float64, 2}(undef, nboxes, tsteps);
    RA = Array{Float64, 3}(undef, nboxes, nb, tsteps);

    for i in range(1,tsteps)
        total_biomass[:,i]=sum(bdom[:,:,i], dims=2)
    end

    for i in range(1,nb)
        for j in range(1,tsteps)
            RA[:,i,j] = bdom[:,i,j] ./ total_biomass[:,j]
        end
    end

    return RA

end

function bloom_heatmaps_B_RA(state_var, zc, filename, xunit, pulse_day)

    parent_folder = "results/plots/blooms/heatmaps/"
    dir = check_subfolder_exists(filename, parent_folder)

    tits = ["Labile Copio", "Semi-Labile Copio", "Refractory Copio","Labile Oligo", "Semi-Labile Oligo", "Refractory Oligo" ]
    
    max_d = 60;
    d = -zc[1:max_d];
    res1, res2 =(800, 350), (900, 700)
    cmap1, cmap2 = :solar, :thermal

    #---------------------------------------------------------------------------------------------------
    bdom = state_var[:,4:end,:];
    RA = calculate_bdom_RA(bdom);
    sum_lab_cop_RA = dropdims(sum(RA[1:max_d, 1:3, :], dims=2), dims=2);
    sum_lab_oli_RA = dropdims(sum(RA[1:max_d, 6:8, :], dims=2), dims=2);

    start=1;
    # step=200; stop=37001 # plots to 180days
    step=200; stop=73001 # plots to 365days
    # step=200; stop=146001 # plots to 730days

    f1 = Figure(size=res2);
    zmax, z97 = set_zmax(RA, size(RA,2));
    joint_limits = (0, zmax);

    data1 = sum_lab_cop_RA;
    daily_data1 = data1[:, start:step:stop];
    days = collect(1:size(daily_data1, 2));
    z1 = get_hmap_z_axis(d, days, daily_data1);
    ax1, hm1 = heatmap(f1[1, 1], days, d, z1, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits)
    # Colorbar(f1[1, 2], hm1)
    ax1.title = tits[1]

    data2 = RA[1:max_d, 4, :];
    daily_data2 = data2[:, start:step:stop];
    z2 = get_hmap_z_axis(d, days, daily_data2);
    ax2, hm2 = heatmap(f1[1, 2], days, d, z2, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits, yformatter=Returns(""))
    # Colorbar(f1[1, 4], hm2)
    ax2.title = tits[2]

    data3 = RA[1:max_d, 5, :];
    daily_data3 = data3[:, start:step:stop];
    z3 = get_hmap_z_axis(d, days, daily_data3);
    ax3, hm3 = heatmap(f1[1, 3], days, d, z3, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits, yformatter=Returns(""))
    # Colorbar(f1[1, 6], hm3)
    ax3.title = tits[3]

    data4 = sum_lab_oli_RA;
    daily_data4 = data4[:, start:step:stop];
    z4 = get_hmap_z_axis(d, days, daily_data4);
    ax4, hm4 = heatmap(f1[2, 1], days, d, z4, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits, yformatter=Returns(""))
    # Colorbar(f1[2, 2], hm4)
    ax4.title = tits[4]

    data5 = RA[1:max_d, 9, :];
    daily_data5 = data5[:, start:step:stop];
    z5 = get_hmap_z_axis(d, days, daily_data5);
    ax5, hm5 = heatmap(f1[2, 2], days, d, z5, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits, yformatter=Returns(""))
    # Colorbar(f1[2, 4], hm5)
    ax5.title = tits[5]

    data6 = RA[1:max_d, 10, :];
    daily_data6 = data6[:, start:step:stop];
    z6 = get_hmap_z_axis(d, days, daily_data6);
    ax6, hm6 = heatmap(f1[2, 3], days, d, z6, xlabelrotation=45, colormap=cmap2, interpolate=true, colorrange=joint_limits, yformatter=Returns(""))
    # Colorbar(f1[2, 6], hm6)
    ax6.title = tits[6]

    Colorbar(f1[:, end+1], colormap=cmap2, colorrange=joint_limits, size=30, label="RA", labelsize=20)

    #---------------------------------------------------------------------------------------------------

    xunit == "days" ? xlab="Days" : xlab="Timesteps"
    ylab = "Depth (m)"

    f1[0, :] = Label(f1, "Relative Abundance over Time", fontsize=25)
    f1[1:end, 0] = Label(f1, ylab, rotation = pi/2, fontsize=20)
    # f1[1:end, 7] = Label(f1, rich("Relative Abundance"), rotation = pi/2, fontsize=20)
    f1[3, :] = Label(f1, xlab, fontsize=20)

    save("test7_1816.png", f1)

end

pulse_day = 90

# fsaven = "results/outfiles/blooms/blm_cont_240204_01:43_Wi50yPP_6P4Z13B8D.nc"
# fsaven = "results/outfiles/blooms/blm_240205_16:56_Wi50yPP_6P4Z13B8D.nc"
# fsaven = "results/outfiles/blooms/blm_240205_18:00_Wi50yPP_6P4Z13B8D.nc"
# fsaven="results/outfiles/cont_blm_240204_18:02_Wi100yNP_6P4Z13B8D.nc"
# fsaven="results/outfiles/blooms/blm_240212_17:19_Wi30yNP_6P3Z13B8D.nc";
# fsaven="results/outfiles/bloom_tests/mlz2_blm_240213_18:16_Wi50yNP_6P3Z13B8D.nc"
# fsaven="results/outfiles/bloom_tests/mlz_blm_240213_18:16_Wi50yNP_6P3Z13B8D.nc"
# fsaven="results/outfiles/bloom_tests/unsteady_blm_240213_18:16_Wi50yNP_6P3Z13B8D.nc"
# fsaven="results/outfiles/blooms/blm_240213_18:16_Wi50yNP_6P3Z13B8D.nc"
# plot_bloom_heatmaps(fsaven, "B", "days", pulse_day)


# function BP_plt(BP_pre, BP_blm, ds)

#     H = 890
#     zc = get_zc(H)

#     tfs = 9 
#     ls = 3
#     ls2 = 7
#     ab = 0.3
#     lfs = 7
#     ls3 = 3
#     xtfs = 8
#     lg = :bottomright

#     pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "mediumpurple3"]
#     bcols = ["teal", "azure4", "red4", "black", "seagreen", "purple4", "maroon", "brown3", "grey", "lime", "orchid", "pink2", "coral"]
#     tls = ["\nPicohet. Prod. (pre-bloom)", "\nPicohet. Prod. (bloom)", "\nTotal BP", "\nPhytoplankton (model)","\nTotal PP"]

#     lab_c = dropdims(sum(BP_pre[1:89, 4:6, :], dims=2), dims=2);
#     lab_o = dropdims(sum(BP_pre[1:89, 9:11, :], dims=2), dims=2);
#     lab_c_blm = dropdims(sum(BP_blm[1:89, 4:6, :], dims=2), dims=2);
#     lab_o_blm = dropdims(sum(BP_blm[1:89, 9:11, :], dims=2), dims=2);


#     p1 = plot(lab_c[:, :], -zc, lw=ls2, lc=bcols[4], label=" LC", legendfontsize=lfs, ylabel="Depth (m)", xlabel=L"cells/mL/day" ,
#             xrotation=45, title =tls[1], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab,
#             foreground_color_legend = nothing)
#             plot!(lab_o[:, :], -zc, lw=ls, lc=bcols[9], label=" LO", legendfontsize=lfs)
#             plot!(BP_pre[:, 7], -zc, lw=ls2, lc=bcols[7], label=" SLC", legendfontsize=lfs, alpha=ab)
#             plot!(BP_pre[:, 12], -zc, lw=ls, lc=bcols[12], label=" SLO", legendfontsize=lfs)
#             plot!(BP_pre[:, 8], -zc, lw=ls2, lc=bcols[8], label=" RC", legendfontsize=lfs, alpha=ab)
#             plot!(BP_pre[:, 13], -zc, lw=ls, lc=bcols[13], label=" RO", legendfontsize=lfs)
#             lens!([100, 11000], [-450, -100], inset=(1, bbox(0.05, 0.22, 0.6, 0.3, :top, :right)), subplot=2, xrotation=45, 
#             border=:box, bg_color_inside=:gray98, yticks=false)

#     p2 = plot(lab_c_blm[:, :], -zc, lw=ls2, lc=bcols[4], label=" LC", legendfontsize=lfs, yformatter=Returns(""), xlabel=L"cells/mL/day" ,
#             xrotation=45, title =tls[2], titlefontsize=tfs, grid=false, border=:box, legend=lg, xtickfontsize=xtfs, alpha=ab,
#             foreground_color_legend = nothing)
#             plot!(lab_o_blm[:, :], -zc, lw=ls, lc=bcols[9], label=" LO", legendfontsize=lfs)
#             plot!(BP_blm[:, 7], -zc, lw=ls2, lc=bcols[7], label=" SLC", legendfontsize=lfs, alpha=ab)
#             plot!(BP_blm[:, 12], -zc, lw=ls, lc=bcols[12], label=" SLO", legendfontsize=lfs)
#             plot!(BP_blm[:, 8], -zc, lw=ls2, lc=bcols[8], label=" RC", legendfontsize=lfs, alpha=ab)
#             plot!(BP_blm[:, 13], -zc, lw=ls, lc=bcols[13],label=" RO", legendfontsize=lfs); 
#             lens!([100, 11000], [-450, -100], inset=(1, bbox(0.05, 0.22, 0.6, 0.3, :top, :right)), subplot=2, xrotation=45, 
#             border=:box, bg_color_inside=:gray98, yticks=false)

#     fig = plot(p1, p2, layout=(1,2))

#     savefig(fig,"bloom_prod.png")

#     return fig
# end
