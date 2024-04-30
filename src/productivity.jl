using NCDatasets, SparseArrays, LinearAlgebra, Statistics
using CSV, DataFrames, LaTeXStrings, Plots, Measures

include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")
include("/home/lee/Dropbox/Development/NPZBD_1D/src/plotting/productivity_plots.jl")


function get_productivity(fsaven, bloom=0)

    ds = NCDataset(fsaven)
    ds.attrib["Season"] == "winter" ? season_num = 1 : season_num = 2

    if ds["pulse"][:] == 1
        N, P, Z, B, D = get_endpoints(["n", "p", "z", "b", "d"], ds)

    elseif bloom == 1
        N, P, Z, B, D = load_bloom(["n", "p", "z", "b", "d"], ds)
        SPOT_BP_meso, SPOT_BP_oli, SPOT_PP_meso, SPOT_PP_oli = get_productivity_spot("dummy", bloom)
        BP_pre, BP_blm = bloom_BP(B, D, ds)
        PP_pre, PP_blm = bloom_PP(P, N, ds)
        fig1 = bloom_productivity_plot(BP_pre, BP_blm, SPOT_BP_meso, SPOT_BP_oli, ds)
        fig2 = bloom_PP_plot(PP_pre, PP_blm, SPOT_PP_meso, SPOT_PP_oli, ds)

        return

    else
        N, P, Z, B, D = mean_over_time(["n", "p", "z", "b", "d"], ds, season_num)
    end

    BP_Leu, BP_Thy, PP = get_productivity_spot(ds.attrib["Season"])
    BP_model = bacterial_productivity(B, D, ds)
    PP_model = primary_productivity(P, N, ds)

    fig, dir, filename = plot_seasonal_BP(BP_model, PP_model, BP_Leu, BP_Thy, PP, ds, fsaven, season_num)
    savefig(fig,"$(dir)/$(filename).pdf")

end

function bacterial_productivity(B, D, ds)

    B[B .< 1e-6] .= 0
    D[D .< 1e-6] .= 0

    CM = ds["CM"][:,:]
    umax_ij = ds["umax_ij"][:,:]
    Km_ij = ds["Km_ij"][:,:]
    temp = ds["temp_fun"][:]

    II, JJ = get_nonzero_axes(CM)
    BP = zeros(89, 13)

    for j = axes(II, 1)
        mu = temp .* umax_ij[II[j],JJ[j]] .* D[:,II[j]] ./ (D[:,II[j]] .+ Km_ij[II[j],JJ[j]])
        BP[:,JJ[j]] = B[:,JJ[j]] .* mu
    end

    BP_out = convert_units(BP)

    return BP_out

end

function primary_productivity(P, N, ds)

    P[P .< 1e-6] .= 0
    N[N .< 1e-6] .= 0

    CMp = ds["CMp"][:,:]
    vmax_ij = ds["vmax_ij"][:,:]
    Kp_ij = ds["Kp_ij"][:,:]
    temp = ds["temp_fun"][:]
    K_I = ds["K_I"][:]

    Iz = light_atten(P)
    II, JJ = get_nonzero_axes(CMp)
    PP = zeros(89, 6)

    for j = axes(II, 1)
        mu = temp .* vmax_ij[II[j],JJ[j]] .* min.(N ./ (N .+ Kp_ij[II[j],JJ[j]]), Iz ./ (Iz .+ K_I))
        PP[:,JJ[j]] = P[:,JJ[j]] .* mu
    end

    PP_out = convert_units_PP(PP)

    return PP_out

end

function bloom_BP(B, D, ds)
    # assumes bloom tracked over 365 days

    prebloom = B[:, :, 1:17999];
    prebloom_D = D[:, :, 1:17999];
    bloom = B[:,:,18000:end];
    bloom_D = D[:, :, 18000:end];

    prebloom_mean = dropdims(mean(prebloom, dims=3), dims=3);
    prebloom_D_mean = dropdims(mean(prebloom_D, dims=3), dims=3);
    bloom_mean = dropdims(mean(bloom, dims=3), dims=3);
    bloom_D_mean = dropdims(mean(bloom_D, dims=3), dims=3);

    BP_prebloom = bacterial_productivity(prebloom_mean, prebloom_D_mean, ds);
    BP_bloom = bacterial_productivity(bloom_mean, bloom_D_mean, ds);

    return BP_prebloom, BP_bloom

end   

function bloom_PP(P, N, ds)
    # assumes bloom tracked over 365 days

    prebloom = P[:, :, 1:17999];
    prebloom_N = N[:, :, 1:17999];
    bloom = P[:,:,18000:end];
    bloom_N = N[:, :, 18000:end];

    prebloom_mean = dropdims(mean(prebloom, dims=3), dims=3);
    prebloom_N_mean = dropdims(mean(prebloom_N, dims=3), dims=3);
    bloom_mean = dropdims(mean(bloom, dims=3), dims=3);
    bloom_N_mean = dropdims(mean(bloom_N, dims=3), dims=3);

    PP_prebloom = primary_productivity(prebloom_mean, prebloom_N_mean, ds);
    PP_bloom = primary_productivity(bloom_mean, bloom_N_mean, ds);

    return PP_prebloom, PP_bloom

end  

function light_atten(P)
    # Following Zakem et al 2015 (using chl2_min as default)

    I_max = 1400                                        # Incident radiation at surface W/m2 #TODO what's the val at SPOT?
    I_in = I_max/2                                      # avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
    a_chlD = 0.04                                       # Absorption coeff incorporating Chl-a and CDOM (m2/mg Chl) 
    # chl2c_max = 0.2                                       # max chlorophyll to carbon ratio (mg Chl/mmol C)
    chl2c_min = 0.02                                    # min chlorophyll to carbon ratio (mg Chl/mmol C)
    kw = 0.04                                           # attenuation coeff of water (m2/mg Chl)
    
    zc, ngrid = get_zc(890), 89                         # height of water column (m), num boxes

    Iz = zeros(ngrid)                         
    chl_tot = sum(P, dims=2) .* chl2c_min .* 6.6        # mmolN/m3 * mgChl/mmolC  (redfield ratio -> 6.6 C for every N)

    for d in range(1, ngrid)
        Iz[d] = I_in * exp(-zc[d]*(kw + sum(chl_tot[1:d]*a_chlD)))
    end

    return Iz

end

function get_productivity_spot(season,bloom=0)

    csv_path = "data/spot_data/seasonal_means/bottle/seasonal_means.csv"
    df = DataFrame(CSV.File(csv_path))

    # Winter, spring, summer, fall
    mean_seasonal_PP_spot = [1018.022752636842, 1903.5313703434347, 1179.4424693131705, 892.4915176739885]

    if bloom==0
        BP_leu, BP_thy = select_cols(df)
        mean_leu, mean_thy = filter_by_season(BP_leu, BP_thy, season)

        return mean_leu, mean_thy, mean_seasonal_PP_spot
    else
        BP_leu, BP_thy = select_cols(df)
        meanBP_meso, meanBP_oli = filter_by_season_bloom(BP_leu, "spring", "fall")
        PP_meso, PP_oli = mean_seasonal_PP_spot[2], mean_seasonal_PP_spot[4]

        return meanBP_meso, meanBP_oli, PP_meso, PP_oli
    end

end

function select_cols(df)

    BP_leu = select(df[:,:], :depth, :mean_Leu, :sd_Leu, :season)
    BP_thy = select(df[:,:], :depth, :mean_Thy, :sd_Thy, :season)

    return BP_leu, BP_thy
end

function filter_by_season(leu_df, thy_df, seson)

    leu = filter(row -> row.season == seson, leu_df)
    thy = filter(row -> row.season == seson, thy_df)

    return leu, thy

end 

function filter_by_season_bloom(leu_df, spring, fall)

    meanBP_meso = filter(row -> row.season == spring, leu_df)
    meanBP_oli = filter(row -> row.season == fall, leu_df)

    return meanBP_meso, meanBP_oli

end 

function convert_units(BP)

    # From - mmol N/m3 to cells/mL
    # Estimate of cell nitrogen quota -> 1 fmol N/cell = 1e^-12 mmol N/cell
    # 1 mmol N/1e^-12 cells = 1e^12 cells per mol N ; 1 m3 = 1e^6 mL
    # 1 mmol N/m3 == 1e^12 / 1e^6 = 1e^6 cells/mL

    return BP .*= 1e6

end

function convert_units_PP(PP)

    # From - mmol N/m3 to mg C/m2
    # (mg = mmol * atomic weight) -- 1 mmol N/m3 = 14 mg N/m3
    # 1 mmol N/m3 = (14 * 6.6) mg C/m3
    # --> 1 mmol N/m3 = 92.4 mg C/m3
    #NOTE cant directly convert volume to area - however, satellite optical depth is approx 0-10m 
    # according to e.g. https://www.mdpi.com/2072-4292/9/3/301 so surface box of model may be considered
    # comparable to PP estimates from satellite

    return PP .*= 920.4

end



# fsaven = "results/outfiles/Wi50y_240103_14:18_6P3Z13B8D.nc" # win steady
# fsaven = "results/outfiles/Wi50y_240108_20:47_6P3Z13B8D.nc" # win pulsed
# fsaven="results/outfiles/Su50y_240104_14:04_6P3Z13B8D.nc" # sum steady
# fsaven="results/outfiles/Su50y_240109_00:03_6P3Z13B8D.nc" #sum pulsed

#changed z mq to 0.05
# fsaven = "results/outfiles/Wi50ySP_240118_13:16_6P3Z13B8D.nc" # win pulsed
# fsaven="results/outfiles/Su50ySP_240118_15:39_6P3Z13B8D.nc" #sum pulsed

# fsaven = "results/outfiles/240124_17:58_Wi50yPP_6P4Z13B8D.nc" # win pulsed
# fsaven="results/outfiles/240124_19:30_Su50yPP_6P4Z13B8D.nc" #sum pulsed

# get_productivity(fsaven)
fsaven="results/outfiles/blooms/blm_240213_18:16_Wi50yNP_6P3Z13B8D.nc"
bloom=1
get_productivity(fsaven, bloom)


