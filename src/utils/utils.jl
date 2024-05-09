
using DataFrames, NCDatasets, Printf

#TODO create savepath that includes pulse type

function message(v::String, nd::Int64=0, nb::Int64=0, nn::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "ST1" => ["Start Prescribed Model Run", "Continue Run", "Track Bloom"],
        "ST2" => "\nChoose run type: ",
        "CT2" => "\nContinue primary run or bloom run?",
        "CT1" => ["Primary run", "Bloom run"],
        "TM1" => ["bloom (365 days)", "2 years (days=732)", "10 years (days=3660)", "30 years (days=10980)", "50 years (days=18300)", "100 years (days=36600)"],
        "TM2" => "\nSelect Simulation Runtime:",
        "P1" => ["None (steady state)", "Periodic pulse at 25 (winter) or 45 (summer) day intervals", "Semi-stochastic pulse at 25 or 40 day intervals"],
        "P2" => "Select nutrient pulsing regime: ",
        "ENV" => "\n SETTING NUTRIENT SUPPLY \n -------------------------- ",
        "SE2" => "\nSimulate winter or summer conditions?",
        "SE1" => ["Winter", "Summer"],
        "GZ2" => "\nInclude grazers?",
        "GZ1" => ["yes", "no"],
        "LY2" => "\nInclude explicit or implicit viral lysis?",
        "LY1" => ["Explicit", "Implicit"],
        "LVP" => "\nSelect params to load: ",
        "SV" => "Saving to: $fsaven",

    )

    return m["$v"]

end


function microbe_num(MSG)
    println(message(MSG))
    input = readline()
    n = parse(Int64, input)
    return n
end


function user_select(run_type=0)


    println(message("DN"))
    input = readline()
    nd = parse(Int64, input) 
    nb = microbe_num("BN")
    np = microbe_num("PN")
    nz = microbe_num("ZN")
    nv = microbe_num("VN")
    nn = 1
    println(message("SUB"))
    uptake = request(message("UP2"), RadioMenu(message("UP1")))
    uptake == 1 ? umax_i = ordered_uptake_arr(nd) : umax_i = random_uptake_arr(nd)
    uptake_p = request(message("UPP2"), RadioMenu(message("UPP1")))
    uptake_p == 1 ? vmax_i = fill(1., np) : vmax_i = random_uptake_arr(np)

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    yield == 1 ? y_i = ones(nd)*0.3 : y_i = rand(nd)*0.5
    supply_weight = request(message("SW2"), RadioMenu(message("SW1")))

    println(message("ENV"))
    pulse = request(message("P2"), RadioMenu(message("P1")))
    season = request(message("SE2"), RadioMenu(message("SE1")))

    @info("User Selections: \n pulse type = $pulse, SW = $supply_weight \n B yield = $y_i \n B uptake = $umax_i \n P uptake = $vmax_i \n Season == $season \n")

    return nd, nb, np, nz, nn, nv, y_i, supply_weight, umax_i, vmax_i, season, pulse

end


function get_previous_params(continue_type=1)

    if continue_type == 1
        files = readdir("results/outfiles/endpoints")
        f = request("\nSelect output file:", RadioMenu(files))
        ds = NCDataset("results/outfiles/endpoints/$(files[f])")
    else
        files = readdir("results/outfiles/blooms")
        f = request("\nSelect output file:", RadioMenu(files))
        ds = NCDataset("results/outfiles/blooms/$(files[f])")
    end

    H = ds["H"][:][1]
    dz = ds["dz"][:][1]
    nIC, pIC, zIC, bIC, dIC, vIC, oIC0 = get_endpoints(["n", "p", "z", "b", "d", "v", "o"], ds)
    oIC = oIC0[:,:]
    nn, np, nz, nb, nd, nv = get_size([nIC, pIC, zIC, bIC, dIC, vIC])
    vmax_i = ds["vmax_i"][:]
    vmax_ij = ds["vmax_ij"][:,:]
    Kp_i = ds["Kp_i"][:]
    Kp_ij = ds["Kp_ij"][:,:]
    m_lp = ds["m_lp"][:]
    m_qp = ds["m_qp"][:]
    light = ds["light"][:]
    temp_fun = ds["temp_fun"][:]
    K_I = ds["K_I"][:][1]
    CMp = ds["CMp"][:,:]
    Fg_p = ds["Fg_p"][:]
    umax_i = ds["umax_i"][:]
    umax_ij = ds["umax_ij"][:,:]
    Km_i = ds["Km_i"][:]
    Km_ij = ds["Km_ij"][:,:]
    m_lb = ds["m_lb"][:]
    m_qb = ds["m_qb"][:]
    y_ij = ds["y_ij"][:,:]
    CM = ds["CM"][:,:]
    Fg_b = ds["Fg_b"][:]
    g_max = ds["g_max"][:]
    K_g = ds["K_g"][:]
    γ = ds["γ"][:]
    m_lz = ds["m_lz"][:]
    m_qz = ds["m_qz"][:]
    GrM = ds["GrM"][:,:]
    kappa_z = ds["kappa_z"][:]
    wd = ds["wd"][:,:]
    ngrid = ds["ngrid"][:][1]
    pulse = ds["pulse"][:][1]
    om_dist_mort = ds["om_dist_mort"][:]
    om_dist_lys = ds["om_dist_lys"][:]
    om_dist_vde = ds["om_dist_vde"][:]
    VM = ds["VM"][:,:]
    vly = ds["vly"][:][1]
    vbs = ds["vbs"][:][1]
    vde = ds["vde"][:][1]
    e_o = ds["e_o"][:][1]
    yo_ij = ds["yo_ij"][:,:]
    koverh = ds["koverh"][:][1]
    o2_sat = ds["o2_sat"][:][1]
    ml_boxes = ds["ml_boxes"][:][1]
    t_o2relax = ds["t_o2relax"][:][1]
    o2_deep = ds["o2_deep"][:][1]
    season = ds.attrib["Season"]
    lysis = ds.attrib["Lysis"]
    prev_fname = files[f]

    return H, dz, nIC, pIC, zIC, bIC, dIC, vIC, oIC, nn, np, nz, nb, nd, nv, 
           vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
           umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, CM, Fg_b,
           g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
           om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde, e_o, yo_ij,
           koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, lysis, prev_fname

end

function print_info(prms)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n days = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn, prms.nd, prms.days)
    println("File will be saved as: ", prms.fsaven)
    println("nt = ", prms.nt)

end


function set_logger(launch_time)

    loginfo = string(Dates.format(launch_time, "yyyymmdd_HHMM"), ".log")
    logger = activate_logger(loginfo)

    return logger

end


function activate_logger(loginfo)

    logger = TeeLogger(
        MinLevelLogger(FileLogger("logs/$loginfo"), Logging.Info),
        MinLevelLogger(FileLogger("logs/error.log"), Logging.Warn),
    ); 
    
    global_logger(logger)

    return logger 

end


function log_params(prms, season, lysis, graze)

    @info(
    """Model Params: 
    days:           $(prms.days)
    ts/day (dt):    $(prms.dt)
    total ts (nt):  $(prms.nt)                    
    recorded ts:    $(prms.nt) 
    nn:             $(prms.nn)                    
    np:             $(prms.np)                     
    nz:             $(prms.nz)        
    nb:             $(prms.nb)    
    nv:             $(prms.nv)      
    nd:             $(prms.nd)    
    nIC:            $(prms.nIC[1])                        
    pIC:            $(prms.pIC[1])
    zIC:            $(prms.zIC[1])
    bIC:            $(prms.bIC[1])
    dIC:            $(prms.dIC[1])
    vIC:            $(prms.vIC[1])
    oIC:            $(prms.oIC[1])
    vmax_i:         $(prms.vmax_i)
    vmax_ij:        $(prms.vmax_ij)
    umax_i:         $(prms.umax_i)
    umax_ij:        $(prms.umax_ij)
    Kp_i:           $(prms.Kp_i)
    Kp_ij:          $(prms.Kp_ij)
    Km_i:           $(prms.Km_i)
    Km_ij:          $(prms.Km_ij)
    y_ij:           $(prms.y_ij[1])
    K_g:            $(prms.K_g[1])
    g_max:          $(prms.g_max[1])
    γ:              $(prms.γ[1])
    m_lp:           $(prms.m_lp[1])
    m_qp:           $(prms.m_qp[1])
    m_lb:           $(prms.m_lb[1])
    m_qb:           $(prms.m_qb[1])
    m_lz:           $(prms.m_lz[1])
    m_qz:           $(prms.m_qz[1])
    wd:             $(prms.wd)
    Fg_p:           $(prms.Fg_p)
    Fg_b:           $(prms.Fg_b)
    light half sat: $(prms.K_I)
    pulse (1=no)    $(prms.pulse)  
    om_dist_mort:   $(prms.om_dist_mort)
    om_dist_lys:    $(prms.om_dist_lys)
    om_dist_vde:    $(prms.om_dist_vde)
    vly:            $(prms.vly)
    vde:            $(prms.vde)
    season (1=win)  $(season) 
    lysis (1=exp)   $(lysis)
    grazing (1=yes)$(graze)
    savefile:       $(prms.fsaven)
    """) 
    
    print_info(prms)


end

function get_matrix(Mtype, nd, nb, nn, np, nz) 
    
    if Mtype == "CM"
        M = build_consumption_matrix(nd, nb)
    elseif Mtype == "CMp"
        M = build_consumption_matrix(nn, np)
    else
        M = build_grazing_matrix(np, nb, nz)
    end

    return M
end


function check_for_negatives(RS)

    for i in eachindex(RS)
        for j in eachindex(RS[i])
            RS[i][j] = ifelse(RS[i][j] < 0 || RS[i][j] > 15, NaN, RS[i][j])
        end
    end

    return RS

end

function group_interactions(Cs, n)

    interactions = Any[]
    for (i, row) in enumerate(eachrow(Cs))
        for (j, col) in enumerate(eachcol(Cs))
            if Cs[i, j] > 0
                push!(interactions, [i, j])
            end
        end
    end

    return get_interaction_dict(interactions, n) 

end


function get_interaction_dict(interactions, n)

    out = Dict(name => Any[] for name in collect(1:1:n))
    for i in interactions
        for j in keys(out) 
            if i[1] == j 
                push!(out[j], i[2]) 
            end 
        end 
    end

    return out

end




function get_endpoints(vars, ds=nothing)

    endpoints = Vector{Any}()

    if ds !== nothing
        for v in vars
            if v != "o"
                append!(endpoints, [ds["$v"][:,:,end]])
            else
                append!(endpoints, [ds["$v"][:,end]])
            end
        end
    else
        for v in vars
            if v != "o"
                append!(endpoints, [v[:,:,end]])
            else
                append!(endpoints, [v[:,end]])
            end
        end
    end

    return endpoints
end

function get_final_year(ds, vars)
    # where year is 12 * 30 days
    final_yr = Vector{Any}()

    for v in vars
        if v != "o"
            append!(final_yr, [ds[v][:, :, end-7199:end]])
        else
            append!(final_yr, [ds[v][:, end-7199:end]])
        end
    end

    return final_yr

end

function get_final_three_cycles(ds, vars, pulse_freq)
    # where year is 12 * 30 days

    final_3 = Vector{Any}()
    ts = 20 * pulse_freq * 3

    for v in vars
        if v != "o"
            append!(final_3, [ds[v][:, :, end-(ts-1):end]])
        else
            append!(final_3, [ds[v][:, end-(ts-1):end]])
        end
    end

    return final_3

end

function get_zc(H)

    dz = 10
    zc = [dz/2 : dz : H - dz/2;] 

    return zc

end

function get_hmap_z_axis(depth, days, daily_data)

    z = Array{Float64}(undef, size(depth, 1), size(daily_data, 2))

    for i in eachindex(depth)
        for j in eachindex(days)
            z[i, j] = daily_data[i, j]
        end
    end

    return z'

end

function check_dir_exists(path)

    if isdir(path) == false
        mkdir(path)
    end

end

function set_zmax(state_var, num_state_var)

    zmax, z97 = 0, 0

    for s in range(1, num_state_var)
        state_var_max = maximum(state_var[:, s, :])
        state_var_max > zmax ? zmax = state_var_max : nothing

        s_arr = vec(state_var[:, s, :])
        state_var_97 = quantile(s_arr, 0.97)
        state_var_97 > z97 ? z97 = state_var_97 : nothing
    end

    return zmax, z97

end

function set_extinct_to_zero(ds)

    dss = copy(ds)
    ex = 10^-6
    dss .= ifelse.(dss .<= ex, 0.0, dss)

    return dss

end

function mean_over_time(state_vars, ds, season_num)

    if season_num == 1
        pulse_freq = 10
        means_over_time = get_cycle_means(state_vars, pulse_freq, ds)
    else
        pulse_freq = 30
        means_over_time = get_cycle_means(state_vars, pulse_freq, ds)
    end

    return means_over_time
    
end

function get_bloom_means(vars, ds)

    bloom_start_ts = 10952

    bloom_means = Vector{Any}()
    for v in vars
        if v != "o"
            period = ds["$v"][:,:,bloom_start_ts:end]
            period_mean = dropdims(mean(period, dims=3), dims=3)
            append!(bloom_means, [period_mean])
        else
            append!(bloom_means, [mean(ds["$v"][:,bloom_start_ts:end], dims=2)])
        end
    end

    return bloom_means

end

function get_cycle_means(vars, pulse_freq, ds)

    #mean of last 10 pulse cycles given 20 recorded ts per day
    cycles = 10
    days = pulse_freq * cycles
    ts = days * 20

    cycle_mean = Vector{Any}()

    for v in vars
        if v != "o"
            period = ds["$v"][:,:,end-ts:end]
            period_mean = dropdims(mean(period, dims=3), dims=3)
            append!(cycle_mean, [period_mean])
        else
            append!(cycle_mean, [mean(ds["$v"][:,end-ts:end], dims=2)])
        end
    end


    return cycle_mean

end

function load_bloom(vars, ds)

    bloom_data = Vector{Any}()

    for v in vars
        if v != "o"
            append!(bloom_data, [ds[v][:, :, :]])
        else
            append!(bloom_data, [ds[v][:, :]])
        end
    end

    return bloom_data

end


function get_size(arr)

    out = Vector{Int}()
    
    for a in arr
        append!(out, size(a, 2))
    end

    return out

end


function get_nonzero_axes(Mat)

    Cs = sparse(Mat)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


#------------------------------------------------------------------------
#                             PLOT UTILS 
#------------------------------------------------------------------------

function get_plot_vars()

    bcols = ["cyan3", "darkorange", "indigo", "coral4", "lightcyan4", "magenta2", "thistle", "seagreen4",
            "darkkhaki", "purple", "black",  "yellow3", "navajowhite4",  "coral4", "orange2", "orangered4", "yellow3", 
            "lightyellow4", "goldenrod4", "slateblue4", "mediumpurple3"]
    dcols = ["blue3", "black", "maroon", "coral", "orange3", "silver", "magenta3", "cyan3", "seagreen", "mediumpurple3"]
    pcols = ["hotpink2", "darkgreen","red4", "cyan4", "gold3", "black", "brown", "wheat2", "mediumpurple3", "darkseagreen" ]
    ncols = ["blue2"]
    zcols = ["black", "slategray4", "deeppink3", "sienna", "mediumpurple3", "darkseagreen", "snow4", 
            "silver", "salmon", "coral4", "orange2", "orangered4", "yellow3", "lightyellow4", "goldenrod4",
            "chartreuse", "lightseagreen", "blueviolet", "slateblue4", "magenta4" ]
    ab = 0.8
    ab_ext = 0.8
    ls = 4
    lfs = 9
    lg = :bottomright
    
    return bcols, dcols, pcols, ncols, zcols, ab, ab_ext, ls, lfs, lg

end




