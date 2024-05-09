 module NPZBDV_1D
    
    import REPL
    using REPL.TerminalMenus
    using Logging, LoggingExtras 
    using Dates, Printf, Parameters
    using SparseArrays, Distributions, LinearAlgebra
    using Statistics, StatsBase, Random, StableRNGs
    using DataFrames, NCDatasets, CSV
    using Plots, ColorSchemes, Colors, LaTeXStrings, Measures
    using Distributed

    include("model.jl")
    include("params.jl")
    include("physics.jl")
    include("consumption_matrix.jl")
    include("grazing_matrix.jl")
    include("traits.jl")
    include("integrate.jl")
    include("rstar.jl")
    include("nutrient_pulse.jl")
    # include("plotting/heatmaps.jl")
    include("plotting/state_var_plots.jl")
    include("plotting/timeseries_plots.jl")
    include("plotting/rstar_plots.jl")
    include("utils/utils.jl")
    include("utils/save_utils.jl")

    GlobalRNG = StableRNG(123)
    addprocs(15, exeflags = "--project=$(Base.active_project())")
    println("\n > Number of cores: ", nprocs())
    println(" > Number of workers: ", nworkers())


    #------------------------------------------------------------------------------------------------------------#
    #   TIMES AND LOGS
    #------------------------------------------------------------------------------------------------------------#
    println(message("START"))
    run_type = request(message("RT1"), RadioMenu(message("RT2")))

    if run_type == 3
        years = 1
        dt = 0.001
    else 
        println(message("TM"))
        user_input = readline()
        years = parse(Int64, user_input)
        dt = 0.01
    end
            
    days = years * 366
    nrec = days * 20
    nt = Int(days/dt)
    logger = set_logger(now())

    ##--------------------------------------------------------------------------------------------------#
    ##   LOAD MODEL
    ##--------------------------------------------------------------------------------------------------#
    if run_type == 1 
        # prompts user to select model from src/prescribed
        bloom=false
        lysis = request(message("LY1"), RadioMenu(message("LY2")))
        graze = request(message("GZ1"), RadioMenu(message("GZ2")))
        files = readdir("src/prescribed")
        f = request("\nSelect model to run:", RadioMenu(files))
        include("prescribed/$(files[f])")

        N, P, Z, B, D, V, O, track_time = run_NPZBDV(params, season)
        log_params(params, season, lysis, graze)
        plot_state_vars(fsaven, season, lysis, graze)
        exit()

    elseif run_type == 2
        # Continues from end of previous runs (prompts user to select .nc file in results/outfiles )
        continue_type = request(message("CT1"), RadioMenu(message("CT2")))

        H, dz, nIC, pIC, zIC, bIC, dIC, vIC, oIC, nn, np, nz, nb, nd, nv,
        vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
        umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, CM, Fg_b,
        g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
        om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde, e_o, yo_ij,
        koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, lysis, prev_fname = get_previous_params(continue_type)

        fsaven = continuation_savefile(prev_fname)
        params = Prms(
                days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, nv, pIC, bIC, zIC, nIC, dIC, vIC, oIC, 
                vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
                umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
                g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
                om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde,
                e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven, lysis, graze
        )

        N, P, Z, B, D, V, O, track_time = run_NPZBDV(params, season)
        log_params(params, season, lysi, graze)
        plot_state_vars(fsaven, season, lysis, graze)
        exit()

    elseif run_type == 3
    # Takes output of previous run and tracks single pulse over 365 days
        H, dz, nIC, pIC, zIC, bIC, dIC, vIC, oIC, nn, np, nz, nb, nd, nv,
        vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
        umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, CM, Fg_b,
        g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
        om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde, e_o, yo_ij,
        koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, prev_fname = get_previous_params(continue_type)

        years = 1
        bloom=true

        fsaven = continuation_savefile(prev_fname, bloom)
        params = Prms(
                days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, nv, pIC, bIC, zIC, nIC, dIC, vIC, oIC, 
                vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
                umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
                g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
                om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde,
                e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven, lysis, graze
        )

        N, P, Z, B, D, V, O, track_time = run_NPZBD(params, season, bloom)  
        log_params(params, season, lysis, graze)  
        plot_state_vars(fsaven, season, lysis, graze, bloom) 
        exit()

    end

end

