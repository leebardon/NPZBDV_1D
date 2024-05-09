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

        run_type = request(message("ST2"), RadioMenu(message("ST1")))

        simulation_time = request(message("TM2"), RadioMenu(message("TM1")))
        if simulation_time == 1
            years = 0
            days = 366
            nrec = 7320
        elseif simulation_time == 2
            years = 2
            days = 732
            nrec = 14640
        elseif simulation_time == 3
            years = 10
            days = 3660
            nrec = 73200
        elseif simulation_time == 4
            years = 30
            days = 10980
            nrec = 219600
        elseif simulation_time == 5
            years = 50
            days = 18300
            nrec = 366000
        else 
            years = 100
            days = 36600
            nrec = 732000
        end

        # nrec = 20 per day when dt = 0.01 (i.e. 100 ts per day, tracking recorded every 5 ts)
        simulation_time == 1 ? dt = 0.001 : dt = 0.01
        nt = Int(days/dt)
        logger = set_logger(now())

        if run_type == 1 
            ##--------------------------------------------------------------------------------------------------#
            ##   LOAD MODEL
            ##--------------------------------------------------------------------------------------------------#
            bloom=false
            lysis = request(message("LY2"), RadioMenu(message("LY1")))
            graze = request(message("GZ2"), RadioMenu(message("GZ1")))
            files = readdir("src/prescribed")
            f = request("\nSelect model to run:", RadioMenu(files))
            include("prescribed/$(files[f])")

            ##--------------------------------------------------------------------------------------------------#
            ##   RUN MODEL
            ##--------------------------------------------------------------------------------------------------#
            N, P, Z, B, D, V, O, track_time = run_NPZBDV(params, season)
            log_params(params, season, lysis, graze)
            plot_state_vars(fsaven, season, lysis, graze)
            exit()

        elseif run_type == 2
            # Continues from end of previous runs (prompts user to select .nc file in results/outfiles )
            continue_type = request(message("CT2"), RadioMenu(message("CT1")))

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

        elseif run_type == 4
            # Takes output of previous run and tracks single pulse over 365 days
            H, dz, nIC, pIC, zIC, bIC, dIC, vIC, oIC, nn, np, nz, nb, nd, nv,
            vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p, 
            umax_i, umax_ij, Km_i, Km_ij, m_lb, m_qb, y_ij, CM, Fg_b,
            g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
            om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde, e_o, yo_ij,
            koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, season, prev_fname = get_previous_params(continue_type)

            years = 0
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

