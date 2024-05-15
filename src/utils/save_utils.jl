
using DataFrames, NCDatasets


function set_savefiles(launch_time, season, pulse, years, np, nz, nb, nd, nc, nv)

    season == 1 ? season_str = "Wi" : season_str = "Su"

    if pulse == 1
        ptype = ""
    elseif pulse == 2
        ptype = "PP"
    else
        ptype = "SP"
    end

    fsave = "results/outfiles/"
    fsaven = string(fsave, Dates.format(launch_time, "yymmdd_HH:MM"), "_$(season_str)$(years)y$(ptype)_$(nn)N$(nc)C$(np)P$(nz)Z$(nb)B$(nd)D$(nv)V.nc")

    return fsaven

end

function continuation_savefile(prev_fname, bloom=false)

    fsave = "results/outfiles/"
    if bloom == true
        fname = replace(prev_fname, "_ep.nc" => ".nc")
        fsaven = string(fsave, "blooms/blm_", fname)
    else
        fname = replace(prev_fname, "_ep.nc" => ".nc")
        fsaven = string(fsave, "cont_", prev_fname)
    end

    return fsaven

end


function check_subfolder_exists(filename, parent_folder)

    rgx = r"_[^_]+$"
    sub_folder = match(rgx, filename)

    isdir(parent_folder * sub_folder.match) ? 
        folder = parent_folder * sub_folder.match : folder = mkdir(parent_folder * sub_folder.match)

    return folder

end


function save_full_run(p, b, z, n, c, d, v, o, timet, tst, tfn, prms, season_num, lysis)

    outdir = "/home/lee/Dropbox/Development/NPZBDV_1D/"
    season_num == 1 ? season = "winter" : season = "summer"

    if prms.pulse == 1
        pulse_type = "None (steady state)"
    elseif prms.pulse == 2
        pulse_type = "Peroidic nutrient pulse"
    else 
        pulse_type = "Semi-stochastic nutrient pulse"
    end

    path = joinpath(outdir, prms.fsaven)
    println("\nSaving output file to: ", path)
    f = NCDataset(path, "c") 

    # define the dim of p, b, z, n, c, d
    defDim(f,"np",prms.np)
    defDim(f,"nb",prms.nb)
    defDim(f,"nz",prms.nz)
    defDim(f,"nn",prms.nn)
    defDim(f,"nc",prms.nc)
    defDim(f,"nd",prms.nd)
    defDim(f,"nv",prms.nv)


    # define the dim of the depth
    defDim(f,"ndepth",prms.ngrid)
    defDim(f,"ndepth1",prms.ngrid+1)

    # define the dim of the time length 
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    defDim(f,"nrec",nrec1)

    nprey = prms.np + prms.nb
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBD 1D model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 
    f.attrib["Season"] = season
    f.attrib["Pulse Type"] = pulse_type
    f.attrib["Lysis"] = lysis

    # simulated results
    w = defVar(f,"p",Float64,("ndepth" ,"np","nrec"))
    w[:,:,:] = p
    w.attrib["units"] = "mmol N/m3 biomass"

    w = defVar(f,"b",Float64,("ndepth" ,"nb","nrec"))
    w[:,:,:] = b
    w.attrib["units"] = "mmol N/m3 biomass"

    w = defVar(f,"z",Float64,("ndepth" ,"nz","nrec"))
    w[:,:,:] = z
    w.attrib["units"] = "mmol N/m3 biomass"
    
    w = defVar(f,"n",Float64,("ndepth" ,"nn","nrec"))
    w[:,:,:] = n
    w.attrib["units"] = "mmol N/m3"

    w = defVar(f,"c",Float64,("ndepth" ,"nc","nrec"))
    w[:,:,:] = c
    w.attrib["units"] = "mmol C/m3"

    w = defVar(f,"d",Float64,("ndepth" ,"nd","nrec"))
    w[:,:,:] = d
    w.attrib["units"] = "mmol N/m3"

    w = defVar(f,"v",Float64,("ndepth" ,"nv","nrec"))
    w[:,:,:] = v
    w.attrib["units"] = "mmol N/m3"
    
    w = defVar(f,"o",Float64,("ndepth" ,"nrec"))
    w[:,:] = o
    w.attrib["units"] = "mmol/m3 O2"

    # --------------------------------------------------
    
    w = defVar(f,"pIC",Float64,("ndepth","np"))
    w[:,:] = prms.pIC
    w.attrib["units"] = "mmol N/m3 biomass"

    w = defVar(f,"bIC",Float64,("ndepth","nb"))
    w[:,:] = prms.bIC
    w.attrib["units"] = "mmol N/m3 biomass"

    w = defVar(f,"zIC",Float64,("ndepth","nz"))
    w[:,:] = prms.zIC
    w.attrib["units"] = "mmol N/m3 biomass"

    w = defVar(f,"nIC",Float64,("ndepth","nn"))
    w[:,:] = prms.nIC
    w.attrib["units"] = "mmol N/m3"

    w = defVar(f,"cIC",Float64,("ndepth","nc"))
    w[:,:] = prms.cIC
    w.attrib["units"] = "mmol C/m3"

    w = defVar(f,"dIC",Float64,("ndepth","nd"))
    w[:,:] = prms.dIC
    w.attrib["units"] = "mmol N/m3"

    w = defVar(f,"vIC",Float64,("ndepth","nv"))
    w[:,:] = prms.vIC
    w.attrib["units"] = "mmol N/m3"

    # --------------------------------------------------
    
    w = defVar(f, "timet", Float64, ("nrec",))
    w[:] = timet
    w.attrib["units"] = "days"

    w = defVar(f, "H", Int, ())
    w[:] = prms.H
    w.attrib["units"] = "m; total height"

    w = defVar(f, "dz", Int, ())
    w[:] = prms.dz
    w.attrib["units"] = "m; box height"

    w = defVar(f, "np", Int, ())
    w[:] = prms.np
    w.attrib["units"] = "num phytoplankton"

    w = defVar(f, "nb", Int, ())
    w[:] = prms.nb
    w.attrib["units"] = "num heterotrophic bacteria"

    w = defVar(f, "nz", Int, ())
    w[:] = prms.nz
    w.attrib["units"] = "num zooplankton"

    w = defVar(f, "nn", Int, ())
    w[:] = prms.nn
    w.attrib["units"] = "num inorganic nutrient pools"

    w = defVar(f, "nd", Int, ())
    w[:] = prms.nd
    w.attrib["units"] = "num organic nutrient pools"

    w = defVar(f, "nv", Int, ())
    w[:] = prms.nv
    w.attrib["units"] = "num virus pools"

    w = defVar(f, "vly", Float64, ())
    w[:] = prms.vly
    w.attrib["units"] = "viral lysis rate"

    w = defVar(f, "vbs", Int, ())
    w[:] = prms.vbs
    w.attrib["units"] = "viral burst size"

    w = defVar(f, "vde", Float64, ())
    w[:] = prms.vde
    w.attrib["units"] = "viral decay rate"

    w = defVar(f, "pulse", Int, ())
    w[:] = prms.pulse
    w.attrib["nutrient_pulse"] = "1: no pulse, 2: Total N,D redistributed along water col"

    w = defVar(f, "kappa_z", Float64, ("ndepth1",))
    w[:] = prms.kappa_z
    w.attrib["units"] = "vertical water velocity"
    
    w = defVar(f, "wd", Float64, ("ndepth1","nd"))
    w[:] = prms.wd
    w.attrib["units"] = "sinking rate"

    w = defVar(f, "temp_fun", Float64, ("ndepth",))
    w[:] = prms.temp_fun
    w.attrib["units"] = "temp mod to metabolic rate"

    w = defVar(f, "light", Float64, ("ndepth",))
    w[:] = prms.light
    w.attrib["units"] = "irradiance profile applied to water col, with eutrophic zone factored in"

    w = defVar(f, "e_o", Float64, ())
    w[:] = prms.e_o
    w.attrib["units"] = "production of O2 (excretion) - mol O2/mol N uptake"

    w = defVar(f, "K_I", Float64, ())
    w[:] = prms.K_I
    w.attrib["units"] = "light half-sat constant"

    w = defVar(f, "ngrid", Int, ())
    w[:] = prms.ngrid
    w.attrib["units"] = "num of horizontal grid cells"

    w = defVar(f, "yo_ij", Float64, ("nd","nb"))
    w[:,:] = prms.yo_ij
    w.attrib["units"] = "oxygen yield - mol B/mol O2"

    w = defVar(f, "koverh", Float64, ())
    w[:] = prms.koverh
    w.attrib["units"] = "gas transfer coefficient for each box comprising the mixed layer"

    w = defVar(f, "o2_sat", Float64, ())
    w[:] = prms.o2_sat
    w.attrib["units"] = "O2 half-sat constant (mmol/m3)"

    w = defVar(f, "t_o2relax", Float64, ())
    w[:] = prms.t_o2relax
    w.attrib["units"] = "deep oxygen relaxation time (1/day)"

    w = defVar(f, "o2_deep", Float64, ())
    w[:] = prms.o2_deep
    w.attrib["units"] = "deep oxygen relaxation mmol/m3"

    w = defVar(f, "ml_boxes", Int, ())
    w[:] = prms.ml_boxes
    w.attrib["units"] = "discrete n of boxes in the mixed layer, close to 100m total sum"

    # --------------------------------------------------

    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = prms.CM
    w.attrib["units"] = "Consumption Matrix (nb x nd)"

    w = defVar(f,"CMp",Float64,("nn","np"))
    w[:,:] = prms.CMp
    w.attrib["units"] = "Consumption Matrix (np x nn)"

    w = defVar(f,"GrM",Float64,("nz","nprey"))
    w[:,:] = prms.GrM
    w.attrib["units"] = "Grazing Matrix"

    w = defVar(f,"VM",Float64,("nv","nb"))
    w[:,:] = prms.VM
    w.attrib["units"] = "Lysing Matrix"

    w = defVar(f,"om_dist_mort",Float64,("nd",))
    w[:] = prms.om_dist_mort 
    w.attrib["units"] = "Distribution of OM from mortality"

    w = defVar(f,"om_dist_lys",Float64,("nd",))
    w[:] = prms.om_dist_lys
    w.attrib["units"] = "Distribution of OM from viral lysis"

    w = defVar(f,"om_dist_vde",Float64,("nd",))
    w[:] = prms.om_dist_vde
    w.attrib["units"] = "Distribution of OM from viral decay"
    

    # --------------------------------------------------
    
    w = defVar(f, "vmax_i", Float64, ("np",))
    w[:] = prms.vmax_i
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "vmax_ij", Float64, ("nn", "np"))
    w[:,:] = prms.vmax_ij
    w.attrib["units"] = "per n; max uptake rate"

    w = defVar(f, "Kp_i", Float64, ("np",))
    w[:] = prms.Kp_i
    w.attrib["units"] = "mmol/m3; intrinsic half-sat of P_i before trade-off applied"

    w = defVar(f, "Kp_ij", Float64, ("nn", "np"))
    w[:,:] = prms.Kp_ij
    w.attrib["units"] = "mmol/m3; half-sat of P_i on N_j"

    w = defVar(f, "m_lp", Float64, ("np",))
    w[:] = prms.m_lp
    w.attrib["units"] = "m3/mmol; linear death rate of p"

    w = defVar(f, "m_qp", Float64, ("np",))
    w[:] = prms.m_qp
    w.attrib["units"] = "m3/mmol; quadratic death rate of p"

    w = defVar(f, "Fg_p", Float64, ("np",))
    w[:] = prms.Fg_p
    w.attrib["units"] = "per p; fraction proteome assigned to growth"

    # --------------------------------------------------

    w = defVar(f, "umax_i", Float64, ("nd",))
    w[:] = prms.umax_i
    w.attrib["units"] = "per d; max uptake rate of d_i before trade-off applied"
    
    w = defVar(f, "umax_ij", Float64, ("nd", "nb"))
    w[:,:] = prms.umax_ij
    w.attrib["units"] = "max uptake rate of b_i on d_j"

    w = defVar(f, "Km_i", Float64, ("nd",))
    w[:] = prms.Km_i
    w.attrib["units"] = "mmol/m3; half-sat of d_i before trade-off applied"
    
    w = defVar(f, "Km_ij", Float64, ("nd", "nb"))
    w[:] = prms.Km_ij
    w.attrib["units"] = "mmol/m3; half-sat of b_i on d_j"
    
    w = defVar(f, "m_lb", Float64, ("nb",))
    w[:] = prms.m_lb
    w.attrib["units"] = "m3/mmol; linear death rate of b"

    w = defVar(f, "m_qb", Float64, ("nb",))
    w[:] = prms.m_qb
    w.attrib["units"] = "m3/mmol; quadratic death rate of b"

    w = defVar(f, "y_ij", Float64, ("nd","nb"))
    w[:,:] = prms.y_ij
    w.attrib["units"] = "per d; max yield rate"

    w = defVar(f, "Fg_b", Float64, ("nb",))
    w[:] = prms.Fg_b
    w.attrib["units"] = "per b; fraction proteome assigned to growth"

    # --------------------------------------------------

    w = defVar(f, "K_g", Float64, ("nz",))
    w[:] = prms.K_g
    w.attrib["units"] = "m3/mmol; half-sat rate of z"

    w = defVar(f, "g_max", Float64, ("nz",))
    w[:] = prms.g_max
    w.attrib["units"] = "m3/mmol; max uptake rate of z"
    
    w = defVar(f, "γ", Float64, ("nz",))
    w[:] = prms.γ
    w.attrib["units"] = "fraction of digestion of z"
    
    w = defVar(f, "m_lz", Float64, ("nz",))
    w[:] = prms.m_lz
    w.attrib["units"] = "m3/mmol; linear death rate of z"

    w = defVar(f, "m_qz", Float64, ("nz",))
    w[:] = prms.m_qz
    w.attrib["units"] = "m3/mmol; quadratic death rate of z"

    # --------------------------------------------------
    
    close(f)

end


function save_endpoints(p, b, z, n, c, d, v, o, timet, tst, tfn, prms, season, lysis)

    outdir = "/home/lee/Dropbox/Development/NPZBDV_1D/"
    ep_path = replace(prms.fsaven, "results/outfiles" => "results/outfiles/endpoints", ".nc" => "_ep.nc")
    path = joinpath(outdir, ep_path) 
    println("\nSaving endpoints to: ", path)

    if prms.pulse == 1
        pulse_type = "None (steady state)"
    elseif prms.pulse == 2
        pulse_type = "Peroidic nutrient pulse"
    else 
        pulse_type = "Semi-stochastic nutrient pulse"
    end

    n, c, p, z, b, d, v, o = get_endpoints([n, c, p, z, b, d, v, o])
    # n, p, z, b, d, o = ep[1], ep[2], ep[3], ep[4], ep[5], ep[6]

    f = NCDataset(path, "c") 

    # define the dim of p, b, z, n, d
    defDim(f, "np", prms.np)
    defDim(f, "nb", prms.nb)
    defDim(f, "nz", prms.nz)
    defDim(f, "nn", prms.nn)
    defDim(f, "nc", prms.nc)
    defDim(f, "nd", prms.nd)
    defDim(f, "no", 1)
    defDim(f,"nv", prms.nv)

    # define the dim of the depth
    defDim(f,"ndepth",89)
    defDim(f,"ndepth1",90)

    # define the dim of the time length
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    nprey = prms.np + prms.nb
    
    defDim(f,"nrec",nrec1)
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBD 1D model endpoints"
    f.attrib["Season"] = season
    f.attrib["Pulse Type"] = pulse_type
    f.attrib["Lysis"] = lysis

    # simulated results
    w = defVar(f,"p",Float64,("ndepth" ,"np"))
    w[:,:] = p
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"b",Float64,("ndepth" ,"nb"))
    w[:,:] = b
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"z",Float64,("ndepth" ,"nz"))
    w[:,:] = z
    w.attrib["units"] = "mmol/m3 N biomass"
    
    w = defVar(f,"n",Float64,("ndepth" ,"nn"))
    w[:,:] = n
    w.attrib["units"] = "mmol/m3 N"

    w = defVar(f,"c",Float64,("ndepth" ,"nc"))
    w[:,:] = c
    w.attrib["units"] = "mmol/m3 C"

    w = defVar(f,"d",Float64,("ndepth" ,"nd"))
    w[:,:] = d
    w.attrib["units"] = "mmol/m3 N"

    w = defVar(f,"v",Float64,("ndepth" ,"nv"))
    w[:,:] = v
    w.attrib["units"] = "mmol/m3 N"

    w = defVar(f, "vly", Float64, ())
    w[:] = prms.vly
    w.attrib["units"] = "viral lysis rate"

    w = defVar(f, "vbs", Int, ())
    w[:] = prms.vbs
    w.attrib["units"] = "viral burst size"

    w = defVar(f, "vde", Float64, ())
    w[:] = prms.vde
    w.attrib["units"] = "viral decay rate"
    
    w = defVar(f,"o",Float64,("ndepth" ,"no"))
    w[:,:] = o
    w.attrib["units"] = "mmol/m3 O2"

  # --------------------------------------------------
    
    w = defVar(f,"pIC",Float64,("ndepth","np"))
    w[:,:] = prms.pIC
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"bIC",Float64,("ndepth","nb"))
    w[:,:] = prms.bIC
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"zIC",Float64,("ndepth","nz"))
    w[:,:] = prms.zIC
    w.attrib["units"] = "mmol/m3 N biomass"

    w = defVar(f,"nIC",Float64,("ndepth","nn"))
    w[:,:] = prms.nIC
    w.attrib["units"] = "mmol/m3 N"

    w = defVar(f,"cIC",Float64,("ndepth","nc"))
    w[:,:] = prms.cIC
    w.attrib["units"] = "mmol/m3 C"

    w = defVar(f,"dIC",Float64,("ndepth","nd"))
    w[:,:] = prms.dIC
    w.attrib["units"] = "mmol/m3 N"

    w = defVar(f,"vIC",Float64,("ndepth","nv"))
    w[:,:] = prms.vIC
    w.attrib["units"] = "mmol/m3 N"

    # --------------------------------------------------
    
    w = defVar(f, "timet", Float64, ("nrec",))
    w[:] = timet
    w.attrib["units"] = "days"

    w = defVar(f, "H", Int, ())
    w[:] = prms.H
    w.attrib["units"] = "m; total height"

    w = defVar(f, "dz", Int, ())
    w[:] = prms.dz
    w.attrib["units"] = "m; box height"

    w = defVar(f, "np", Int, ())
    w[:] = prms.np
    w.attrib["units"] = "num phytoplankton"

    w = defVar(f, "nb", Int, ())
    w[:] = prms.nb
    w.attrib["units"] = "num heterotrophic bacteria"

    w = defVar(f, "nz", Int, ())
    w[:] = prms.nz
    w.attrib["units"] = "num zooplankton"

    w = defVar(f, "nn", Int, ())
    w[:] = prms.nn
    w.attrib["units"] = "num inorganic nitrogen pools"

    w = defVar(f, "nc", Int, ())
    w[:] = prms.nc
    w.attrib["units"] = "num inorganic carbon pools"

    w = defVar(f, "nd", Int, ())
    w[:] = prms.nd
    w.attrib["units"] = "num organic DOM pools"

    w = defVar(f, "nv", Int, ())
    w[:] = prms.nv
    w.attrib["units"] = "num virus pools"

    w = defVar(f, "pulse", Int, ())
    w[:] = prms.pulse
    w.attrib["nutrient_pulse"] = "1: no pulse, 2: Total N,D redistributed along water col"

    w = defVar(f, "kappa_z", Float64, ("ndepth1",))
    w[:] = prms.kappa_z
    w.attrib["units"] = "vertical water velocity"
    
    w = defVar(f, "wd", Float64, ("ndepth1","nd"))
    w[:] = prms.wd
    w.attrib["units"] = "sinking rate"

    w = defVar(f, "temp_fun", Float64, ("ndepth",))
    w[:] = prms.temp_fun
    w.attrib["units"] = "temp mod to metabolic rate"

    w = defVar(f, "light", Float64, ("ndepth",))
    w[:] = prms.light
    w.attrib["units"] = "irradiance profile applied to water col, with eutrophic zone factored in"

    w = defVar(f, "e_o", Float64, ())
    w[:] = prms.e_o
    w.attrib["units"] = "production of O2 (excretion) - mol O2/mol N uptake"

    w = defVar(f, "K_I", Float64, ())
    w[:] = prms.K_I
    w.attrib["units"] = "light half-sat constant"

    w = defVar(f, "ngrid", Int, ())
    w[:] = prms.ngrid
    w.attrib["units"] = "num of horizontal grid cells"

    w = defVar(f, "yo_ij", Float64, ("nd","nb"))
    w[:,:] = prms.yo_ij
    w.attrib["units"] = "oxygen yield - mol B/mol O2"

    w = defVar(f, "koverh", Float64, ())
    w[:] = prms.koverh
    w.attrib["units"] = "gas transfer coefficient for each box comprising the mixed layer"

    w = defVar(f, "o2_sat", Float64, ())
    w[:] = prms.o2_sat
    w.attrib["units"] = "O2 half-sat constant (mmol/m3)"

    w = defVar(f, "t_o2relax", Float64, ())
    w[:] = prms.t_o2relax
    w.attrib["units"] = "deep oxygen relaxation time (1/day)"

    w = defVar(f, "o2_deep", Float64, ())
    w[:] = prms.o2_deep
    w.attrib["units"] = "deep oxygen relaxation mmol/m3"

    w = defVar(f, "ml_boxes", Int, ())
    w[:] = prms.ml_boxes
    w.attrib["units"] = "discrete n of boxes in the mixed layer, close to 100m total sum"

    # --------------------------------------------------

    w = defVar(f,"CM",Float64,("nd","nb"))
    w[:,:] = prms.CM
    w.attrib["units"] = "Consumption Matrix (nb x nd)"

    w = defVar(f,"CMp",Float64,("nn","np"))
    w[:,:] = prms.CMp
    w.attrib["units"] = "Consumption Matrix (np x nn)"

    w = defVar(f,"GrM",Float64,("nz","nprey"))
    w[:,:] = prms.GrM
    w.attrib["units"] = "Grazing Matrix"

    w = defVar(f,"VM",Float64,("nv","nb"))
    w[:,:] = prms.VM
    w.attrib["units"] = "Lysing Matrix"

    w = defVar(f,"om_dist_mort",Float64,("nd",))
    w[:] = prms.om_dist_mort 
    w.attrib["units"] = "Distribution of OM from mortality"

    w = defVar(f,"om_dist_lys",Float64,("nd",))
    w[:] = prms.om_dist_lys
    w.attrib["units"] = "Distribution of OM from viral lysis"

    w = defVar(f,"om_dist_vde",Float64,("nd",))
    w[:] = prms.om_dist_vde
    w.attrib["units"] = "Distribution of OM from viral decay"

    # --------------------------------------------------
    
    w = defVar(f, "vmax_i", Float64, ("np",))
    w[:] = prms.vmax_i
    w.attrib["units"] = "m3/mmol/d; max growth rate of p"

    w = defVar(f, "vmax_ij", Float64, ("nn", "np"))
    w[:,:] = prms.vmax_ij
    w.attrib["units"] = "per n; max uptake rate"

    w = defVar(f, "Kp_i", Float64, ("np",))
    w[:] = prms.Kp_i
    w.attrib["units"] = "mmol/m3; intrinsic half-sat of P_i before trade-off applied"

    w = defVar(f, "Kp_ij", Float64, ("nn", "np"))
    w[:,:] = prms.Kp_ij
    w.attrib["units"] = "mmol/m3; half-sat of P_i on N_j"

    w = defVar(f, "m_lp", Float64, ("np",))
    w[:] = prms.m_lp
    w.attrib["units"] = "m3/mmol; linear death rate of p"

    w = defVar(f, "m_qp", Float64, ("np",))
    w[:] = prms.m_qp
    w.attrib["units"] = "m3/mmol; quadratic death rate of p"

    w = defVar(f, "Fg_p", Float64, ("np",))
    w[:] = prms.Fg_p
    w.attrib["units"] = "per p; fraction proteome assigned to growth"

    # --------------------------------------------------

    w = defVar(f, "umax_i", Float64, ("nd",))
    w[:] = prms.umax_i
    w.attrib["units"] = "per d; max uptake rate of d_i before trade-off applied"
    
    w = defVar(f, "umax_ij", Float64, ("nd", "nb"))
    w[:,:] = prms.umax_ij
    w.attrib["units"] = "max uptake rate of b_i on d_j"

    w = defVar(f, "Km_i", Float64, ("nd",))
    w[:] = prms.Km_i
    w.attrib["units"] = "mmol/m3; half-sat of d_i before trade-off applied"
    
    w = defVar(f, "Km_ij", Float64, ("nd", "nb"))
    w[:] = prms.Km_ij
    w.attrib["units"] = "mmol/m3; half-sat of b_i on d_j"
    
    w = defVar(f, "m_lb", Float64, ("nb",))
    w[:] = prms.m_lb
    w.attrib["units"] = "m3/mmol; linear death rate of b"

    w = defVar(f, "m_qb", Float64, ("nb",))
    w[:] = prms.m_qb
    w.attrib["units"] = "m3/mmol; quadratic death rate of b"

    w = defVar(f, "y_ij", Float64, ("nd","nb"))
    w[:,:] = prms.y_ij
    w.attrib["units"] = "per d; max yield rate"

    w = defVar(f, "Fg_b", Float64, ("nb",))
    w[:] = prms.Fg_b
    w.attrib["units"] = "per b; fraction proteome assigned to growth"

    # --------------------------------------------------

    w = defVar(f, "K_g", Float64, ("nz",))
    w[:] = prms.K_g
    w.attrib["units"] = "m3/mmol; half-sat rate of z"

    w = defVar(f, "g_max", Float64, ("nz",))
    w[:] = prms.g_max
    w.attrib["units"] = "m3/mmol; max uptake rate of z"
    
    w = defVar(f, "γ", Float64, ("nz",))
    w[:] = prms.γ
    w.attrib["units"] = "fraction of digestion of z"
    
    w = defVar(f, "m_lz", Float64, ("nz",))
    w[:] = prms.m_lz
    w.attrib["units"] = "m3/mmol; linear death rate of z"

    w = defVar(f, "m_qz", Float64, ("nz",))
    w[:] = prms.m_qz
    w.attrib["units"] = "m3/mmol; quadratic death rate of z"
    
    close(f)


end


# ds = NCDataset("/home/lee/Dropbox/Development/NPZBD_1D/results/outfiles/Wi100y_8P6Z13B5D_230915_00:03.nc")
# save_endpoints(ds,"summer", 100, 8, 13, 6, 1, 5, 1, "Wi100y_8P6Z13B5D_230915_00:03_ep.nc")