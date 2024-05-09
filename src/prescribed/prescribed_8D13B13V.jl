

# include("../model.jl")
# include("../params.jl")
# include("../physics.jl")
# include("../consumption_matrix.jl")
# include("../grazing_matrix.jl")
# include("../traits.jl")
# include("../integrate.jl")
# include("../rstar.jl")
# include("../nutrient_pulse.jl")

nd = 8
nb = 13
np = 6
nz = 3
nn = 1
nv = 13

println(message("ENV"))
pulse = request(message("P2"), RadioMenu(message("P1")))
season = request(message("SE2"), RadioMenu(message("SE1")))

fsaven = set_savefiles(now(), season, pulse, years, np, nz, nb, nd, nv)

#------------------------------------------------------------------------------------------------------------#
#                                           GRID SETUP
#------------------------------------------------------------------------------------------------------------#
H = 890                         # depth at SPOT (m)
dz = 10                         # height per box
ngrid = Int(H/dz)               # number of boxes
zc = [dz/2 : dz : H - dz/2;]    # centered depth 
zf = [0 : dz : H;]              # face depth; 1 longer than zc


# -----------------------------------------------------------------------------------------------------------#
#                                       PHYTOPLANKTON PARAMS 
#------------------------------------------------------------------------------------------------------------
CMp = get_matrix("CMp", nd, nb, nn, np, nz)
K_I = 10.0         
e_o = 150/16   

Fg_p = [0.1, 0.25, 0.5, 0.68, 0.79, 0.91]      # fraction of proteome optimized to growth
Fa_p = 1. .- Fg_p                              # fraction optimized to substrate affintiy

vmax_i = [0.5, 2.0, 4.0, 7.0, 10.0, 15.0]      # max growth rates (per day)
Kp_i = vmax_i./10                              # half saturation of P_i

vmax_ij = set_vmax_ij(nn, np, vmax_i, Fg_p)    # growth rate of P_i on N
Kp_ij = set_Kp_ij(nn, np, Fa_p, CMp, vmax_ij)  # half saturation of P_i on N


# -----------------------------------------------------------------------------------------------------------#
#                                    HETEROTROPHIC BACTERIA PARAMS
#------------------------------------------------------------------------------------------------------------#
CM = [1  0  0  0  0  0  0  0  0  0  0  0  0   
      0  1  0  0  0  0  0  0  0  0  0  0  0    
      0  0  1  0  0  0  0  0  0  0  0  0  0  
      0  0  0  1  0  0  0  0  1  0  0  0  0 
      0  0  0  0  1  0  0  0  0  1  0  0  0 
      0  0  0  0  0  1  0  0  0  0  1  0  0 
      0  0  0  0  0  0  1  0  0  0  0  1  0 
      0  0  0  0  0  0  0  1  0  0  0  0  1 ] 
    
y_i = ones(nd)*0.3
y_ij = broadcast(*, y_i, CM) 
yo_ij = y_ij*10                                # PLACEHOLDER VALUE mol B/mol O2. not realistic

Fg_b = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
Fa_b = 1. .- Fg_b      

umax_i = [1., 1., 1., 1., 1., 1., 1., 1.]      # Fg_b and umax_i are dummy vals - already provided from Emily's previous work (trade-off applied)
Km_i = umax_i./10 

umax_ij =  [10.0  0  0  0  0  0  0  0  0  0  0  0  0    # first 3 are POM, next 10 are DOM
            0  1.0  0  0  0  0  0  0  0  0  0  0  0    
            0  0  0.1  0  0  0  0  0  0  0  0  0  0  
            0  0  0  32.0  0  0  0  0  22.0  0  0  0  0 
            0  0  0  0  5.6  0  0  0  0  4.0  0  0  0 
            0  0  0  0  0  1.0  0  0  0  0  0.71  0  0 
            0  0  0  0  0  0  0.18  0  0  0  0  0.13  0 
            0  0  0  0  0  0  0  0.029  0  0  0  0  0.022 ] 

Km_ij =[0.54  0  0  0  0  0  0  0  0  0  0  0  0    # in units of K_DON, calculated from K_DOC in zakem paper by dividing by 5
        0  0.054  0  0  0  0  0  0  0  0  0  0  0    
        0  0  0.0054  0  0  0  0  0  0  0  0  0  0  
        0  0  0  1.54  0  0  0  0  0.7  0  0  0  0 
        0  0  0  0  0.28  0  0  0  0  0.124  0  0  0 
        0  0  0  0  0  0.048  0  0  0  0  0.022  0  0 
        0  0  0  0  0  0  0.0086  0  0  0  0  0.004  0 
        0  0  0  0  0  0  0  0.00154  0  0  0  0  0.0007 ]  



# -----------------------------------------------------------------------------------------------------------#
#                                       ZOOPLANKTON PARAMS 
#------------------------------------------------------------------------------------------------------------#
if graze == 1
# GrM - first 6 cols are phyto, next 3 are POM and then 10 dom consuming bacteria, rows are zoo
    GrM =  [1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0 ] ;

    g_max = ones(nz)*1.6
    K_g = ones(nz)*1.0
    γ = ones(nz)*0.3
    zIC = ones(Float64, ngrid, nz) * 0.01
else
    GrM = fill(0, (nz, (np + nb)))
    g_max = zeros(nz)
    K_g = zeros(nz)
    γ = zeros(nz)
    zIC = zeros(Float64, ngrid, nz) 

end

# -----------------------------------------------------------------------------------------------------------#
#                                       VIRUS PARAMS 
#------------------------------------------------------------------------------------------------------------#
if lysis == 1
    # rows are viruses, cols are bacteria
    VM = [  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0 
            0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  ];
        
    # vly = 4*10e−6         # lysis rate 8*10e-12 L/virus/day (weitz et al 2015) or 2.17 * 10^-11 (2.17*10e−14 m3/virus/day - Xie et al 2022)
    vly = 1.6
    vbs = 23                # burst size 
    vde = 0.1               # decay rate (day-1)  0.17
    vIC = ones(Float64, ngrid, nv) * 0.1

    # how much nitrogen per virus? get a virus quota per nitrogen i.e. mmol N / virus then convert back to mmol/day
else 
    VM = fill(0, (nv, nb))
    vly = 0.0
    vbs = 0                
    vde = 0.0   
    vIC = zeros(Float64, ngrid, nv)       
end


# -----------------------------------------------------------------------------------------------------------#
#                                     MORTALITY RATES (mmol/m3/day)
#------------------------------------------------------------------------------------------------------------#
m_lp = ones(np) * 0.01 
m_qp = ones(np) * 0.1  # (.1 if explicit grazers, if not, 1)

m_lb = ones(nb) * 0.01
m_qb = ones(nb) * 0.05

m_lz = ones(nz) * 0.01
m_qz = ones(nz) * 0.5


# -----------------------------------------------------------------------------------------------------------#
#                                           ORGANIC MATTER
#------------------------------------------------------------------------------------------------------------#
om_dist_mort = [0.06, 0.19, 0.06, 0.02, 0.13, 0.39, 0.13, 0.02 ]  
om_dist_lys = [0.0, 0.0, 0.0, 0.4, 0.3, 0.1, 0.0, 0.2]   # viral shunt reduces POM, enhances lab, some slab and ref from cell walls
om_dist_vde = [0.0, 0.0, 0.0, 0.7, 0.3, 0.0, 0.0, 0.0]  # viral decay mostly amino acids/labile DOM


# Sinking rate for POM  
ws = zeros(nd)                  
ws[1], ws[2], ws[3] = 10.0, 10.0, 10.0


#------------------------------------------------------------------------------------------------------------#
#                                          PHYSICAL ENVIRONMENT
#------------------------------------------------------------------------------------------------------------#
# VERTICAL VELOCITY
    w = zeros(ngrid + 1)           
    wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
    wd[1,:] .= 0                    # no flux boundary at surface 
    wd[end,:] .= 0                  # no flux boundary (bottom box accumulates D)

# VERTICAL MIXING 
    if season == 1 
        mlz = 30                    # Mixing lengthscale 
        kappazmin = 1e-4            # min mixing coeff -value for most of the deep ocean (higher at top and bottom)
    else
        mlz = 15
        kappazmin = 1e-5
    end
                     
    kappazmax = 1e-2              # max mixing coeff -value at top of mixed layer (and bottom boundary mixed layer)
    kappa_z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
    kappa_z[1] = 0
    kappa_z[end] = 0

# LIGHT (Irradiance, I) 
    euz = 25                    # euphotic zone lengthscale (m)
    light_top = 700             # avg incoming PAR = total/2  = 1400 kW/m2/year / 2 ??! Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
    light = light_top .* exp.( -zc ./ euz)

# OXYGEN (air-sea exchange)
    o2_sat = 212.1              # mmol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
    Kgast = 3e-5*86400          # m/d
    ml_boxes = 100/dz           # discrete n of boxes in the mixed layer, close to 100m total sum
    koverh = Kgast/ml_boxes     # gas transfer coeff for each of the n boxes comprising the ml. 

# OXYGEN (deep oxygen relaxation)
    o2_deep = 200.0             # mmol/m3, avg (for ~7 C) and 35
    t_o2relax = 0.01            # 1/day, range from 0.01 to 0.1. Set to 0 to turn off.


# TEMPERATURE (SPOT along water column)
#fit to SPOT data (approx 20 to 4, approx 16 to 4)
    if season == 1 
        temp = 6.5 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 3
    else
        temp = 10 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 2.9
    end

# TEMPERATURE (modification to metabolic rates)
    temp_coeff_arr = 0.8
    temp_ae_arr = -4000
    temp_ref_arr = 293.15   
    t_kel = 273.15
    temp_fun = temp_coeff_arr .* exp.(temp_ae_arr .*(1 ./ (temp .+ t_kel) .- 1 ./ temp_ref_arr))


# -----------------------------------------------------------------------------------------------------------#
#                                   INITIAL CONDITIONS (mmol N/m3)
#------------------------------------------------------------------------------------------------------------# 
nIC = ones(Float64, ngrid, nn) * 20.0
pIC = ones(Float64, ngrid, np) * 0.01 
dIC = ones(Float64, ngrid, nd) * 0.1
bIC = ones(Float64, ngrid, nb) * 0.01 
oIC = ones(Float64, ngrid, 1)  * 100.0


# -----------------------------------------------------------------------------------------------------------#
#                                   INSTANTIATE PARAMS 
#------------------------------------------------------------------------------------------------------------#
params = Prms(
            days, dt, nt, nrec, H, dz, np, nb, nz, nn, nd, nv, pIC, bIC, zIC, nIC, dIC, vIC, oIC, 
            vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, light, temp_fun, K_I, CMp, Fg_p,
            umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
            g_max, K_g, γ, m_lz, m_qz, GrM, kappa_z, wd, ngrid, pulse, 
            om_dist_mort, om_dist_lys, om_dist_vde, VM, vly, vbs, vde,
            e_o, yo_ij, koverh, o2_sat, ml_boxes, t_o2relax, o2_deep, fsaven, lysis, graze
        )

    


