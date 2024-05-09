
# using Parameters

# @autodoc

@with_kw struct Prms
    days::Int64                    # num days (run time)
    dt::Float64                    # length of one time-step
    nt::Int64                      # num days / num timesteps
    nrec::Int64                    # num timepoints to record
    H::Int64                       # total depth (m)
    dz::Int64                      # height of each box/parcel of water (m)
    np::Int64                      # number of phytoplankton
    nb::Int64                      # number of bacteria
    nz::Int64                      # number of zooplankton
    nn::Int64                      # number of inorganic matter pools
    nd::Int64                      # number of organic matter pools
    nv::Int64                      # number of virus pools
    pIC::Array{Float64,2}          # initial condition for np p at time 0
    bIC::Array{Float64,2}          # initial condition for nb b at time 0
    zIC::Array{Float64,2}          # initial condition for nz z at time 0
    nIC::Array{Float64,2}          # initial condition for nn n at time 0
    dIC::Array{Float64,2}          # initial condition for nd d at time 0
    vIC::Array{Float64,2}          # initial condition for nv v at time 0
    oIC::Array{Float64,2}          # initial condition for no o at time 0
    vmax_i::Array{Float64,1}       # max uptake rate overall for p_i
    vmax_ij::Array{Float64,2}      # max uptake rate of phyto j on n
    Kp_i::Array{Float64,1}         # half saturation rate for p_i 
    Kp_ij::Array{Float64,2}        # half saturation rate of phyto j on n 
    m_lp::Array{Float64,1}         # linear mort of p (1/d) 
    m_qp::Array{Float64,1}         # quadratic mort of p (m3/mmol/d)
    light::Array{Float64, 1}       # light (irradiance) W/m2
    temp_fun::Array{Float64, 1}    # temperature function (modification to metabolic rates)
    K_I::Float64                   # light half-saturation constant
    CMp::Array{Bool,2}             # consumption matrix: nn X np
    Fg_p::Array{Float64, 1}        # fraction of proteome devoted to growth for each P
    umax_i::Array{Float64,1}       # max uptake rate overall for d_i
    umax_ij::Array{Float64,2}      # max uptake rate of bacteria j on d_i
    Km_i::Array{Float64,1}         # half saturation rate overall for d_i 
    Km_ij::Array{Float64,2}        # half saturation rate of bacteria j on d_i 
    y_ij::Array{Float64,2}         # yield rate of bacteria j on d_i
    m_lb::Array{Float64,1}         # linear mort of b
    m_qb::Array{Float64,1}         # quadratic mort of b
    CM::Array{Bool,2}              # consumption matrix: nd X nb
    Fg_b::Array{Float64, 1}        # fraction of proteome devoted to growth for each B
    g_max::Array{Float64,1}        # max grazing rate of z
    K_g::Array{Float64,1}          # half saturation rate of z
    γ::Array{Float64,1}            # fraction of assimilation (assimilation efficiency) for z 
    m_lz::Array{Float64,1}         # linear mort of z
    m_qz::Array{Float64,1}         # quadratic mort of z
    GrM::Array{Bool,2}             # grazing matrix: nz x (np + nb)
    kappa_z::Array{Float64, 1}     # vertical eddy diffusivities
    wd::Array{Float64, 2}          # sinking rate for POM
    ngrid::Int64                   # number of boxes at depth
    pulse::Int64                   # 1 = no nutrient pulsing, 2 = periodic pulse, 3 = semi-stochastc pulse
    om_dist_mort::Array{Float64, 1} # distribution of OM to each d pool from mortality
    om_dist_lys::Array{Float64, 1} # distribution of OM to each d pool from viral lysis
    om_dist_vde::Array{Float64, 1} # distribution of OM to each d pool from viral decay
    VM::Array{Bool,2}              # viral lysis matrix: nv X nb
    vly::Float64                   # viral lysis rate
    vbs::Float64                   # viral burst size
    vde::Float64                   # viral decay rate
    e_o::Float64                   # production of O2 (excretion) (mol O2 per mol N uptake)
    yo_ij::Array{Float64, 2}       # O2 yield rate of bacteria j on d_i (mol B/mol O2)
    koverh::Float64                # gas transfer coefficient for each box comprising the mixed layer
    o2_sat::Float64                # O2 half-sat constant (mmol/m3)
    ml_boxes::Int64                # num boxes in mixed layer, close to 100m total sum
    t_o2relax::Float64             # deep oxygen relaxation (1/day)
    o2_deep::Float64               # mmol/m3
    fsaven::String                 # save file name
    lysis::Int64                   # Explicit viral lysis = 1, Implicit (quadratic mortality term) = 2
    graze::Int64                   # Include grazers = 1, Exclude grazers = 2
end

