
function advection(c, wd, dz)
    c = vcat(zeros(2,size(c)[2]), c, zeros(2, size(c)[2]))
    positive_wd = (wd + abs.(wd))./2. 
    negative_wd = (wd - abs.(wd))./2.
    f_up = positive_wd[2:end,:].*c[3:end-2,:] + negative_wd[2:end,:].*c[4:end-1,:]
    f_down = positive_wd[1:end-1,:].*c[2:end-3,:] + negative_wd[1:end-1,:].*c[3:end-2,:]
    adv = (f_up .- f_down)./dz 
    return adv 
end


function diffusion(c, kappa_z, dz)
    c = vcat(transpose(zeros(size(c)[2])), c, transpose(zeros(size(c)[2])))
    kappa_z_len = size(kappa_z)[1]
    f_up = kappa_z[1:kappa_z_len-1] .* (c[2:kappa_z_len,:] .- c[1:kappa_z_len-1,:]) ./dz 
    f_down = kappa_z[2:kappa_z_len] .* (c[3:kappa_z_len+1,:] .- c[2:kappa_z_len,:]) ./dz
    diff = (f_down .- f_up) ./ dz
    return diff
end


# function get_pom_sinking_rates(npom, nd)
#     # particle size follows power law distribution, but since we are dealing with pools, 1 large may contain e.g. same volume as 100 small
#     # sinking rate range informed by fig 3 of Omand et al 2020 as guide https://www.nature.com/articles/s41598-020-60424-5
#     # here we use 1 m/d for sinking rate of smallest particle sizes and 100 m/d as large particles / small aggregates
  
#     ws = zeros(nd)
#     dist = (10).^(range(0.1, stop=2, length=1000))    
    
#     if npom == 1
#         ws[1] = median(dist)
#     elseif npom == 2
#         ws[1] = minimum(dist)
#         ws[2] = median(dist)
#     elseif npom == 3
#         ws[1] = minimum(dist)
#         ws[2] = median(dist)
#         ws[3] = maximum(dist)
#     else
#         ws[1] = popfirst!(dist)
#         ws[npom] = pop!(dist)

#         # sample remaining pom sinking speeds evenly from IQR
#         deleteat!(dist, findall(x -> x<quantile(dist, 0.25) || x>quantile(dist, 0.75), dist))
#         smple = LinRange(minimum(dist), maximum(dist), npom-2)

#         for i in range(1, npom)
#             if i == 1 || i == npom
#                 :nothing
#             else
#                 ws[i] = smple[i-1]
#             end
#         end

#     end

#     return ws

# end 


# function get_physical_params(ngrid, nd, season, zc, zf, H, dz)

#     # WATER VERTICAL VELOCITY
#         w = zeros(ngrid + 1)            
#         wd = transpose(repeat(ws, 1, ngrid + 1)) + repeat(w, 1, nd) # ngrid+1 x nd
#         wd[1,:] .= 0                        # no flux boundary at surface 
#         wd[end,:] .= 0                      # no flux boundary (bottom box accumulates D)

#     # VERTICAL MIXING 
#         season == 1 ? mlz = 25 : mlz = 15   # mixed layer lengthscale
#         kappazmin = 1e-4                    # min mixing coeff -value for most of the deep ocean (higher at top and bottom)
#         kappazmax = 1e-2                    # max mixing coeff -value at top of mixed layer (and bottom boundary mixed layer)
#         kappa_z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
#         kappa_z[1] = 0
#         kappa_z[end] = 0

#     # LIGHT (Irradiance, I) 
#         euz = 25                            # euphotic zone lengthscale #NOTE scales with amount of biomass but VERY sensitive
#         light_top = 700                     # avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
#         light = light_top .* exp.( -zc ./ euz)

#     # OXYGEN (air-sea exchange)
#         o2_sat = 212.1                      # mmol/m3 from calc_oxsat(25+273,35) in matlab. WOCE clim-avg surf T at 10S, E. Pac.
#         Kgast = 3e-5*86400                  # m/d
#         ml_boxes = 100/dz                   # discrete n of boxes in the mixed layer, close to 100m total sum
#         koverh = Kgast/ml_boxes             # gas transfer coeff for each of the n boxes comprising the ml. 

#     # OXYGEN (deep oxygen relaxation)
#         o2_deep = 200.0                     # mmol/m3, avg (for ~7 C) and 35
#         t_o2relax = 0.01                    # 1/day, range from 0.01 to 0.1. Set to 0 to turn off.


#     # TEMPERATURE (SPOT along water column - fit to SPOT data (approx 20 to 4 summer, approx 16 to 4 winter))
#         if season == 1 
#             temp = 6.5 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 3
#         else
#             temp = 10 .* exp.(-zc ./ 150) .+ 9 .* exp.(-zc ./ 500) .+ 2.9
#         end

#     # TEMPERATURE (modification to metabolic rates)
#         temp_coeff_arr = 0.8
#         temp_ae_arr = -4000
#         temp_ref_arr = 293.15   
#         t_kel = 273.15
#         temp_fun = temp_coeff_arr .* exp.(temp_ae_arr .*(1 ./ (temp .+ t_kel) .- 1 ./ temp_ref_arr))


#     return wd, mlz, kappa_z, light, o2_sat, ml_boxes, koverh, o2_deep, t_o2relax, temp_fun

# end