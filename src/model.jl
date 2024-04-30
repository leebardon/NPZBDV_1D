# using Printf
# using Dates
# using NCDatasets
# using SparseArrays, LinearAlgebra


#make sum drop dimension automatically and set subnormals to zero to avoid slowdown
sumd(x,dims) = dropdims(sum(x,dims=dims),dims=dims)
set_zero_subnormals(true)


#TODO add a sinking (dsink) to remove some detritus or it wont equilibriate 
function run_NPZBDV(prms, season, lysis, bloom=false)

    trec = prms.nt/prms.nrec # frequency of recording (no. timesteps / no. recordings - usually 5)
    start_time = now()

    # Generate empty arrays
    nrec1 = Int(prms.nrec + 1)
    track_n = Array{Float64, 3}(undef, prms.ngrid, prms.nn, nrec1) 
    track_p = Array{Float64, 3}(undef, prms.ngrid, prms.np, nrec1) 
    track_z = Array{Float64, 3}(undef, prms.ngrid, prms.nz, nrec1) 
    track_b = Array{Float64, 3}(undef, prms.ngrid, prms.nb, nrec1) 
    track_d = Array{Float64, 3}(undef, prms.ngrid, prms.nd, nrec1)
    track_v = Array{Float64, 3}(undef, prms.ngrid, prms.nv, nrec1)
    track_o = Array{Float64, 3}(undef, prms.ngrid, 1, nrec1) 
    track_time = Array{Float64,1}(undef, nrec1)   
 
    track_n[:,:,1] .= prms.nIC
    track_p[:,:,1] .= prms.pIC
    track_z[:,:,1] .= prms.zIC
    track_b[:,:,1] .= prms.bIC
    track_d[:,:,1] .= prms.dIC
    track_v[:,:,1] .= prms.vIC
    track_o[:,:,1] .= prms.oIC
    track_time[1] = 0

    #--------------------------------------
    #Initial conditions at time t = 0
    ntemp = copy(prms.nIC) 
    ptemp = copy(prms.pIC) 
    ztemp = copy(prms.zIC) 
    btemp = copy(prms.bIC)
    dtemp = copy(prms.dIC) 
    vtemp = copy(prms.vIC) 
    otemp = copy(prms.oIC) 


    # Create quasi-random vector of unique timesteps on which nutrient pulses occur, if pulsing enabled
    season == 1 ? npulses = floor(prms.days/25) : npulses = floor(prms.days/45) # nutrient pulse on average every 10 days in winter, 40 days in summer
    pulse_vec = sample(1:Int(prms.nt), Int(npulses), replace=false)


    # @time for t = 1:prms.nt 
    for t = 1:prms.nt
        # Runge-Kutta 4th order 
        ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp = rk4(ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, prms, t, bloom, season, lysis)

        if mod(t, trec) == 0 #i.e. if t divisible by 5

            track_n, track_p, track_z, track_b, track_d, track_v, track_o = update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_v, track_o, 
                                                                            track_time, ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, t, trec, prms)
            
            println("Total N: ", sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp) + sum(vtemp))
            
        end 

        # Nutrient redistribution pulsing routine (every 25 or 45 days)
        if bloom == false
            if prms.pulse == 2
                if t < 36478000  # no pulse in final 2 months
                    if season == 1
                        if t % 2500 == 0 
                            println("PULSED at t=$t")
                            ntemp, dtemp = pulse_nutrients(ntemp, dtemp, prms, prms.pulse)      
                        end
                    else
                        if t % 4500 == 0 
                            println("PULSED at t=$t")
                            ntemp, dtemp = pulse_nutrients(ntemp, dtemp, prms, prms.pulse)  
                        end
                    end
                end
            end 
        end

        # Save outputs
        if t == prms.nt
            end_time = now() 
            save_full_run(track_p, track_b, track_z, track_n, track_d, track_v, track_o, track_time, start_time, end_time, prms, season, lysis)
            if bloom == false
                save_endpoints(track_p, track_b, track_z, track_n, track_d, track_v, track_o, track_time, start_time, end_time, prms, season, lysis)
            end
        end
    end 


    return ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, track_time

end 


function model_functions(N, P, Z, B, D, V, O, prms, t, bloom, season, lysis)

    if bloom == true
        if t == 90000
            D[:,:,:] .+= 1.0
        end
        if 90000 < t < 104000
            season == 1 ? mlz = 80 : mlz = 60 
            zf = [0 : 10 : H;] 
            kappazmin = 1e-4              # min mixing coeff -value for most of the deep ocean (higher at top and bottom)
            kappazmax = 1e-2              # max mixing coeff -value at top of mixed layer (and bottom boundary mixed layer)
            kappa_Z = (kappazmax .* exp.(-zf/mlz) .+ kappazmin .+ kappazmax .* exp.((zf .- H) / 100.)) .* 3600 .* 24 
            kappa_Z[1] = 0
            kappa_Z[end] = 0
            # kappa_Z = prms.kappa_z
        else
            kappa_Z = prms.kappa_z
        end
    else
        kappa_Z = prms.kappa_z
    end

    #Transport
    dNdt = diffusion(N, kappa_Z, prms.dz)
    dPdt = diffusion(P, kappa_Z, prms.dz)
    dZdt = diffusion(Z, kappa_Z, prms.dz)
    dBdt = diffusion(B, kappa_Z, prms.dz)
    dDdt = diffusion(D, kappa_Z, prms.dz) .- advection(D, prms.wd, prms.dz)
    dVdt = diffusion(V, kappa_Z, prms.dz)
    dOdt = diffusion(O, kappa_Z, prms.dz)

    d_gain_mort = zeros(prms.ngrid)
    d_gain_vly = zeros(prms.ngrid)
    d_gain_vde = zeros(prms.ngrid)

    # phyto uptake 
    dPdt, dNdt, dOdt = phyto_uptake(prms, N, P, dNdt, dPdt, dOdt, t)

    # bacteria uptake
    dDdt, dBdt, dNdt, dOdt = bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, dOdt)

    # grazing
    dZdt, dNdt, dPdt, dBdt = grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt)

    #phytoplankton mortality 
    dPdt, d_gain_mort = phyto_mortality(prms, P, dPdt, d_gain_mort)

    #bacterial mortality 
    dBdt, d_gain_mort = bacterial_mortality(prms, B, dBdt, d_gain_mort, lysis)
    
    #zooplankton mortality 
    dZdt, d_gain_mort = zoo_mortality(prms, Z, dZdt, d_gain_mort)

    # viral lysis
    if lysis == 1
        # viral lysis (B only for now)
        dVdt, dBdt, d_gain_vly = viral_lysis(prms, B, V, dVdt, dBdt, d_gain_vly)
        #viral decay
        dVdt, d_gain_vde = viral_decay(prms, V, dVdt, d_gain_vde)
    else
    end

    #split accumulated OM into nd pools according to probability of generation
    dDdt = distribute_d(prms, dDdt, d_gain_mort, d_gain_vly, d_gain_vde)

    #oxygen
    dOdt = change_in_o2(prms, O, dOdt)


    return dNdt, dPdt, dZdt, dBdt, dDdt, dVdt, dOdt

end 


function phyto_uptake(prms, N, P, dNdt, dPdt, dOdt, t)
    # uptake calculated as function of most limiting factor - light supply ot nutrient supply

    II, JJ = get_nonzero_axes(prms.CMp)

    Iz = calc_light_attenuation(P, prms, t)

    for j = axes(II, 1)
        uptake = P[:,JJ[j]] .* prms.temp_fun .* prms.vmax_ij[II[j],JJ[j]] .* min.(N ./ (N .+ prms.Kp_ij[II[j],JJ[j]]), Iz ./ (Iz .+ prms.K_I))
        dNdt += -uptake
        dOdt += uptake * prms.e_o
        dPdt[:,JJ[j]] += uptake 
    end

    return dPdt, dNdt, dOdt

end 


function calc_light_attenuation(P, prms, t)
    # Following Zakem et al 2015

    I_max = 1400                                        # Incident radiation at surface W/m2 #TODO what's the val at SPOT?
    # I_in = I_max/2                                    # avg incoming PAR = (1400/2)  Light_avg*(cos(t*dt*2*3.1416)+1) for light daily cycle
    I_in = I_max/2*(cos(t*prms.dt*2*3.1416)+1)          # for daily light cycling
    a_chlD = 0.04                                       # Absorption coeff incorporating Chl-a and CDOM (m2/mg Chl) 
    chl2c_max = 0.2                                     # max chlorophyll to carbon ratio (mg Chl/mmol C)
    chl2c_min = 0.02                                    # max chlorophyll to carbon ratio (mg Chl/mmol C)
    kw = 0.04                                           # attenuation coeff of water (m2/mg Chl)
    zc = [prms.dz/2 : prms.dz : prms.H - prms.dz/2;]    # centered depth 

    Iz = zeros(prms.ngrid)                         
    chl_tot = sum(P, dims=2) .* chl2c_min .* 6.6       # mmolN/m3 * mgChl/mmolC  (redfield ratio -> 6.6 C for every N)

    for d in range(1, prms.ngrid)
        Iz[d] = I_in * exp(-zc[d]*(kw + sum(chl_tot[1:d]*a_chlD)))
    end

    return Iz

end


function bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, dOdt)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)
        growth_rate = prms.temp_fun .* prms.umax_ij[II[j],JJ[j]] .* D[:,II[j]] ./ (D[:,II[j]] .+ prms.Km_ij[II[j],JJ[j]])
        uptake = B[:,JJ[j]] .* growth_rate
        yield = prms.y_ij[II[j],JJ[j]]
        respired = (1 - yield)
        dDdt[:,II[j]] += -uptake
        dBdt[:,JJ[j]] += uptake .* yield
        dNdt += uptake .* respired
        dOdt +=  -uptake .* (yield ./ prms.yo_ij[II[j],JJ[j]])
    end

    return dDdt, dBdt, dNdt, dOdt

end


function grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt)

    GrM = copy(prms.GrM)
    for k = 1:prms.nz
        if sum(GrM[k,1:prms.np]) > 0 
            dZdt, dNdt, dPdt = phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k)
        end
        if sum(GrM[k,prms.np+1:end]) > 0 
            dZdt, dNdt, dBdt = bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k)
        end
    end

    return dZdt, dNdt, dPdt, dBdt

end

        function phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k)

            prey = sum(GrM[k,1:prms.np]' .*P, dims=2)
            gp = prms.g_max[k] .* prey ./ (prey .+ prms.K_g[k])
            dZdt[:,k] += prms.γ[k] .* gp .* Z[:,k]
            dNdt += (1 - prms.γ[k]) .* gp .* Z[:,k]
            dPdt += -gp .* Z[:,k] .* GrM[k,1:prms.np]' .* P ./ prey 
            
            return dZdt, dNdt, dPdt

        end

        function bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k)

            prey = sum(GrM[k,prms.np+1:end]' .*B, dims=2)
            gb = prms.g_max[k] .* prey ./ (prey .+ prms.K_g[k])
            dZdt[:,k] += prms.γ[k] .* gb .* Z[:,k]
            dNdt += (1 - prms.γ[k]) .* gb .* Z[:,k]
            dBdt +=  -gb .* Z[:,k] .* GrM[k,prms.np+1:end]' .* B ./ prey

            return dZdt, dNdt, dBdt

        end


function viral_lysis(prms, B, V, dVdt, dBdt, d_gain_vly)

    II, JJ = get_nonzero_axes(prms.VM)
        
    for j = axes(II, 1)
        lysis_Bi = prms.vly .* VM[II[j], JJ[j]] .* B[:, JJ[j]] .* V[:, II[j]]
        v_growth = lysis_Bi .* 0.3
        d_gain_vly += lysis_Bi * 0.7
        dBdt[:, JJ[j]]  += -lysis_Bi
        dVdt[:, II[j]] += v_growth
    end
    # NOTE as a first approximation, 30% of the lysed B goes into V growth, and 70% is returned to D 
    # leaving out burst size for now as it was just killing all B fast
        
    return dVdt, dBdt, d_gain_vly
        
end


function phyto_mortality(prms, P, dPdt, d_gain_mort)

    pmort = (transpose(prms.m_lp) .+ transpose(prms.m_qp) .* P) .* P
    dPdt += -pmort
    d_gain_mort += sum(pmort, dims=2)

    return dPdt, d_gain_mort

end


function bacterial_mortality(prms, B, dBdt, d_gain_mort, lysis)

    if lysis == 1
        bmort = transpose(prms.m_lb) .* B
    else
        bmort = (transpose(prms.m_lb) .+ transpose(prms.m_qb) .* B) .* B
    end

    dBdt += -bmort
    d_gain_mort += sum(bmort, dims=2)

    return dBdt, d_gain_mort

end


function zoo_mortality(prms, Z, dZdt, d_gain_mort)

    zmort = (transpose(prms.m_lz) .+ transpose(prms.m_qz) .* Z) .* Z
    dZdt += -zmort
    d_gain_mort += sum(zmort, dims=2)

    return dZdt, d_gain_mort

end


function viral_decay(prms, V, dVdt, d_gain_vde)

    decay = transpose(prms.vde) .* V
    dVdt += -decay
    d_gain_vde += sum(decay, dims=2)

    return dVdt, d_gain_vde

end


function change_in_o2(prms, O, dOdt)

    dOdt[1:Int(prms.ml_boxes)] += prms.koverh .* (prms.o2_sat .- O[1:Int(prms.ml_boxes)])
    dOdt += prms.t_o2relax .* (prms.o2_deep .- O) #relaxation at depth (lateral flux)

    return dOdt

end


function distribute_d(prms, dDdt, d_gain_mort, d_gain_vly, d_gain_vde)

    dDdt += d_gain_mort .* transpose(prms.om_dist_mort)
    dDdt += d_gain_vly .* transpose(prms.om_dist_lys)
    dDdt += d_gain_vde .* transpose(prms.om_dist_vde)

    return dDdt

end


function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


function update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_v, track_o, track_time, 
    ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, t, trec, prms)

    j = Int(t÷trec + 1)
    t_id = t.*prms.dt
    track_p[:,:,j] .= ptemp
    track_b[:,:,j] .= btemp 
    track_z[:,:,j] .= ztemp 
    track_n[:,:,j] .= ntemp 
    track_d[:,:,j] .= dtemp
    track_v[:,:,j] .= vtemp
    track_o[:,:,j] .= otemp
    track_time[j] = t_id 

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.days, t_id/prms.days*100, now())

    return track_n, track_p, track_z, track_b, track_d, track_v, track_o, track_time

end


