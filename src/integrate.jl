
function rk4(ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, prms, t, bloom, season, lysis)

    dNdt1, dPdt1, dZdt1, dBdt1, dDdt1, dVdt1, dOdt1 = model_functions(ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp, prms, t, bloom, season, lysis)
    
    track_n1 = ntemp .+ prms.dt/2 .* dNdt1
    track_p1 = ptemp .+ prms.dt/2 .* dPdt1
    track_z1 = ztemp .+ prms.dt/2 .* dZdt1
    track_b1 = btemp .+ prms.dt/2 .* dBdt1
    track_d1 = dtemp .+ prms.dt/2 .* dDdt1
    track_v1 = vtemp .+ prms.dt/2 .* dVdt1
    track_o1 = otemp .+ prms.dt/2 .* dOdt1

    N1tot = sum(track_p1) + sum(track_b1) + sum(track_z1) + sum(track_n1) + sum(track_d1)
 
    
    dNdt2, dPdt2, dZdt2, dBdt2, dDdt2, dVdt2, dOdt2 = model_functions(track_n1, track_p1, track_z1, track_b1, track_d1, track_v1, track_o1, prms, t, bloom, season, lysis)

    track_n2 = ntemp .+ prms.dt/2 .* dNdt2
    track_p2 = ptemp .+ prms.dt/2 .* dPdt2
    track_z2 = ztemp .+ prms.dt/2 .* dZdt2
    track_b2 = btemp .+ prms.dt/2 .* dBdt2
    track_d2 = dtemp .+ prms.dt/2 .* dDdt2
    track_v2 = vtemp .+ prms.dt/2 .* dVdt2
    track_o2 = otemp .+ prms.dt/2 .* dOdt2

    N2tot = sum(track_p2) + sum(track_b2) + sum(track_z2) + sum(track_n2) + sum(track_d2)

    
    dNdt3, dPdt3, dZdt3, dBdt3, dDdt3, dVdt3, dOdt3  = model_functions(track_n2, track_p2, track_z2, track_b2, track_d2, track_v2, track_o2, prms, t, bloom, season, lysis)

    track_n3 = ntemp .+ prms.dt .* dNdt3
    track_p3 = ptemp .+ prms.dt .* dPdt3
    track_z3 = ztemp .+ prms.dt .* dZdt3
    track_b3 = btemp .+ prms.dt .* dBdt3
    track_d3 = dtemp .+ prms.dt .* dDdt3
    track_v3 = vtemp .+ prms.dt .* dVdt3
    track_o3 = otemp .+ prms.dt .* dOdt3

    N3tot = sum(track_p3) + sum(track_b3) + sum(track_z3) + sum(track_n3) + sum(track_d3)

    
    dNdt4, dPdt4, dZdt4, dBdt4, dDdt4, dVdt4, dOdt4 = model_functions(track_n3, track_p3, track_z3, track_b3, track_d3, track_v3, track_o3, prms, t, bloom, season, lysis)

    ntemp .+= (dNdt1 .+ 2 .* dNdt2 .+ 2 .* dNdt3 .+ dNdt4) .* (prms.dt / 6)
    ptemp .+= (dPdt1 .+ 2 .* dPdt2 .+ 2 .* dPdt3 .+ dPdt4) .* (prms.dt / 6)
    ztemp .+= (dZdt1 .+ 2 .* dZdt2 .+ 2 .* dZdt3 .+ dZdt4) .* (prms.dt / 6)
    btemp .+= (dBdt1 .+ 2 .* dBdt2 .+ 2 .* dBdt3 .+ dBdt4) .* (prms.dt / 6)
    dtemp .+= (dDdt1 .+ 2 .* dDdt2 .+ 2 .* dDdt3 .+ dDdt4) .* (prms.dt / 6)
    vtemp .+= (dVdt1 .+ 2 .* dVdt2 .+ 2 .* dVdt3 .+ dVdt4) .* (prms.dt / 6)
    otemp .+= (dOdt1 .+ 2 .* dOdt2 .+ 2 .* dOdt3 .+ dOdt4) .* (prms.dt / 6)
     
    Ntot = sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp) + sum(vtemp)


    return ntemp, ptemp, ztemp, btemp, dtemp, vtemp, otemp

end