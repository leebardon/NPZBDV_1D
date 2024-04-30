
function pulse_nutrients(track_n, track_d, prms, pulse_type)

    pulsed_n, pulsed_d = redistributed_mean(track_n, track_d, prms)

    return pulsed_n, pulsed_d

end


function redistributed_mean(track_n, track_d, prms)

    pulsed_n = Array{Float64, 2}(undef, prms.ngrid, prms.nn) 
    grid_mean_n = mean(track_n)
    pulsed_n .= grid_mean_n

    pulsed_d = Array{Float64, 2}(undef, prms.ngrid, prms.nd) 
    for d in 1:prms.nd
        grid_mean_d = mean(track_d[:, d])
        pulsed_d[:, d] .= grid_mean_d
    end

    return pulsed_n, pulsed_d

end


# function weighted_mean(track_n, track_d, prms)
#NOTE - step change is too abrupt, not realistic
#     # 50% of total goes to top 200m, the rest to bottom 690m
#     pulsed_n = Array{Float64, 2}(undef, prms.ngrid, prms.nn)  
#     total_n = sum(track_n)
#     half_n = total_n/2
#     pulsed_n[1:20, :] .= half_n/20
#     pulsed_n[21:end, :] .= half_n/(prms.ngrid - 20)

#     pulsed_d = Array{Float64, 2}(undef, prms.ngrid, prms.nd)
#     for d in 1:nd
#         total_d = sum(track_d[:, d])
#         half_d = total_d/2
#         pulsed_d[1:20, d] .= half_d/20
#         pulsed_d[21:end, d] .= half_d/(prms.ngrid - 20)
#     end

#     return pulsed_n, pulsed_d
# end

