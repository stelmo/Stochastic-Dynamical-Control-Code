# Plot the Linear Model KF MC results
using PyPlot

mcN = 200
include("nonlin_mod_pf_lin_mpc_mean_mc.jl")
mcN = 300
include("nonlin_mod_pf_lin_mpc_var_conf_90_mc.jl")

# rc("font", family="serif", serif="Computer Modern", size=24)
# rc("text", usetex=true)
#
# mc1 = abs(readcsv("nonlinmod_kf_mean.csv"))
# mc2 = abs(readcsv("nonlinmod_kf_var90.csv"))
# mc3 = abs(readcsv("nonlinmod_pf_mean.csv"))
# mc4 = abs(readcsv("nonlinmod_pf_var90.csv"))
#
# # Magenta = mean KF
# # Red = 90% KF
# # Green = mean PF
# # Blue = 90% PF
# cs1 = scatter(mc1[1,:][:], mc1[2,:][:], c="m")
# cs2 = scatter(mc2[1,:][:], mc2[2,:][:], c="r")
# cs3 = scatter(mc3[1,:][:], mc3[2,:][:], c="g")
# cs4 = scatter(mc4[1,:][:], mc4[2,:][:], c="b")
# xlabel("Mahalanobis Area")
# ylabel("Time in Violation")
# legend([cs1, cs2, cs3, cs4],["Expected Value KF",L"90$\%$ Chance KF", "Expected Value PF", L"90$\%$ Chance PF"])
# xlim([0, maximum([mc1[1,:] mc2[1,:] mc3[1,:] mc4[1,:]])])
# ylim([0, maximum([mc1[2,:] mc2[2,:] mc3[2,:] mc4[2,:]])])
