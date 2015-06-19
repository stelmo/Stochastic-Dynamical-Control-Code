# Plot the Linear Model KF MC results
using PyPlot

mcN = 200
include("nonlin_mod_kf_lin_mpc_mean_mc.jl")
mcN = 300
include("nonlin_mod_kf_lin_mpc_var_conf_90_mc.jl")
mcN = 400
include("nonlin_mod_kf_lin_mpc_var_conf_99_mc.jl")
mcN = 500
include("nonlin_mod_kf_lin_mpc_var_conf_999_mc.jl")

# rc("font", family="serif", serif="Computer Modern", size=24)
# rc("text", usetex=true)
#
# mc1 = abs(readcsv("nonlinmod_kf_mean.csv"))
# mc2 = abs(readcsv("nonlinmod_kf_var90.csv"))
# mc3 = abs(readcsv("nonlinmod_kf_var99.csv"))
# mc4 = abs(readcsv("nonlinmod_kf_var999.csv"))
#
# # Magenta = mean
# # Red = 90%
# # Green = 99%
# # Blue = 99.9%
# cs1 = scatter(mc1[1,:][:], mc1[2,:][:], c="m")
# cs2 = scatter(mc2[1,:][:], mc2[2,:][:], c="r")
# cs3 = scatter(mc3[1,:][:], mc3[2,:][:], c="g")
# cs4 = scatter(mc4[1,:][:], mc4[2,:][:], c="b")
# xlabel("Mahalanobis Area")
# ylabel("Time in Violation")
# legend([cs1, cs2, cs3, cs4],["Expected Value",L"90$\%$ Chance", L"99$\%$ Chance", L"99.9$\%$ Chance"])
# xlim([0, maximum([mc1[1,:] mc2[1,:] mc3[1,:] mc4[1,:]])])
# ylim([0, maximum([mc1[2,:] mc2[2,:] mc3[2,:] mc4[2,:]])])
