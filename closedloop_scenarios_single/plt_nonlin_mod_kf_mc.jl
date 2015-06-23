# Plot the Linear Model KF MC results
using PyPlot

# mcN = 500
# include("nonlin_mod_kf_lin_mpc_mean_mc.jl")
# mcN = 700
# include("nonlin_mod_kf_lin_mpc_var_conf_90_mc.jl")
# mcN = 900
# include("nonlin_mod_kf_lin_mpc_var_conf_99_mc.jl")
# mcN = 1100
# include("nonlin_mod_kf_lin_mpc_var_conf_999_mc.jl")


rc("font", family="serif", serif="Computer Modern", size=32)
rc("text", usetex=true)

mc1 = abs(readcsv("nonlinmod_kf_mean.csv"))
mc2 = abs(readcsv("nonlinmod_kf_var90.csv"))
mc3 = abs(readcsv("nonlinmod_kf_var99.csv"))
mc4 = abs(readcsv("nonlinmod_kf_var999.csv"))

mc1 = Auxiliary.removeOutliers(mc1, 3)
mc2 = Auxiliary.removeOutliers(mc2, 3)
mc3 = Auxiliary.removeOutliers(mc3, 3)
mc4 = Auxiliary.removeOutliers(mc4, 3)

mmc1 = mean(mc1, 2)
mmc2 = mean(mc2, 2)
mmc3 = mean(mc3, 2)
mmc4 = mean(mc4, 2)

cmc1 = cov(mc1')
cmc2 = cov(mc2')
cmc3 = cov(mc3')
cmc4 = cov(mc4')

# Now plot 90 % confidence regions!
a = 0.5
xs1, ys1 = Ellipse.ellipse(mmc1, cmc1)
cs1 = fill(xs1, ys1, "m", alpha=a, edgecolor="none")
plot(mmc1[1], mmc1[2], "mo", markersize=10)

xs2, ys2 = Ellipse.ellipse(mmc2, cmc2)
cs2 = fill(xs2, ys2, "r", alpha=a, edgecolor="none")
plot(mmc2[1], mmc2[2], "ro", markersize=10)

xs3, ys3 = Ellipse.ellipse(mmc3, cmc3)
cs3 = fill(xs3, ys3, "g", alpha=a, edgecolor="none")
plot(mmc3[1], mmc3[2], "go", markersize=10)

xs4, ys4 = Ellipse.ellipse(mmc4, cmc4)
cs4 = fill(xs4, ys4, "b", alpha=a, edgecolor="none")
plot(mmc4[1], mmc4[2], "bo", markersize=10)

axis(ymin=0.0, xmin=0.0)
#
# # Magenta = mean
# # Red = 90%
# # Green = 99%
# # Blue = 99.9%
xlabel("Mahalanobis area in violation")
ylabel("Time in violation [min]")
legend([cs1, cs2, cs3, cs4],["Expected value constraint",L"90$\%$ Chance constraint", L"99$\%$ Chance constraint", L"99.9$\%$ Chance constraint"], loc="best")
