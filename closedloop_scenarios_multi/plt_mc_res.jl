# Plot the Linear Model KF MC results
using PyPlot
using KernelDensity

rc("font", family="serif", serif="Computer Modern", size=24)
rc("text", usetex=true)

mc1 = abs(readcsv("spf_mean.csv"))
mc2 = abs(readcsv("spf_var90.csv"))
mc3 = abs(readcsv("spf_var99.csv"))

mc1 = Auxiliary.removeOutliers(mc1, 3)
mc2 = Auxiliary.removeOutliers(mc2, 3)
mc3 = Auxiliary.removeOutliers(mc3, 3)

mmc1 = mean(mc1, 2)
mmc2 = mean(mc2, 2)
mmc3 = mean(mc3, 2)

cmc1 = cov(mc1')
cmc2 = cov(mc2')
cmc3 = cov(mc3')

# Now plot 90 % confidence regions!
rc("text", usetex=true)
rc("font", family="serif", serif="Computer Modern", size=24)
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

axis(ymin=0.0, xmin=0.0)
#
# # Magenta = mean
# # Red = 90%
# # Green = 99%
# # Blue = 99.9%
xlabel("Mahalanobis area in violation")
ylabel("Time in violation [min]")
legend([cs1, cs2, cs3],["Expected Value Constraint",L"90$\%$ Chance Constraint", L"99$\%$ Chance Constraint"], loc="best")
