# Plot the Linear Model KF MC results
using PyPlot
using KernelDensity

include("LQG_SPF_P2_mc.jl")
include("MPC_SPF_MEAN_mc.jl")
include("MPC_SPF_VAR_conf90_mc.jl")
include("MPC_SPF_VAR_conf99_mc")

mc1 = readcsv("spf_lqg_mc2.csv")
mc2 = readcsv("spf_mean_mc2.csv")
mc3 = readcsv("spf_var90_mc2.csv")
mc4 = readcsv("spf_var99_mc2.csv")

rows, cols = size(mc1) # all will have the same dimension
ts = [0.0:0.1:200]

# Now plot 90 % confidence regions!
rc("text", usetex=true)
rc("font", family="serif", serif="Computer Modern", size=24)
figure()
subplot(4, 1, 1) # mean
p1 = plot(ts, mc1[:, 1], "k-", linewidth=0.5)
for k=2:cols
  plot(ts, mc1[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ kmol.m$^{-3}$")
legend([p1],["LQG"], loc="best")

subplot(4, 1, 2) # 90%
p2 = plot(ts, mc2[:, 1], "k-", linewidth=0.5)
for k=2:cols
  plot(ts, mc2[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ kmol.m$^{-3}$")
legend([p2],["Expected Value"], loc="best")

subplot(4, 1, 3) # 99%
p3 = plot(ts, mc3[:, 1], "k-", linewidth=0.5)
for k=1:cols
  plot(ts, mc3[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ kmol.m$^{-3}$")
legend([p3],[L"90$\%$ Chance"], loc="best")

subplot(4, 1, 4) # 99.9%
p4 = plot(ts, mc4[:, 1], "k-", linewidth=0.5)
for k=2:cols
  plot(ts, mc4[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ kmol.m$^{-3}$")
legend([p4],[L"99$\%$ Chance"], loc="best")
xlabel("Time [min]")
