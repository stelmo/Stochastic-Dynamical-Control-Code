# Plot the Linear Model KF MC results
using PyPlot
using KernelDensity

# mcN = 50
# include("lin_mod_kf_lin_mpc_mean_mc.jl")
# include("lin_mod_kf_lin_mpc_var_conf_90_mc.jl")
# include("lin_mod_kf_lin_mpc_var_conf_99_mc.jl")
# include("lin_mod_kf_lin_mpc_var_conf_999_mc.jl")

mc1 = readcsv("linmod_kf_mean_mc2.csv")
mc2 = readcsv("linmod_kf_var90_mc2.csv")
mc3 = readcsv("linmod_kf_var99_mc2.csv")
mc4 = readcsv("linmod_kf_var999_mc2.csv")


rows, cols = size(mc1) # all will have the same dimension
ts = [0.0:0.1:80]

# Now plot 90 % confidence regions!
rc("text", usetex=true)
rc("font", family="serif", serif="Computer Modern", size=32)
figure()
subplot(4, 1, 1) # mean
for k=1:cols
  plot(ts, mc1[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (I)")
locator_params(nbins=4)

subplot(4, 1, 2) # 90%
for k=1:cols
  plot(ts, mc2[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (II)")
locator_params(nbins=4)

subplot(4, 1, 3) # 99%
for k=1:cols
  plot(ts, mc3[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (III)")
locator_params(nbins=4)

subplot(4, 1, 4) # 99.9%
for k=1:cols
  plot(ts, mc4[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (IV)")
locator_params(nbins=4)
xlabel("Time [min]")

mcerr1 = 0
for k=1:cols
  mcerr1 +=  abs(Results.calcError3(mc1[end-100:end, k], ysp+b[1]))
end
println("The average MC error is:", mcerr1/cols)

mcerr2 = 0
for k=1:cols
  mcerr2 +=  abs(Results.calcError3(mc2[end-100:end, k], ysp+b[1]))
end
println("The average MC error is:", mcerr2/cols)

mcerr3 = 0
for k=1:cols
  mcerr3 +=  abs(Results.calcError3(mc3[end-100:end, k], ysp+b[1]))
end
println("The average MC error is:", mcerr3/cols)

mcerr4 = 0
for k=1:cols
  mcerr4 +=  abs(Results.calcError3(mc4[end-100:end, k], ysp+b[1]))
end
println("The average MC error is:", mcerr4/cols)
