# Plot the Linear Model KF MC results
using PyPlot
using KernelDensity

# mcN = 50
# include("nonlin_mod_pf_lin_mpc_mean_mc.jl")
# include("nonlin_mod_pf_lin_mpc_var_conf_90_mc.jl")

mc1 = readcsv("nonlinmod_pf_mean_mc2.csv")
mc2 = readcsv("nonlinmod_pf_var90_mc2.csv")

rows, cols = size(mc1) # all will have the same dimension
ts = [0.0:0.1:80]


# Now plot 90 % confidence regions!
rc("text", usetex=true)
rc("font", family="serif", serif="Computer Modern", size=32)
figure()
subplot(2, 1, 1) # mean
for k=1:cols
  plot(ts, mc1[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (I)")
locator_params(nbins=4)

subplot(2, 1, 2) # 90%
for k=1:cols
  plot(ts, mc2[:, k], "k-", linewidth=0.5)
end
plot(ts, ones(rows)*0.49, "g-", linewidth=3.0)
ylabel(L"C$_A$ (II)")
locator_params(nbins=4)

mcerr1 = 0
for k=1:cols
  mcerr1 +=  abs(Results.calcError3(mc1[:, k], ysp+b[1]))
end
println("The average MC error is:", mcerr1/cols)

mcerr2 = 0
for k=1:cols
  mcerr2 +=  abs(Results.calcError3(mc2[:, k], ysp+b[1]))
end
println("The average MC error is:", mcerr2/cols)
