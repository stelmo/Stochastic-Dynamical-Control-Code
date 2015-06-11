# Plotting and results analysis module
module Results

using PyPlot
using Ellipse
using LaTeXStrings

function plotTracking(ts, xs, ys, fmeans, us, obs, setpoint)

  tend = ts[end]
  setpoints = ones(length(ts))*setpoint

  umax = maximum(abs(us))
  if umax == 0.0
    subplt = 2
  else
    subplt = 3
  end
  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  skipmeas = int(length(ts)/20)
  skipmean = int(length(ts)/20)
  figure()
  subplot(subplt,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  if obs == 2 # plot second measurement
    y2, = plot(ts[1:skipmeas:end], ys[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k1, = plot(ts[1:skipmean:end], fmeans[1, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ksp = plot(ts, setpoints, "g-", linewidth=3)
  ylabel(L"C_A~[kmol.m^{-3}]")
  legend([x1, ksp],[L"Underlying~Model", L"Set~Point"], loc="lower right", ncol=2)
  xlim([0, tend])
  # ylim([0, 1])

  subplot(subplt,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  if obs == 1
    y2, = plot(ts[1:skipmeas:end], ys[1:skipmeas:end], "kx", markersize=5, markeredgewidth=1)
  else
    y2, = plot(ts[1:skipmeas:end], ys[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k2, = plot(ts[1:skipmean:end], fmeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel(L"T_R~[K]")
  legend([k2, y2],[L"Filtered~Mean", L"Observation"], loc="upper right", ncol=2)
  xlim([0, tend])
  # ylim([minimum(xs[2,:]), maximum(xs[2,:])])
  if subplt == 3
    subplot(subplt,1,3)
    plot(ts, us)
    xlim([0, tend])
    ylabel(L"Q~[kJ.min^{-1}]")
  end
  xlabel(L"Time~[min]")
end

function plotTracking(ts, xs, ys, fmeans, us, obs)

  tend = ts[end]

  umax = maximum(abs(us))
  if umax == 0.0
    subplt = 2
  else
    subplt = 3
  end
  rc("font", family="serif", size=24)
  rc("text", usetex=true)

  skipmeas = int(length(ts)/40)
  skipmean = int(length(ts)/40)
  figure()
  subplot(subplt,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  if obs == 2 # plot second measurement
    y2, = plot(ts[1:skipmeas:end], ys[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k1, = plot(ts[1:skipmean:end], fmeans[1, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel(L"C_A~[kmol.m^{-3}]")
  legend([x1],[L"Underlying~Model"], loc="best")
  xlim([0, tend])
  # ylim([0, 1])

  subplot(subplt,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  if obs == 1
    y2, = plot(ts[1:skipmeas:end], ys[1:skipmeas:end], "kx", markersize=5, markeredgewidth=1)
  else
    y2, = plot(ts[1:skipmeas:end], ys[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k2, = plot(ts[1:skipmean:end], fmeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel(L"T_R~[K]")
  legend([k2, y2],[L"Filtered~Mean", L"Observation"], loc="best")
  xlim([0, tend])
  # ylim([minimum(xs[2,:]), maximum(xs[2,:])])
  if subplt == 3
    subplot(subplt,1,3)
    plot(ts, us)
    xlim([0, tend])
    ylabel("Heat Input [kJ/min]")
  end
  xlabel(L"Time~[min]")
end

function plotStateSpaceSwitch(linsystems, xs)

  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  figure() # Model and state space
  for k=1:length(linsystems)
    plot(linsystems[k].op[1],linsystems[k].op[2],"kx",markersize=5, markeredgewidth=1)
  annotate(latexstring("M_",k),
        xy=[linsystems[k].op[1],linsystems[k].op[2]],
        xytext=[linsystems[k].op[1],linsystems[k].op[2]],
        fontsize=22.0,
        ha="center",
        va="bottom")
  end
  plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  xlim([-0.1, 1.1])
  xlabel(L"C_A~[kmol.m^{-3}]")
  ylabel(L"T_R~[K]")
end

function plotSwitchSelection(numSwitches, strack, ts, cbaron)

  figure() # Model selection
  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  axes = Array(Any, numSwitches)
  im = 0
  width = 500
  for k=1:numSwitches
    ax = subplot(numSwitches, 1, k)
    axes[k] = ax
    im = imshow(repeat(strack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto")
    tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
    tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
    ylabel(latexstring("M_",k))
  end

  tick_params(axis="x", labelbottom = "on")
  tempts = [1:int(length(ts)/10.0):length(ts)]
  temp = Array(String, length(tempts))
  for lat=1:length(tempts)
    temp[lat] = latexstring(ts[tempts[lat]])
  end

  xticks(tempts, temp)

  if cbaron == true
    colorbar(im, ax=axes)
  end
  xlabel(L"Time~[min]")
end

function plotEllipses(ts, xs, fmeans, fcovars, fname, legloc)

  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  N = length(ts)
  skip = int(length(ts)/20)
  figure()
  b1 = 0.0
  for k=1:N
    p1, p2 = Ellipse.ellipse(fmeans[:,k], fcovars[:,:, k])
    # b1, = plot(p1, p2, "b")
    b1, = fill(p1, p2, "b", edgecolor="none")
  end
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "mx", markersize=5, markeredgewidth = 2)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  ylabel(L"T_R~[K]")
  xlabel(L"C_A~[kmol.m^{-3}]")
  legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", L"90\%~Confidence~Region"], loc=legloc)
end

function plotEllipses(ts, xs, fmeans, fcovars, fname, line, sp, nf, sigma, pick, legloc)

  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  N = length(ts)
  skip = int(length(ts)/30)

  nf && figure() # only create a new figure if required
  b1 = 0.0
  for k=1:N
    p1, p2 = Ellipse.ellipse(fmeans[:,k], fcovars[:,:, k], sigma)
    # b1, = plot(p1, p2, "b")
    b1, = fill(p1, p2, "b", edgecolor="none")
  end
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "mx", markersize=5, markeredgewidth = 2)
  #plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)

  # line = [b,c] => y + bx + c = 0
  # line => y = - bx - c

  lxs = [-0.1:0.05:1.1]
  lys = -line[1].*lxs .- line[2]
  xlim([0.0, 1.0])
  ylim([minimum(xs[2,:]-10), maximum(xs[2, :]+10)])
  plot(lxs, lys, "r-")

  plot(sp[1], sp[2], "gx",markersize=8, markeredgewidth = 4)

  ylabel(L"T_R~[K]")
  xlabel(L"C_A~[kmol.m^{-3}]")
  # conf = round((1.0 - exp(-sigma/2.0))*100.0, 3)
  # temp = latexstring(conf,"\%", "Confidence~Region")
  # legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", temp], loc="best")
  if pick==1
    legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", L"90\%~Confidence~Region"], loc=legloc)
  elseif pick == 2
    legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", L"99\%~Confidence~Region"], loc=legloc)
  elseif pick == 3
    legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", L"99.9\%~Confidence~Region"], loc=legloc)
  else
    legend([x1,f1, b1],[L"Underlying~Model",L"Filtered~Mean", L"99.99\%~Confidence~Region"], loc=legloc)
  end
end

function plotEllipses(fmeans, fcovars, fstart, pmeans, pcovars, pskip::Int64)

  rc("font", family="serif", size=24)
  # rc("text", usetex=true)
  figure() # dont create a new figure
  b1 = 0.0
  b2 = 0.0
  N = size(pmeans)
  for k=1:pskip:N[2]
    p1, p2 = Ellipse.ellipse(pmeans[:,k], pcovars[:,:, k])
    # b1, = plot(p1, p2, "b")
    b1, = fill(p1, p2, "r", edgecolor="none")

    p1, p2 = Ellipse.ellipse(fmeans[:,fstart+k], fcovars[:,:, fstart+k])
    # b1, = plot(p1, p2, "b")
    b2, = fill(p1, p2, "b", edgecolor="none")
  end

  f1, = plot(pmeans[1, 1], pmeans[2, 1], "ro", markersize=5, markeredgewidth = 2)
  f1, = plot(pmeans[1, end], pmeans[2, end], "rx", markersize=5, markeredgewidth = 2)
  f1, = plot(pmeans[1, 2:pskip:end][:], pmeans[2, 2:pskip:end][:], "rx", markersize=5, markeredgewidth = 2)

  # f2, = plot(fmeans[1, 1], fmeans[2, 1], "bo", markersize=5, markeredgewidth = 2)
  # f2, = plot(fmeans[1, end], fmeans[2, end], "bx", markersize=5, markeredgewidth = 2)
  f2, = plot(fmeans[1, fstart+1:pskip:fstart+N[2]][:], fmeans[2, fstart+1:pskip:fstart+N[2]][:], "bx", markersize=5, markeredgewidth = 2)

  ylabel("Temperature [K]")
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  legend([b1, b2],["Prediction", "Filter"], loc="best")
end

function plotEllipseComp(f1means, f1covars, f1name, f2means, f2covars, f2name, xs, ts, line, sp, sigma=4.605)

  N = length(ts)
  skip = 1 # int(length(ts)/20)
  figure()
  # rc("text", usetex=true)
  rc("font", family="serif", size=24)
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  x11, = plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
  f1, = plot(f1means[1, 1:skip:end][:], f1means[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
  f2, = plot(f2means[1, 1:skip:end][:], f2means[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
  b1 = 0.0
  b2 = 0.0
  for k=1:skip:N
    p1, p2 = Ellipse.ellipse(f1means[:,k], f1covars[:,:, k], sigma)
    b1, = plot(p1, p2, "r")

    p3, p4 = Ellipse.ellipse(f2means[:,k], f2covars[:,:, k], sigma)
    b2, = plot(p3, p4, "b")
  end

  # line = [b,c] => y + bx + c = 0
  # line => y = - bx - c

  lxs = [-0.1:0.05:1.1]
  lys = -line[1].*lxs .- line[2]
  plot(lxs, lys, "r-")
  plot(sp[1], sp[2], "gx",markersize=8, markeredgewidth = 4)


  plot(xs[1,:][:], xs[2,:][:], "k", linewidth=3)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  ylabel("Temperature [K]")
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  conf = round((1.0 - exp(-sigma/2.0))*100.0, 3)
  temp1 = string("$(f1name) ", conf,"\% Confidence Ellipse")
  temp2 = string("$(f2name) ", conf,"\% Confidence Ellipse")
  legend([x1,f1,f2, b1, b2],["Nonlinear Model","$(f1name) Mean","$(f2name) Mean", temp1, temp2], loc="best")
end

function plotEllipseComp(f1means, f1covars, f2means, f2covars, xs, ts, sigma=4.605)

  N = length(ts)
  skip = int(length(ts)/30)
  figure()
  rc("text", usetex=true)
  rc("font", family="serif", size=24)
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  f1, = plot(f1means[1, 1:skip:end][:], f1means[2, 1:skip:end][:], "yx", markersize=5, markeredgewidth = 2)
  f2, = plot(f2means[1, 1:skip:end][:], f2means[2, 1:skip:end][:], "gx", markersize=5, markeredgewidth = 2)
  b1 = 0.0
  b2 = 0.0
  for k=1:N
    p1, p2 = Ellipse.ellipse(f1means[:,k], f1covars[:,:, k], sigma)
    # b1, = plot(p1, p2, "r")
    b1, = fill(p1, p2, "r", edgecolor="none")

    p3, p4 = Ellipse.ellipse(f2means[:,k], f2covars[:,:, k], sigma)
    # b2, = plot(p3, p4, "b")
    b2, = fill(p3, p4, "b", edgecolor="none")
  end

  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  ylabel(L"T_R~[K]")
  xlabel(L"C_A~[kmol.m^{-3}]")
  legend([x1,f1,f2, b1, b2],[L"Underlying~Model",L"Particle~Filter",L"Kalman~Filter", L"PF~90\%~Confidence~Region", L"KF~90\%~Confidence~Region"], loc="best")
end

function plotTrackingBreak(ts, xs, xsb, ys, fmeans, obs)

  N = length(ts)
  tend = ts[end]
  skipm = int(length(ts)/20)
  figure() # Plot filtered results
  rc("text", usetex=true)
  rc("font", family="serif", size=24)
  subplot(2,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  x1nf, = plot(ts, xsb[1,:]', "g--", linewidth=3)
  if obs == 2
    y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
  ylabel(L"C_A~[kmol.m^{-3}]")
  legend([x1, k1],[L"Underlying~Model",L"Filtered~Mean"], loc="best")
  xlim([0, tend])
  subplot(2,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  x2nf, = plot(ts, xsb[2,:]', "g--", linewidth=3)
  if obs == 2
    y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  else
    y2, = plot(ts[1:skipm:end], ys[1:skipm:end], "kx", markersize=5, markeredgewidth=1)
  end
  k2, = plot(ts, fmeans[2,:]', "r--", linewidth=3)
  ylabel(L"T_R~[K]")
  xlabel(L"Time~[min]")
  legend([y2, x2nf],[L"Observation",L"Underlying~Model~w/o~Break"], loc="best")
  xlim([0, tend])
end

function plotTrackingTwoFilters(ts, xs, ys, f1means, f2means)

  skipm = int(length(ts)/30)
  skip = int(length(ts)/30)
  tend = ts[end]
  figure() # Plot filtered results
  rc("font", family="serif", size=24)
  rc("text", usetex=true)
  subplot(2,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  k1, = plot(ts[1:skip:end], f1means[1,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
  y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  k12, = plot(ts[1:skip:end], f2means[1, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
  ylabel(L"C_A~[kmol.m^{-3}]")
  legend([x1, k1],[L"Underlying~Model",L"Particle~Filter"], loc="best", ncol=2)
  xlim([0, tend])
  subplot(2,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  k2, = plot(ts[1:skip:end], f1means[2,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
  k22, = plot(ts[1:skip:end], f2means[2, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
  ylabel(L"T_R~[K]")
  xlabel(L"Time~[min]")
  legend([y2, k22],[L"Observation", L"Kalman~Filter"], loc="best", ncol=2)
  xlim([0, tend])
end

function plotTrackingComparison(ts, xs1, us1, xs2, us2, setpoint)

    tend = ts[end]

    setpoints = ones(length(ts))*setpoint

    rc("font", family="serif", size=24)
    # rc("text", usetex=true)

    figure()
    subplot(3,1,1)
    x11, = plot(ts, xs1[1,:]', "r", linewidth=3)
    x12, = plot(ts, xs2[1,:]', "b", linewidth=3)
    ksp = plot(ts, setpoints, "g-", linewidth=3)
    ylabel(L"Concentration [kmol.m$^{-3}$]")
    legend([x11],["Switching Controller"], loc="best")
    xlim([0, tend])
    ylim([0, 1])

    subplot(3,1,2)
    x12, = plot(ts, xs1[2,:]', "r", linewidth=3)
    x22, = plot(ts, xs2[2,:]', "b", linewidth=3)
    ylabel("Temperature [K]")
    legend([x22],["Static Controller"], loc="best")
    xlim([0, tend])
    ylim([minimum(xs2[2,:]), maximum(xs2[2,:])])

    subplot(3,1,3)
    u1, = plot(ts, us1, "r", linewidth=3)
    u2, = plot(ts, us2, "b", linewidth=3)
    xlim([0, tend])
    ylabel("Controller Input")
    xlabel("Time [min]")
end

function plotKLdiv(ts, kldiv, basediv, unidiv, logged)
  rc("text", usetex=true)
  rc("font", family="serif", size=24)

  figure()
  kl, = plot(ts, kldiv, "r", linewidth=3)
  gd, = plot(ts, basediv, "b", linewidth=3)
  ud, = plot(ts, unidiv, "g", linewidth=3)
  xlabel(L"Time~[min]")
  if logged
    ylabel(L"Divergence~[ln(Nats)]")
  else
    ylabel(L"Divergence~[Nats]")
  end
  legend([kl, gd, ud],[L"Approximation",L"Baseline", L"Uniform"], loc="upper right")
end

function calcError(x, y::Array{Float64, 2})

  r, N = size(x)
  avediff1 = (1.0/N)*sum(abs((x[1, :].-y[1, :])./x[1,:]))*100.0
  avediff2 = (1.0/N)*sum(abs((x[2, :].-y[2, :])./x[2,:]))*100.0

  println("Average Concentration Error: ", round(avediff1, 4),  "%")
  println("Average Temperature Error: ", round(avediff2, 4), "%")
  return avediff1, avediff2
end

function calcError(x, y::Float64)

  r, N = size(x)
  avediff1 = (1.0/N)*sum(abs((x[1, :] - y)./y))*100.0

  println("Average Concentration Error: ", round(avediff1, 4),  "%")
  return avediff1
end

function calcEnergy(us, uss, h)
  N = length(us)
  avecost = (h/N)*sum(abs(us-uss))
  println("Average Input: ", avecost)
  return avecost
end

function checkConstraint(ts, xs, line)
  # line = [b,c] => y + bx + c = 0
  # line => y = - bx - c
  r, N = size(xs)
  conmargin = zeros(N)
  minneg = Inf
  minpos = Inf
  for k=1:N
    temp = xs[2, k] + xs[1, k]*line[1] + line[2]
    if temp < 0.0
      conmargin[k] = -abs(temp)/sqrt(line[1]^2 + 1.0)
      if minneg > abs(temp)/sqrt(line[1]^2 + 1.0)
        minneg = abs(temp)/sqrt(line[1]^2 + 1.0)
      end
    else
      conmargin[k] = abs(temp)/sqrt(line[1]^2 + 1.0)
      if minpos > abs(temp)/sqrt(line[1]^2 + 1.0)
        minpos = abs(temp)/sqrt(line[1]^2 + 1.0)
      end
    end
  end


  println("Minimum Positive Clearance: ", minpos)
  println("Minimum Negative Clearance: ", minneg)

  figure()
  rc("font", family="serif", size=24)
  rc("text", usetex=true)

  plot(ts, zeros(N), "r", linewidth=1)
  plot(ts, conmargin, "k", linewidth=3)
  xlabel(L"Time [min]")
  ylabel(L"Clearance")
end

end #module
