# Plotting and results analysis module
module Results

using PyPlot
using Ellipse

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
  ylabel(L"Concentration [kmol.m$^{-3}$]")
  legend([x1],["Nonlinear Model"], loc="best")
  xlim([0, tend])
  ylim([0, 1])

  subplot(subplt,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  if obs == 1
    y2, = plot(ts[1:skipmeas:end], ys[1:skipmeas:end], "kx", markersize=5, markeredgewidth=1)
  else
    y2, = plot(ts[1:skipmeas:end], ys[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k2, = plot(ts[1:skipmean:end], fmeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel("Temperature [K]")
  legend([k2, y2],["Filtered Mean Estimate", "Nonlinear Model Measured"], loc="best")
  xlim([0, tend])
  ylim([minimum(xs[2,:]), maximum(xs[2,:])])
  if subplt == 3
    subplot(subplt,1,3)
    plot(ts, us)
    xlim([0, tend])
    ylabel("Controller Input")
  end
  xlabel("Time [min]")
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

  skipmeas = int(length(ts)/20)
  skipmean = int(length(ts)/20)
  figure()
  subplot(subplt,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  if obs == 2 # plot second measurement
    y2, = plot(ts[1:skipmeas:end], ys[1, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k1, = plot(ts[1:skipmean:end], fmeans[1, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel(L"Concentration [kmol.m$^{-3}$]")
  legend([x1],["Nonlinear Model"], loc="best")
  xlim([0, tend])
  ylim([0, 1])

  subplot(subplt,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  if obs == 1
    y2, = plot(ts[1:skipmeas:end], ys[1:skipmeas:end], "kx", markersize=5, markeredgewidth=1)
  else
    y2, = plot(ts[1:skipmeas:end], ys[2, 1:skipmeas:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k2, = plot(ts[1:skipmean:end], fmeans[2, 1:skipmean:end]', "bx", markersize=5, markeredgewidth = 2)
  ylabel("Temperature [K]")
  legend([k2, y2],["Filtered Mean Estimate", "Nonlinear Model Measured"], loc="best")
  xlim([0, tend])
  ylim([minimum(xs[2,:]), maximum(xs[2,:])])
  if subplt == 3
    subplot(subplt,1,3)
    plot(ts, us)
    xlim([0, tend])
    ylabel("Controller Input")
  end
  xlabel("Time [min]")
end

function plotStateSpaceSwitch(linsystems, xs)

  rc("font", family="serif", size=24)

  figure() # Model and state space
  for k=1:length(linsystems)
    plot(linsystems[k].op[1],linsystems[k].op[2],"kx",markersize=5, markeredgewidth=1)
  annotate(string("Switch: ", k),
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
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  ylabel("Temperature [K]")
end

function plotSwitchSelection(numSwitches, strack, ts, cbaron)
  figure() # Model selection
  axes = Array(Any, numSwitches)
  im = 0
  width = 500
  for k=1:numSwitches
    ax = subplot(numSwitches, 1, k)
    axes[k] = ax
    im = imshow(repeat(strack[k,:], outer=[width, 1]), cmap="cubehelix",vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto")
    tick_params(axis="y", which="both",left="off",right="off", labelleft = "off")
    tick_params(axis="x", which="both",bottom="off", labelbottom = "off")
    ylabel(string("S::",k))
  end

  tick_params(axis="x", labelbottom = "on")
  xticks([1:int(length(ts)/10.0):length(ts)], ts[1:int(length(ts)/10.0):end])

  if cbaron == true
    colorbar(im, ax=axes)
  end
  xlabel("Time [min]")
end

function plotEllipses(ts, xs, fmeans, fcovars, fname)

  rc("font", family="serif", size=24)
  N = length(ts)
  skip = int(length(ts)/20)
  figure()
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
  b1 = 0.0
  for k=1:skip:N
    p1, p2 = Ellipse.ellipse(fmeans[:,k], fcovars[:,:, k])
    b1, = plot(p1, p2, "b")
  end
  plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
  plot(xs[1,1], xs[2,1], "ko", markersize=10, markeredgewidth = 4)
  plot(xs[1,end], xs[2,end], "kx", markersize=10, markeredgewidth = 4)
  ylabel("Temperature [K]")
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  temp = string("$(fname) ", L"$\sigma$-Ellipse")
  legend([x1,f1, b1],["Nonlinear Model","$(fname) Mean", temp], loc="best")
end

function plotEllipses(ts, xs, fmeans, fcovars, fname, line, sp, nf, sigma=4.605)

  rc("font", family="serif", size=24)
  N = length(ts)
  skip = 1

  nf && figure() # only create a new figure if required

  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  f1, = plot(fmeans[1, 1:skip:end][:], fmeans[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
  b1 = 0.0
  for k=1:skip:N
    p1, p2 = Ellipse.ellipse(fmeans[:,k], fcovars[:,:, k], sigma)
    b1, = plot(p1, p2, "b")
  end
  plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
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

  ylabel("Temperature [K]")
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  temp = string("$(fname) ", L"$\sigma$-Ellipse")
  legend([x1,f1, b1],["Nonlinear Model","$(fname) Mean", temp], loc="best")
end

function plotEllipses(fmeans, fcovars, fstart, pmeans, pcovars, pskip::Int64)

  rc("font", family="serif", size=24)
  figure() # dont create a new figure
  N = size(pmeans)
  f1, = plot(pmeans[1, 1], pmeans[2, 1], "ro", markersize=5, markeredgewidth = 2)
  f1, = plot(pmeans[1, end], pmeans[2, end], "rx", markersize=5, markeredgewidth = 2)
  f1, = plot(pmeans[1, 2:pskip:end][:], pmeans[2, 2:pskip:end][:], "rx", markersize=5, markeredgewidth = 2)

  # f2, = plot(fmeans[1, 1], fmeans[2, 1], "bo", markersize=5, markeredgewidth = 2)
  # f2, = plot(fmeans[1, end], fmeans[2, end], "bx", markersize=5, markeredgewidth = 2)
  f2, = plot(fmeans[1, fstart+1:pskip:fstart+N[2]][:], fmeans[2, fstart+1:pskip:fstart+N[2]][:], "bx", markersize=5, markeredgewidth = 2)

  b1 = 0.0
  b2 = 0.0
  for k=1:pskip:N[2]
    p1, p2 = Ellipse.ellipse(pmeans[:,k], pcovars[:,:, k])
    b1, = plot(p1, p2, "r")

    p1, p2 = Ellipse.ellipse(fmeans[:,fstart+k], fcovars[:,:, fstart+k])
    b2, = plot(p1, p2, "b")
  end
  ylabel("Temperature [K]")
  xlabel(L"Concentration [kmol.m$^{-3}$]")
  legend([b1, b2],["Prediction", "Filter"], loc="best")
end

function plotEllipseComp(f1means, f1covars, f1name, f2means, f2covars, f2name, xs, ts, line, sp)

  N = length(ts)
  skip = int(length(ts)/20)
  figure()
  x1, = plot(xs[1,:][:], xs[2,:][:], "k",linewidth=3)
  x11, = plot(xs[1, 1:skip:end][:], xs[2, 1:skip:end][:], "kx", markersize=5, markeredgewidth = 2)
  f1, = plot(f1means[1, 1:skip:end][:], f1means[2, 1:skip:end][:], "rx", markersize=5, markeredgewidth = 2)
  f2, = plot(f2means[1, 1:skip:end][:], f2means[2, 1:skip:end][:], "bx", markersize=5, markeredgewidth = 2)
  b1 = 0.0
  b2 = 0.0
  for k=1:skip:N
    p1, p2 = Ellipse.ellipse(f1means[:,k], f1covars[:,:, k])
    b1, = plot(p1, p2, "r")

    p3, p4 = Ellipse.ellipse(f2means[:,k], f2covars[:,:, k])
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
  temp1 = string("$(f1name) ", L"$\sigma$-Ellipse")
  temp2 = string("$(f2name) ", L"$\sigma$-Ellipse")
  legend([x1,f1,f2, b1, b2],["Nonlinear Model","$(f1name) Mean","$(f2name) Mean", temp1, temp2], loc="best")
end

function plotTrackingBreak(ts, xs, xsb, ys, fmeans, obs)

  N = length(ts)
  tend = ts[end]
  skipm = int(length(ts)/20)
  figure() # Plot filtered results
  subplot(2,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  x1nf, = plot(ts, xsb[1,:]', "g--", linewidth=3)
  if obs == 2
    y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  end
  k1, = plot(ts, fmeans[1,:]', "r--", linewidth=3)
  ylabel(L"Concentration [kmol.m$^{-3}$]")
  legend([x1, k1],["Nonlinear Model","Filtered Mean"], loc="best")
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
  ylabel("Temperature [K]")
  xlabel("Time [min]")
  legend([y2, x2nf],["Nonlinear Model Measured","Nonlinear Model No Switch"], loc="best")
  xlim([0, tend])
end

function plotTrackingTwoFilters(ts, xs, ys, f1means, f2means, f1name, f2name)

  skipm = int(length(ts)/20)
  skip = int(length(ts)/20)
  tend = ts[end]
  figure() # Plot filtered results
  subplot(2,1,1)
  x1, = plot(ts, xs[1,:]', "k", linewidth=3)
  k1, = plot(ts[1:skip:end], f1means[1,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
  y2, = plot(ts[1:skipm:end], ys[1, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  k12, = plot(ts[1:skip:end], f2means[1, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
  ylabel(L"Concentration [kmol.m$^{-3}$]")
  legend([x1, k1],["Nonlinear Model","$(f1name)"], loc="best")
  xlim([0, tend])
  subplot(2,1,2)
  x2, = plot(ts, xs[2,:]', "k", linewidth=3)
  y2, = plot(ts[1:skipm:end], ys[2, 1:skipm:end][:], "kx", markersize=5, markeredgewidth=1)
  k2, = plot(ts[1:skip:end], f1means[2,1:skip:end]', "rx", markersize=5, markeredgewidth=2)
  k22, = plot(ts[1:skip:end], f2means[2, 1:skip:end]', "bx", markersize=5, markeredgewidth=2)
  ylabel("Temperature [K]")
  xlabel("Time [min]")
  legend([y2, k22],["Nonlinear Model Measured", "$(f2name)"], loc="best")
  xlim([0, tend])
end

function plotTrackingComparison(ts, xs1, us1, xs2, us2, setpoint)

    tend = ts[end]

    setpoints = ones(length(ts))*setpoint

    rc("font", family="serif", size=24)

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

end #module
