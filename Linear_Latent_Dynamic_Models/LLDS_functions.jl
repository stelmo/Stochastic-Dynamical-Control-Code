module LLDS_functions

using Distributions

immutable LLDS
  # Assume zero mean transition and emission functions.
  # The linear latent dynamical system should have the
  # state space form:
  # x(t+1) = A*x(t) + Bu(t+1) + Q
  # y(t+1) = C*x(t+1) + Du(t+1) + R
  A :: Array{Float64, 2}
  B :: Array{Float64, 2}
  C :: Array{Float64, 2}
  D :: Array{Float64, 2}
  Q :: Array{Float64, 2} # Process Noise
  R :: Array{Float64, 2} # Measurement Noise
end

function step(xprev::Array{Float64,1}, unow::Array{Float64,1}, model::LLDS)
  # Controlled, move multivariate model one time step forward.
  dprocess = MvNormal(model.Q)
  dmeasure = MvNormal(model.R)
  xnow = model.A*xprev + model.B*unow + rand(dprocess)
  ynow = model.C*xnow +model.D*unow + rand(dmeasure)

  return xnow, ynow
end

function init_filter(initmean::Array{Float64, 1}, initvar::Array{Float64, 2}, unow::Array{Float64, 1}, ynow::Array{Float64, 1}, model::LLDS)
  # Initialise the filter. No prediction step, only a measurement update step.
  updatedMean ::Array{Float64, 1}, updatedVar :: Array{Float64, 2} = step_update(initmean, initvar, unow, ynow, model)
  return updatedMean, updatedVar
end

function step_filter(prevmean::Array{Float64, 1}, prevvar::Array{Float64, 2}, unow::Array{Float64,1},ynow::Array{Float64, 1}, model::LLDS)
  # Return the posterior over the current state given the observation and previous
  # filter result.
  pmean :: Array{Float64, 1}, pvar :: Array{Float64, 2} = step_predict(prevmean, prevvar, unow, model)
  updatedMean ::Array{Float64, 1}, updatedVar :: Array{Float64, 2} = step_update(pmean, pvar, unow, ynow, model)
  return updatedMean, updatedVar
end

function step_predict(xprev::Array{Float64,1}, varprev::Array{Float64, 2}, unow::Array{Float64,1}, model::LLDS)
  # Return the one step ahead predicted mean and covariance.
  pmean = model.A*xprev + model.B*unow
  pvar =  model.Q + model.A*varprev*transpose(model.A)
  return pmean, pvar
end

function step_update(pmean::Array{Float64,1}, pvar::Array{Float64, 2}, unow::Array{Float64,1}, ymeas::Array{Float64,1}, model::LLDS)
  # Return the one step ahead measurement updated mean and covar.
  kalmanGain = pvar*transpose(model.C)*inv(model.C*pvar*transpose(model.C) + model.R)
  ypred = model.C*pmean + model.D*unow #predicted measurement
  updatedMean = pmean + kalmanGain*(ymeas - ypred)
  rows, cols = size(pvar)
  updatedVar = (eye(rows) - kalmanGain*model.C)*pvar
  return updatedMean, updatedVar
end

function smooth(kmeans::Array{Float64, 2}, kcovars::Array{Float64, 2}, us::Array{Float64, 2}, model::LLDS)
  # Returns the smoothed means and covariances
  rows, cols = size(kmeans)
  smoothedmeans = zeros(rows, cols)
  smoothedvars = zeros(rows, rows*cols)
  smoothedmeans[:, end] = kmeans[:, end]
  smoothedvars[:, end-5:end] = kcovars[:, end-5:end]

  for t=cols-1:-1:1
    Pt = model.A*kcovars[:, (1+(t-1)*6):t*6]*transpose(model.A) + model.Q
    Jt = kcovars[:, (1+(t-1)*6):t*6]*transpose(model.A)*inv(Pt)
    smoothedmeans[:, t] = kmeans[:, t] + Jt*(smoothedmeans[:, t+1] - model.A*kmeans[:, t] - model.B*us[:, t])
    smoothedvars[:, (1+(t-1)*6):t*6] = kcovars[:, (1+(t-1)*6):t*6] + Jt*(smoothedvars[:, (1+(t)*6):(t+1)*6] - Pt)*transpose(Jt)
  end

  return smoothedmeans, smoothedvars
end


end #module
