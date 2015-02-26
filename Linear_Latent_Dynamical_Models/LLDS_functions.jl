module LLDS_functions

using Distributions

immutable LLDS
  # Assume zero mean transition and emission functions.
  # The linear latent dynamical system should have the
  # state space form:
  # x(t+1) = A*x(t) + Bu(t) + b + Q
  # y(t+1) = C*x(t+1) + Du(t) + R
  A :: Array{Float64, 2}
  B :: Array{Float64, 2}
  b :: Array{Float64, 1}
  C :: Array{Float64, 2}
  D :: Array{Float64, 2}
  Q :: Array{Float64, 2} # Process Noise
  R :: Array{Float64, 2} # Measurement Noise
end

function step(xprev_noisy::Array{Float64,1}, unow::Array{Float64,1}, model::LLDS)
  # Controlled, move multivariate model one time step forward.
  dprocess = MvNormal(model.Q)
  dmeasure = MvNormal(model.R)

  xnow_noisy = model.A*xprev_noisy + model.B*unow + model.b + rand(dprocess)
  ynow_noisy = model.C*xnow_noisy + model.D*unow + rand(dmeasure)
  ynow = model.C*xnow_noisy +model.D*unow

  return xnow_noisy, ynow_noisy, ynow
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
  pmean = model.A*xprev + model.B*unow + model.b
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

function smooth(kmeans::Array{Float64, 2}, kcovars::Array{Float64, 3}, us::Array{Float64, 2}, model::LLDS)
  # Returns the smoothed means and covariances
  rows, cols = size(kmeans)
  smoothedmeans = zeros(rows, cols)
  smoothedvars = zeros(rows, rows, cols)
  smoothedmeans[:, end] = kmeans[:, end]
  smoothedvars[:, :, end] = kcovars[:, :, end]

  for t=cols-1:-1:1
    Pt = model.A*kcovars[:, :, t]*transpose(model.A) + model.Q
    Jt = kcovars[:, :, t]*transpose(model.A)*inv(Pt)
    smoothedmeans[:, t] = kmeans[:, t] + Jt*(smoothedmeans[:, t+1] - model.A*kmeans[:, t] - model.B*us[:, t])
    smoothedvars[:,:, t] = kcovars[:,:, t] + Jt*(smoothedvars[:,:, t+1] - Pt)*transpose(Jt)
  end

  return smoothedmeans, smoothedvars
end

function predict_visible(kmean::Array{Float64, 1}, kcovar::Array{Float64, 2}, us::Array{Float64, 2}, model::LLDS)
  # Predict the visible states n steps into the future given the controller action.
  rows, = size(kmean)
  rus, n = size(us)
  predicted_means = zeros(rows, n)
  predicted_covars = zeros(rows, rows, n)

  predicted_means[:, 1] = model.A*kmean + model.B*us[:,1] + model.b
  predicted_covars[:, :, 1] = model.Q + model.A*kcovar*transpose(model.A)

  for k=2:n #cast the state forward
    predicted_means[:, k], predicted_covars[:, :, k] = step_predict(predicted_means[:,k-1], predicted_covars[:, :, k-1],us[:,k], model)
  end

  rows, cols = size(model.R)
  predicted_vis_means = zeros(rows, n)
  predicted_vis_covars = zeros(rows, cols, n)

  for k=1:n # convert the hidden state to the observed state
    predicted_vis_means[:, k] = model.C*predicted_means[:,k] + model.D*us[:,k]

    predicted_vis_covars[:, :, k] = model.R + model.C*predicted_covars[:, :, k]*transpose(model.C)
  end

  return predicted_vis_means, predicted_vis_covars
end

function predict_hidden(kmean::Array{Float64, 1}, kcovar::Array{Float64, 2}, us::Array{Float64, 2}, model::LLDS)
  # Predict the hidden states n steps into the future given the controller action.
  rows, = size(kmean)
  rus, n = size(us)
  predicted_means = zeros(rows, n)
  predicted_covars = zeros(rows, rows, n)

  predicted_means[:, 1] = model.A*kmean + model.B*us[:,1] + model.b
  predicted_covars[:, :, 1] = model.Q + model.A*kcovar*transpose(model.A)

  for k=2:n #cast the state forward
    predicted_means[:, k], predicted_covars[:, :, k] = step_predict(predicted_means[:,k-1], predicted_covars[:, :, k-1],us[:,k], model)
  end

  return predicted_means, predicted_covars
end

end #module
