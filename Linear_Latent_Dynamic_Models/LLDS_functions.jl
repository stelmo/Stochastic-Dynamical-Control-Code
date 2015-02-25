module LLDS_functions

using Distributions

immutable LLDS{T}
  # T = Float64 or Array{Float64, 2} so that univariate
  # and multivariate systems can be handled easily.
  # NB: all matrices need to have two dimensions.
  # Assume zero mean transition and emission functions.
  # The linear latent dynamical system should have the
  # state space form:
  # x(t+1) = A*x(t) + Bu(t+1) + Q
  # y(t+1) = C*x(t+1) + Du(t+1) + R
  A :: T
  B :: T
  C :: T
  D :: T
  Q :: T # Process Noise
  R :: T # Measurement Noise
end

function step(xprev::Float64, model::LLDS{Float64})
  # No control, move univariate model one time step forward.
  dprocess = Normal(model.Q)
  dmeasure = Normal(model.R)
  xnow = model.A*xprev + rand(dprocess)
  ynow = model.C*xnow + rand(dmeasure)

  return xnow, ynow
end

function step(xprev::Float64, unow::Float64, model::LLDS{Float64})
  # Controlled, move univariate model one time step forward.
  dprocess = Normal(model.Q)
  dmeasure = Normal(model.R)
  xnow = model.A*xprev + model.B*unow + rand(dprocess)
  ynow = model.C*xnow +model.D*unow + rand(dmeasure)

  return xnow, ynow
end

function step(xprev::Array{Float64,1}, unow::Array{Float64,1}, model::LLDS{Array{Float64,2}})
  # Controlled, move multivariate model one time step forward.
  dprocess = MvNormal(model.Q)
  dmeasure = MvNormal(model.R)
  xnow = model.A*xprev + model.B*unow + rand(dprocess)
  ynow = model.C*xnow +model.D*unow + rand(dmeasure)

  return xnow, ynow
end

function step(xprev::Array{Float64,1}, model::LLDS{Array{Float64,2}})
  # No control, move multivariate model one time step forward.
  dprocess = MvNormal(model.Q)
  dmeasure = MvNormal(model.R)
  xnow = model.A*xprev + rand(dprocess)
  ynow = model.C*xnow + rand(dmeasure)

  return xnow, ynow
end

function filter()
  # Step the filter by one time step
  
end

end #module
