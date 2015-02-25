module LLDS_functions

using Distributions

immutable LLDS{T} # so that univariate and multivariate systems can be handled easily
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

function step(xprev::Float64, yprev::Float64, model::LLDS)

  dprocess = Normal(model.Q)
  dmeasure = Normal(model.R)

  return xnow, ynow
end



end #module
