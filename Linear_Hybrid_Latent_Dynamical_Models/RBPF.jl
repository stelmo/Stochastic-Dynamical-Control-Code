# Rao Blackwellised Particle Filter
module RBPF
# WARNING: this is made specifically for the system I am investigating

using Distributions
cd("..\\CSTR_Model")
using Reactor_functions
cd("..\\Linear_Latent_Dynamical_Models")
using LLDS_functions
cd("..\\Linear_Hybrid_Latent_Dynamical_Models")
using SPF

type Particles
  mus :: Array{Float64, 1} # mean
  sigmas :: Array{Float64, 2} # covariance
  ss :: Array{Int64, 1} # switches
  ws :: Array{Float64, 1} # weights
end

type Model
  A
  B
  b
  C
  Q
  R
end

function setup_RBPF(linsystems, C, Q, R)
  # Setup each switch
  N = length(linsystems)
  models = Array(Model, N)
  for k=1:N
    models[k] = Model(linsystems[k].A, linsystems[k].B, linsystems[k].b, C, Q, R)
  end
  A = SPF.calcA(linsystems)
  return models, A
end

function init_RBPF()

end


end #module
