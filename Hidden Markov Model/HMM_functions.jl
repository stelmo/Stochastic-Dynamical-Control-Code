# Functions which implement inference algorithms for HMMs.
# Discrete, time invariant transition and emission probabilities are assumed.
module HMM_functions

# Add export functions
export  HMM, forward, normalise, baum_welch, viterbi, prediction

immutable HMM #time invariant
  tp :: Array{Float64, 2} #transmission probability table
  ep :: Array{Float64, 2} #emission probability table

  HMM(tp, ep) = new(tp, ep)
end

function forward(model::HMM, initial::Array{Float64, 1}, evidence::Array{Int64,1})
  # Forwards algorithm for HMM
  ns, = size(model.tp) #number of states
  ne, = size(evidence) #number of observations
  alpha = zeros(ns, ne) #matrix of forward probabilities

  for ks=1:ns # first observation
    alpha[ks, 1] = model.ep[evidence[1], ks]*initial[ks]
  end
  alpha[:,1] = normalise(alpha[:,1]) # normalise probabilities

  for ke=2:ne # loop through each evidence observation

    for ks=1:ns # loop through each state to build up the joint given observations

      prediction = 0.0 # the prediction part given the old joint
      for ks_pred=1:ns
        prediction += model.tp[ks, ks_pred]*alpha[ks_pred, ke-1] # can be made faster by looping column major!
      end

      alpha[ks, ke] = model.ep[evidence[ke], ks]*prediction
    end

    alpha[:, ke] = normalise(alpha[:,ke]) # normalise probabilities
  end

  return alpha
end

function normalise(vec::Array{Float64,1})
  # Normalise mat such that sum(vec) = 1
  return vec./sum(vec)
end

function backward(model::HMM, evidence::Array{Int64,1})
  # Backwards algorithm for HMM.
  ns, = size(model.tp) #number of states
  ne, = size(evidence) #number of observations
  beta = zeros(ns, ne) #matrix of forward probabilities

  #initialise
  beta[:,end] = 1.0 # for correct Bayes

  for ke=ne:-1:2 # iterate backwards over evidence, evidence at t does not matter

    for ks=1:ns # iterate over states

      recur = 0.0
      for ks_next=1:ns #sum over next state
        recur += model.ep[evidence[ke], ks_next]*model.tp[ks_next, ks]*beta[ks_next, ke]
      end
      beta[ks, ke-1] = recur

    end
    beta[:, ke-1] = normalise(beta[:, ke-1])
  end
  return beta
end

function baum_welch(model::HMM, initial::Array{Float64,1}, evidence::Array{Int64,1}, timeLocation::Int64)
  # Forwards-Backwards algorithm. Note that it is required to split the evidence
  # accordingly.
  forwardEvidence = evidence[1:timeLocation]
  backwardEvidence = evidence[timeLocation:end]

  alpha = forward(model, initial, forwardEvidence)[:, end]
  beta = backward(model, backwardEvidence)[:,1]

  smoothed = normalise(alpha.*beta)

  return smoothed
end

function viterbi(model::HMM, initial::Array{Float64, 1}, evidence::Array{Int64, 1})
  # The Viterbi algorithm for the most likely joint sequence of states given observations.
  ns, = size(model.tp) #number of states
  ne, = size(evidence) #number of observations
  mu = zeros(ns, ne) #matrix of most likely values per state

  ## Find mu
  #initialise recursion
  mu[:, end] = 1.0

  for ke=ne:-1:2 # iterate backwards over evidence

    for ks=1:ns #iterate over each state
      mu[ks,ke-1] = max_viterbi(model, ks, evidence[ke], mu[:, ke])
    end
    mu[:, ke-1] = normalise(mu[:, ke-1])

  end

  ## Find sequence of states using mu
  mlss = zeros(Int64, ne)
  mlss[1] = arg_viterbi(model, initial, evidence[1], mu[:,1])

  for ke=2:ne # find the most likely sequence of states
    mlss[ke] = arg_viterbi(model, evidence[ke], mu[:,ke], mlss[ke-1])
  end

  return mlss
end

function max_viterbi(model::HMM, state::Int64, evidence::Int64, mu_before::Array{Float64,1})
  # Finds the maximum mu factor given the evidence and previous state
  ns, = size(model.tp) #number of states
  vmax = 0.0

  for k=1:ns
    tmax = model.tp[k, state]*model.ep[evidence, k]*mu_before[k] :: Float64
    tmax > vmax && (vmax = tmax)
  end

  return vmax
end

function arg_viterbi(model::HMM, initial::Array{Float64,1}, evidence::Int64, mu_vec::Array{Float64,1})
  # Finds the most likely state given the first observation

  ns, = size(model.tp) #number of states
  mls = 0 #most likely state
  smax = 0.0 #state value maximum

  for k=1:ns
    tmax = initial[k]*model.ep[evidence, k]*mu_vec[k]

    if tmax > smax
       mls = k
       smax = tmax
    end
  end

  return mls
end

function arg_viterbi(model::HMM, evidence::Int64, mu_vec::Array{Float64,1}, prev_mls::Int64)
  # Finds the most likely sequence of states using mu
  ns, = size(model.tp) #number of states
  mls = 0 #most likely state
  smax = 0.0 #state value maximum

  for ks=1:ns
    tmax = model.tp[ks, prev_mls]*model.ep[evidence, ks]*mu_vec[ks]

    if tmax > smax
       mls = ks
       smax = tmax
    end
  end

  return mls
end

function prediction(model::HMM, initial::Array{Float64, 1}, evidence::Array{Int64, 1})
  # One step ahead hidden state and observation estimator
  ns, = size(model.tp)
  ne, = size(model.ep)

  ## State prediction
  pred_states = zeros(ns)
  filter = forward(model, initial, evidence)[:, end]
  for ks=1:ns
    temp = 0.0
    for k=1:ns
      temp += filter[k]*model.tp[ks, k]
    end
    pred_states[ks] = temp
  end

  ## Evidence Prediction
  pred_evidence = zeros(ne)
  for ke=1:ne
    temp = 0.0
    for ks1=1:ns #h_t
      for ks2=1:ns #h_t+1
        temp += filter[ks1]*model.tp[ks2, ks1]*model.ep[ke, ks2]
      end
    end
  pred_evidence[ke] = temp
  end

  return pred_states, pred_evidence
end

end # Module
