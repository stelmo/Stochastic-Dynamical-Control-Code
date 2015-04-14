#  LQR Controller
module LQR

function lqr(A, B, Q, R)
  # Returns the infinite horizon LQR solution.
  # Don't confuse the weighting matrices Q and R with
  # the noise covariance matrices!
  P = dare(A, B, Q, R)
  F = inv(R+B'*P*B)*B'*P*A
  return F
end

function dare(A, B, Q, R)
  # Solves the discrete algebraic ricatti equation (dare)

  nr,nc = size(A)
  P = eye(nr)
  Pnow = eye(nr)
  counter = 0
  while true
    Pnow = Q + A'*P*A - A'*P*B*inv(B'*P*B+R)*B'*P*A
    if norm(Pnow-P, Inf) < 1E-06
      P = copy(Pnow)
      return P
    else
      P = copy(Pnow)
    end
    (counter > 1000000) && (error("DARE did not converge..."))
    counter += 1
  end
end

function inv(x::Array{Float64, 1})
  if length(x) == 1
    return 1./x[1]
  else
    error("Cannot invert a vector with more than one entry!")
  end
end

function offset(A, B, C, H, ysp)
  # Returns the state and controller offset.
  rA, cA = getSize(A)
  rB, cB = getSize(B)
  rC, cC = getSize(C)
  rH, cH = getSize(H)
  lenysp = length(ysp)
  if cA+cB-cC == 0
    z1 = zeros(rA+rH-rB)
  else
    z1 = zeros(rA+rH-rB,cA+cB-cC)
  end
  z2 = zeros(rA+rH-lenysp)
  ssvec = [z2; ysp]
  ssmat = [eye(rA)-A -B;H*C z1]
  ss = ssmat\ssvec

  x_off = ss[1:rA]
  u_off = ss[rA+1:end]
  return x_off, u_off
end

function getSize(A)
  # wraps size() but with added robustness
  x = size(A)
  r = 0
  c = 0
  if length(x) == 1
    r = x[1]
    c = 0
  else
    r = x[1]
    c = x[2]
  end

  return r, c
end


end #module
