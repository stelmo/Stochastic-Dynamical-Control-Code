module MPC

using LQR # for the more robust getSize
# Auxiliary functions for the MPC

function getQuad(Q, R, nQ, nR)
  # Constructs the quadratic term of a quadratic objective function
  sQ = LQR.getSize(Q)
  sR = LQR.getSize(R)

  trows = nQ*sQ[1] + nR*sR[1]
  tcols = nQ*sQ[2] + nR*sR[2]
  quadmat = zeros(trows, tcols) # create the big matrix

  # we know the matrix must be square and symmetric
  diagonal = 1
  counter = 1
    for k=1:nQ+nR
      if counter <= nQ
        quadmat[diagonal:diagonal+sQ[1]-1, diagonal:diagonal+sQ[1]-1] = Q
        diagonal += sQ[2]
      else
        quadmat[diagonal:diagonal+sR[1]-1, diagonal:diagonal+sR[2]-1] = R
        diagonal += sR[2]
      end
      counter += 1
  end

  return quadmat
end


end
