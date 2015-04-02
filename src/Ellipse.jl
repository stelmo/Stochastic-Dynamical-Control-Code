# Plot 1-Sigma Ellipses
module Ellipse

function ellipse(mean, covariance)
  # Plots a confidence ellipse around 2D Gaussians.
  # Converted to Julia from the Matlab function:
  # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

  # Calculate the eigenvectors and eigenvalues
  eigenval, eigenvec = eig(covariance);

  if eigenval[1] > eigenval[2]
    largest_eigenval = eigenval[1]
    largest_eigenvec = eigenvec[:, 1]
    smallest_eigenval = eigenval[2]
    smallest_eigenvec = eigenvec[:,2]
  else
    largest_eigenval = eigenval[2]
    largest_eigenvec = eigenvec[:, 2]
    smallest_eigenval = eigenval[1]
    smallest_eigenvec = eigenvec[:,1]
  end

  angle = atan2(largest_eigenvec[2], largest_eigenvec[1])

  if angle < 0.0
      angle = angle + 2*pi
  end

  chisquare_val = sqrt(2.2788) # one sigma ellipse 68% confidence
  theta_grid = linspace(0.0,2.0*pi)
  phi = angle
  X0=mean[1]
  Y0=mean[2]
  a=chisquare_val*sqrt(largest_eigenval)
  b=chisquare_val*sqrt(smallest_eigenval)

  # the ellipse in x and y coordinates
  ellipse_x_r  = a.*cos( theta_grid )
  ellipse_y_r  = b.*sin( theta_grid )

  #Define a rotation matrix
  R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ]

  # rotate the ellipse to some angle phi
  r_ellipse = zeros(2, length(ellipse_x_r))
  for k=1:length(ellipse_x_r)
    r_ellipse[:, k] = R * [ellipse_x_r[k]; ellipse_y_r[k]] + mean
  end
  return r_ellipse[1,:][:], r_ellipse[2,:][:]
end

end
