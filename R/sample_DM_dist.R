
# sample from a DM distributions
#output is a list with both x- and y-samples

sampleDM_dist = function(assessments, weights_, DMx, L, U, quantiles, numSamplePoints) {

  Nexperts = length(assessments)      # number of experts
  N = length(assessments[[1]][,1]) # number of calibration questions
  #numSamplePoints = 50000

  # get x coordinates if needed
  if (TRUE == is.null(DMx)) {
    DMx = list()
    for (q in 1:N) {
      x = list()
      for (e in 1:Nexperts) {
        for (quantile in 1:3) {
          x = c(x, assessments[[e]][q, quantile])
        }
      }
      DMx[[q]] = c(L[q], unique(sort(unlist(x))), U[q])
    }
  }

  tmpDMy = list()
  isItemWeights = NCOL(weights_) > 1

  for (q in 1:N) {
    tmpDMy[[q]] = vector()
    # DMx[[q]] also contains bounds of intrinsic range, those will have y=0 and y=1
    for (x in DMx[[q]][c(-1,-length(DMx[[q]]))]) {
      y = 0
      for (e in 1:Nexperts) {
        tmp = assessments[[e]][q,]
        tmp = unname(tmp)
        if (isItemWeights) {
          weight = weights_[q, e]
        } else {
          weight = weights_[e]
        }
        if(any(diff(t(tmp))<=0))browser()
        y = y + unname(weight) * approx(c(L[q], tmp, U[q]), c(0, quantiles, 1), x, n=numSamplePoints)$y
        #we are just aggregating at interpolation knots - where the cdf changes slope; later on we can sample
        #as many points we want once the interpolation knots of the DM are known
      }
      tmpDMy[[q]] = c(tmpDMy[[q]], y)
    }
  }

  # add boundary values to y
  for (q in 1:N) {
    tmpDMy[[q]] = c(0, tmpDMy[[q]], 1)
  }
  
  #now sample from the DM distribution - swap x and y to sample from x
  samples_DM_x = approx(samples_DM_y, samples_x, n=10000)$y
  samples_DM_y_plot = approx(samples_DM_y,samples_x, n=10000)$x
  
  # get DMs' distribution
  # use approx with x and y swapped to input the y and sample the x 
  
  DM_distr = list()
  
  for (q in 1:N) {
    
  DM = matrix(0, nrow=numSamplePoints, ncol=2)
  
      if(any(diff(tmpDMy[[q]])<=0))browser()
      #x samples
      DM[, 1] = approx(tmpDMy[[q]], DMx[[q]], n=numSamplePoints)$y
      #y samples
      DM[, 2] = approx(tmpDMy[[q]], DMx[[q]], n=numSamplePoints)$x

  DM_distr[[q]] = DM
  }

  return (DM_distr)
}


