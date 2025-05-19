
# this function calculates the calibration score and returns it
calculateCalibrationScore = function (assessments, realizations) {

  Nexperts = length(assessments)      # number of experts
  Ncal = length(which(!is.na(realizations)))# number of calibration questions
  #length(realizations) # this is not working, as there are also the QoI in the realization file
#  cat("Nexperts: ", Nexperts, ", Ncal: ", Ncal, "\n")

  # check realization against assessments using bucketing
  df = matrix(0, nrow=Nexperts, ncol=4)

  for (e in 1:Nexperts) {
    for (q in 1:Ncal) {
      if (realizations[q] <= assessments[[e]][q,][1]) {
        df[e, 1] = df[e, 1] + 1
      } else if (realizations[q] <= assessments[[e]][q,][2]) {
        df[e, 2] = df[e, 2] + 1
      } else if (realizations[q] <= assessments[[e]][q,][3]) {
        df[e, 3] = df[e, 3] + 1
      } else {
        df[e, 4] = df[e, 4] + 1
      }
    }
  }

  # A+B+C might not equal C+B+A with doubles
  # so swap columns (1,4) and (2,3) such that the lowest value has the lowest index
  t = pmin(df[,1], df[,4])
  df[,4] = pmax(df[,1], df[,4])
  df[,1] = t
  t = pmin(df[,2], df[,3])
  df[,3] = pmax(df[,2], df[,3])
  df[,2] = t

  # empirical probability vector, normalized
  se = list(array(0, dim=Nexperts))
  for (e in 1:Nexperts) {
    se[[e]] = c(as.vector(df[e,]))
    se[[e]] = se[[e]] / sum(se[[e]])
  }

  # theoretical probability vector
  p = c(0.05, 0.45, 0.45, 0.05)

  # Kullback-Leibler divergence
  # l(s,p) = sum i=1:4 s_i * ln(s_i / p_i)

  l = array(0, dim=Nexperts)
  for (e in 1:Nexperts) {
    for (i in 1:4) {
      # log(0) does not exist, so check argument first
      if (se[[e]][i] > 1e-6) {
        l[e] = l[e] + se[[e]][i] * log(se[[e]][i] / p[i])
      }
    }
  }

  # calibration score
  # Cal(e) = 1 - F(2 * m * l(s,p))

  calibrationScores = 1 - pchisq(2 * Ncal * l, df=3)

  #print(calibrationScores)
  return(calibrationScores)
}

calculateCalibrationScoreForExpert <- function(assessments, realizations) {
  # Number of calibration questions
  Ncal <- length(which(!is.na(realizations)))

  # Initialize count buckets for realizations
  counts <- c(0, 0, 0, 0)

  # Check realization against assessments using bucketing
  for (q in 1:Ncal) {
    if (realizations[q] <= assessments[q, 1]) {
      counts[1] <- counts[1] + 1
    } else if (realizations[q] <= assessments[q, 2]) {
      counts[2] <- counts[2] + 1
    } else if (realizations[q] <= assessments[q, 3]) {
      counts[3] <- counts[3] + 1
    } else {
      counts[4] <- counts[4] + 1
    }
  }

  # Adjust counts for ordering
  counts[c(1, 4)] <- sort(counts[c(1, 4)])
  counts[c(2, 3)] <- sort(counts[c(2, 3)])

  # Empirical probability vector, normalized
  se <- counts / sum(counts)

  # Theoretical probability vector
  p <- c(0.05, 0.45, 0.45, 0.05)

  # Kullback-Leibler divergence
  kl_divergence <- sum(ifelse(se > 1e-6, se * log(se / p), 0))

  # Calibration score
  calibrationScore <- 1 - pchisq(2 * Ncal * kl_divergence, df = 3)

  return(calibrationScore)
}

