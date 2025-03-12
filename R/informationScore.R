
# this script calculates the information score
# only for 3 given quantiles (5%, 50%, 95%)

# input:
# assessments: list of matrices contained expert assessments
# realizations: list of realizations of calibration questions
#Ncal: the number of considered calibration questions - this is for the dynamic case, when subsets of calibration questions can be considered
#for validation or optimization
# k: overshoot for intrinsic range, default: 0.1
# bounds: lower and upper bounds for intrinsic range for each question

# output: list(infoScores, L, U)
# infoScores: matrix of information scores for each expert and each question
# L: lower bound of intrinsic range per question
# U: upper bound of intrinsic range per question

calculateInformationScore = function (assessments, realizations, k = 0.10, bounds = NULL, Larg = NULL, Uarg = NULL) {

  Nexperts = length(assessments)      # number of experts
  N = NROW(assessments[[1]])
  Ncal = length(which(!is.na(realizations)))
  realizations1 = realizations[which(!is.na(realizations))] #for dynamic case, there might be realizations available for questions which we want to regard as QoI,
  #to investigate performance, for example.
#  cat("Nexperts: ", Nexperts, ", N: ", N, ", k: ", k, "\n")

  # intrinsic range
  if (TRUE == is.null(Larg) && TRUE == is.null(Uarg)) {
    L = c(realizations1, rep(Inf, N-length(realizations1)))
    U = c(realizations1, rep(-Inf, N-length(realizations1)))

    # expert assessments of particular question
    for (q in 1:N) {
      for (e in 1:Nexperts) {
        L[q] = min(L[q], assessments[[e]][q, 1])
        U[q] = max(U[q], assessments[[e]][q, 3])
      }
    }

    R = U - L
    L = L - k*R
    U = U + k*R

    # apply limits to intrinsic range
    if (FALSE == is.null(bounds)) {
      L = pmax(L, bounds[,2][1:N])
      U = pmin(U, bounds[,3][1:N])
    }
  } else {
    L = Larg
    U = Uarg
  }

  # calculate information scores for all experts given
  informationScores = matrix(0, nrow=N, ncol=Nexperts)
  for (e in 1:Nexperts) {
    for (q in 1:N) {
      d1 = assessments[[e]]
      A = 0.05*log(0.05/(d1[q, 1]-L[q]))
      B = 0.45*log(0.45/(d1[q, 2]-d1[q, 1]))
      C = 0.45*log(0.45/(d1[q, 3]-d1[q, 2]))
      D = 0.05*log(0.05/(U[q]-d1[q, 3]))
      E = log(U[q] - L[q])
      informationScores[q, e] = A+B+C+D+E
      if(is.nan(informationScores[q, e])){
        browser()
        informationScores[q, e] = 0.0
      }
    }
  }

  return(list(informationScores, L, U))
}
