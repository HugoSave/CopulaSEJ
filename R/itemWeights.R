# item based weights
#note that this function is different - is is customed to SEJ data what has all experts assessments on one sheet, and where
#only calibration questions assessments are being analyzed 
#######################################################
itemWeights<-function(assessments, realizations)
{
  Nexperts = length(assessments)      # number of experts
  N = NROW(assessments[[1]])

  #which are calibration variables
  Cal = which(!is.na(realizations))
  # which are the variables of interest
  #QoI = which(is.na(realizations))
  #nQoI = N-length(Cal)
  
  # create labels for all questions
  labelsQuestions = c()
  for (q in Cal) {
    labelsQuestions = c(labelsQuestions, sprintf("CAL%d", q))
  }
  
  # create labels for all experts
  labelsExperts = c()
  for (e in 1:Nexperts) {
    labelsExperts = c(labelsExperts, sprintf("Expert %d", e))
  }  

  
  for (e in 1:Nexperts) {
    dimnames(assessments[[e]][Cal,]) = list(labelsQuestions, quantilesPercent)
  }
  
  calibrationScores = calculateCalibrationScore(assessments,realizations)
  
  tmp = calculateInformationScore(assessments, realizations, k=0.1, bounds=NULL, Larg = NULL, Uarg = NULL)
  
  informationScores = tmp[[1]]
  L = tmp[[2]]
  U = tmp[[3]]
  
  itemWeights = matrix(0, nrow=N, ncol=Nexperts)
  dimnames(itemWeights) = list(labelsQuestions, labelsExperts)
  for (q in 1:N) {
    itemWeights[q,] = calibrationScores * informationScores[q,]
    itemWeights[q,] = itemWeights[q,] / sum(itemWeights[q,])
  }
  
  return(itemWeights)
}
