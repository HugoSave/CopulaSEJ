#performance-based item weights function that would provide the list of weights corresponding to the
#best performing (combined score) DM

itemWeights_opt<-function(assessments,realizations)
{
  
  # check all alpha cutoff values
  bestAlpha = 0.0
  bestCalScore = 0.0
  bestInfoScore = 0.0
  bestCombScore = 0.0
  bestWeights = 0
  alpha = 0.0   # cutoff value
  
  Nexperts = length(assessments)      # number of experts
  N = NROW(assessments[[1]])
  Ncal = length(which(!is.na(realizations)))
  realizations1 = realizations[which(!is.na(realizations))] # we need to account only for cal vbls in computing the intrinsic range
  Nquantiles=ncol(assessments[[1]])
  
  calibrationScores = calculateCalibrationScore(assessments, realizations)
  tmp = calculateInformationScore(assessments, realizations, k=0.1, bounds=NULL)
  informationScores = tmp[[1]]#[which(!is.na(realizations)),]
  L = tmp[[2]]#[which(!is.na(realizations))]
  U = tmp[[3]]#[which(!is.na(realizations))]
  
  combinedScores = matrix(0, nrow=N, ncol=Nexperts)
  for (q in 1:N) {
    combinedScores[q,] = calibrationScores * informationScores[q,]
  }
  
  itemWeights = matrix(0, nrow=N, ncol=Nexperts)
  for (q in 1:N) {
    itemWeights[q,] = calibrationScores * informationScores[q,]
    itemWeights[q,] = itemWeights[q,] / sum(itemWeights[q,])
  }
  
  
  orderedCalScores = unique(sort(unname(calibrationScores)))
  if (orderedCalScores[1]==0)
  {
    orderedCalScores = orderedCalScores[-1]
  }
  numAlphas = length(orderedCalScores)
  
  for (a in 1:numAlphas) {
    alpha = orderedCalScores[a]
    expertsForDM = calibrationScores >= alpha
    tmpWeights = sapply(1:ncol(combinedScores),function(x) combinedScores[,x] * (calibrationScores >= alpha)[x] )
    tmpWeights = as.matrix(t(sapply(1:nrow(tmpWeights),function(x) tmpWeights[x,]/sum(tmpWeights[x,]))))
    
    # create optimized DM using best Alpha
    tmpDM = constructDM(assessments[expertsForDM], tmpWeights[, colSums(tmpWeights != 0) > 0], NULL, L, U, quantiles)
    
    # check DM scores
    tmpCalScore = calculateCalibrationScore(list(tmpDM), realizations)
    tmpInfoScore = tmpInfoScore = mean(calculateInformationScore(list(tmpDM), realizations, k=0.1, bounds=NULL, L, U)[[1]])
    
    tmpCombScore = tmpCalScore * tmpInfoScore
    
    if (tmpCombScore > bestCombScore) {
      bestAlpha = alpha
      bestCalScore = tmpCalScore
      bestInfoScore = tmpInfoScore
      bestWeights = tmpWeights
      bestCombScore = tmpCombScore
      bestDM = tmpDM
    }}
  
  return(list(bestWeights,bestDM))
}