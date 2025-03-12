
# Create a list of all possible folds of the data
# data: a data frame
# fold_size: the size of each fold
# return: a list of lists with elements "training" and "test" that are data 
# frames. The outer list has length choose(n, fold_size).
ExhaustiveFolds <- function(data, fold_size){
  
  n = nrow(data)
  indices = seq(n)
  
  # select all possible combinations of fold_size indices
  fold_combinations = combn(indices, fold_size)
  
  # create a list of data frames
  training_sets = list()
  test_sets = list()
  for (i in 1:ncol(fold_combinations)){
    test_sets[[i]] = data[fold_combinations[,i],]
    training_sets[[i]] = data[-fold_combinations[,i],]
  }
  
  # create tibble from list with column names training and test
  folds = tibble(training = training_sets, test = test_sets)
  
  
  return(folds)
}
