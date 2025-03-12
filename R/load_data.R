
#######################################################################
##specify the directory where R works##################################
#######################################################################

#import the data in excel format, all studies in different sheets - data has 49 studies
#sheet_names <- excel_sheets("Data/Expert_data_60_studies.xlsx")
# Get sheet names
#new data with 60 studies
#the data format is changed- there are 2 new empty columns (3&5) introduced
#sheet_names = excel_sheets("Data/Expert Data 2006_2022.xlsx")
#sheet_names

load_data_formatted <- function(excel_file_name) {
  sheet_names <- readxl::excel_sheets(excel_file_name)

  list_all <- lapply(sheet_names, function(x) {
    as.data.frame(readxl::read_excel(excel_file_name, sheet = x)) } )

  list_all<-list_all[-1] # Remove _xltb_storage_ sheet
  sheet_names<-sheet_names[-1]
  names(list_all) <- sheet_names

  Nquantiles = 3
  quantilesPercent = c("5%","50%","95%")
  quantiles = c(0.05, 0.50, 0.95)

  check_data(list_all, sheet_names, Nquantiles)

  # lapply format_data to all studies
  lapply(seq_along(sheet_names), function(x) format_data(list_all, x))
}

#' Combine the nested formatted data into one list of data frames
#'
#' @param formatted_data A list of lists of data frames.
#'
#' @returns A single list of data frames.
#' @export
#'
combine_expert_lists <- function(formatted_data) {
  # combine all the data into one data frame
  formatted_data |> purrr::map(dplyr::bind_rows)
}

#################################################################################################################
###get result for all the studies################################################################################
#################################################################################################################
#make sure the data is in the correct format
check_data <- function(my_data, sheet_names, Nquantiles) {
  check_all <- rep(0, 3)

  for (i in 1:length(unique(sheet_names)))
  {
    data_i <- my_data[[i]]

    # get the vector of realizations for 49 studies
    realizations <- data_i[, 6]
    # for 60 studies
    # realizations = data_i[,8]
    # get the list of questions
    questions <- unique(data_i[, 2])

    m <- list() # create empty list of matrices for assessment data

    # read experts assessments
    badReads <- 0
    goodReads <- 0

    expertMappingRev <- array()

    for (j in unique(data_i[, 1]))
    {
      tryCatch(
        {
          # get expert ID
          expertID <- as.integer(unique(data_i[, 1])[j])

          expertMappingRev[goodReads + 1] <- expertID

          # get each expert's assessments for 49 studies
          data <- data_i[, 3:5]
          # get each expert's assessments for 60 studies
          # data = data_i[,4:6]
          rows <- which(data_i[, 1] == expertID)

          # check if assessments are strictly increasing
          for (row in rows) {
            if (any(diff(t(data[row, ])) <= 0)) {
              errorString <- sprintf(paste0("Error assessments for %s question %d are not strictly increasing in study", sheet_names[i]), expertID, row)
              stop(errorString)
            }
          }

          # for 49 studies
          mm <- as.matrix(data_i[which(data_i[, 1] == expertID), 3:5], nrow = NROW(data_i[which(data_i[, 1] == expertID), 3:5]), ncol = Nquantiles)
          # for 60 studies
          # mm = as.matrix(data_i[which(data_i[,1]==expertID),4:6], nrow=NROW(data_i[which(data_i[,1]==expertID),4:6]), ncol=Nquantiles)
          if (typeof(mm) == "double" || typeof(mm) == "integer") {
            m[[goodReads + 1]] <- mm
            goodReads <- goodReads + 1
          } else {
            print(paste0("Error cannot interpret data as numeric in study ", sheet_names[i]))
            badReads <- badReads + 1
          }
        },
        warning = function(cond) {
          print(cond)
          message("Warning cannot read csv file")
          badReads <- badReads + 1
        }
      )
    }

    check_i <- c(sheet_names[i], goodReads, badReads)

    check_all <- rbind(check_all, check_i)
  }

  colnames(check_all) <- c("Study", "Good reads", "Bad reads")
  return(check_all[-1, ])
}


##data in the right format
#function to do that per study
format_data<-function(my_data,i)
{

  data_i = my_data[[i]]

  #get the vector of realizations for 49 studies
  realizations = data_i[,6]
  # for 60 studies
  #realizations = data_i[,8]
  #get the list of questions
  questions = unique(data_i[,2])

  m = list() # create empty list of matrices for assessment data

  for (j in unique(data_i[,1]))
  {
    # get expert ID
    expertID = as.integer(unique(data_i[,1])[j])

    rows = which(data_i[,1]==expertID)

    if (expertID == 1) { # only the first experts has the realizations
      realizations = data_i[rows,6]
      # for 60 studies
      #realizations = data_i[,8]
    }

    #get each expert's assessments for 49 studies and the realisation
    data = cbind(data_i[rows,3:5], realizations, data_i[rows,1:2], i)
    #for 60 studies
    #data = data_i[,4:6]

    colnames(data) <- c("5th percentile", "50th percentile", "95th percentile", "realization", "expert_id", "question_id", "study_id")
    #for 49 studies
    #for 60 studies
    #mm = as.matrix(data_i[which(data_i[,1]==expertID),4:6], nrow=NROW(data_i[which(data_i[,1]==expertID),4:6]), ncol=Nquantiles)
    m[[j]] = data
  }

  return(m)
}
