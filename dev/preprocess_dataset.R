library(CopulaSEJ)
library(dplyr)
source("dev/load_data.R")


file_name <- "dev/data/Expert Data Nov 24 2021.xlsx"
data_list_form <- load_data_formatted(file_name)
studies <- combine_expert_lists(data_list_form)
saveRDS(studies, "dev/output/data49_nov24.rds")
studies <- load_data_49(relative_dev_folder = FALSE)

remove_high_magintude_difference <- function(studies, max_mag_ratio=100) {
  studies_filtered <- lapply(studies, function(study) {
    study_ranges <- study |> mutate(expert_range = `95th percentile` -`5th percentile`)
    filtered_questions <- study_ranges |> group_by(question_id) |>
      mutate(range_mag_diff = max(expert_range) / min(expert_range)) |> ungroup() |> filter(range_mag_diff < max_mag_ratio)
    return(filtered_questions)
  })
  return(studies_filtered)
}

find_high_magnitude_differences <- function(studies, max_mag_ratio = 100) {
  to_remove <- lapply(studies, function(study) {
    study_id_val <- unique(study$study_id)

    study_ranges <- study %>%
      mutate(expert_range = `95th percentile` - `5th percentile`) %>%
      group_by(question_id) %>%
      mutate(range_mag_diff = max(expert_range) / min(expert_range)) %>%
      ungroup() %>%
      filter(range_mag_diff >= max_mag_ratio) %>%
      distinct(question_id) %>%
      mutate(study_id = study_id_val)

    return(study_ranges %>% select(study_id, question_id))
  })

  # Combine all the to-be-removed entries into a single tibble
  bind_rows(to_remove)
}
