library(CopulaSEJ)

file_name <- "dev/data/Expert Data Nov 24 2021.xlsx"
data_list_form <- load_data_formatted(file_name)
studies <- combine_expert_lists(data_list_form)
saveRDS(studies, "dev/output/data49_nov24.rds")

