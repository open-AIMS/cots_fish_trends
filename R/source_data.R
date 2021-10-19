source("_drake.R")
library(writexl)

drake::loadd(qdaf_models)
fig_1_data <- make_fig_1(qdaf_models)$data
rm(qdaf_models)

drake::loadd(fish_models)
fig_2_data <- make_fig_2(fish_models, return_data = TRUE)
rm(fish_models)

drake::loadd(cots_model)
fig_3_data <- make_fig_3(cots_model, return_data = TRUE)
rm(cots_model)

dir.create("output/source_data", recursive = TRUE, showWarnings = FALSE)
writexl::write_xlsx(fig_1_data, "output/source_data/fig_1_data.xlsx")
writexl::write_xlsx(fig_2_data, "output/source_data/fig_2_data.xlsx")
writexl::write_xlsx(fig_3_data, "output/source_data/fig_3_data.xlsx")
zip(zipfile = "output/source_data", files = "output/source_data")
