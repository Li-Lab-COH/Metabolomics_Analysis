library(dplyr)
F3_data <- read.csv("~/Roselab/Metabolite/results/difference/limma/intraArm/limma_intraArm_F3_IF_vs_NIF.csv")
D5_data <- read.csv("~/Roselab/Metabolite/results/difference/limma/intraArm/limma_intraArm_D5_IF_vs_NIF.csv")

F3_data <- subset(F3_data, P.Value <= 0.05)
D5_data <- subset(D5_data, P.Value <= 0.05)

common_indices <- intersect(F3_data$Index, D5_data$Index)


F3_data_common <- F3_data %>% 
  filter(Index %in% common_indices)

D5_data_common <- D5_data %>% 
  filter(Index %in% common_indices)

write.csv(F3_data_common, "~/Roselab/Metabolite/results/difference/common_D5_F3/F3_common_metabolite_results.csv")
write.csv(D5_data_common, "~/Roselab/Metabolite/results/difference/common_D5_F3/D5_common_metabolite_results.csv")
