library(dplyr)
library(readxl)

# B0 - baseline
# D1 - new start
# D5 - week 5
# F3 - month 3
# F6 - month 6
# F9 - month 9
# F12 - moth 12

# Need to remove D1A and D1B
# Change names of clean df
# add a IF and NIF columns to the metadata

metaData <- read.csv("~/Roselab/Metabolite/data/data_for_analysis/Original/meta_data_human.csv")
metlData <- read_excel("~/Roselab/Metabolite/data/data_for_analysis/Original/all_sample_data_clean_77.xlsx")

labels_to_remove <- c("D1A", "D1B", "D7")
cols_to_remove <- metaData[metaData$Time %in% labels_to_remove, ]$Sample_ID

metaData <- metaData %>%
  filter(!Time %in% labels_to_remove) %>%
  mutate(Arm = ifelse(Arm == "Intermittent Fasting", "IF", "NIF"))

metlData <- metlData %>%
  rename_with(~ gsub("-", "_", .x)) %>%
  mutate(ID = paste(Index, Compounds, sep = "|")) %>%
  select(-Index, -Compounds, -all_of(cols_to_remove)) %>%
  select(ID, everything())


# setdiff(colnames(metlData), metaData$Sample_ID) # only ID is the diff, good

write.csv(metaData, "~/Roselab/Metabolite/data/data_for_analysis/meta_data_74.csv", row.names = FALSE)
write.csv(metlData, "~/Roselab/Metabolite/data/data_for_analysis/metabolite_data_74.csv", row.names = FALSE)
