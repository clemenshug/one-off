library(tidyverse)
library(readxl)

index_barcodes <- read_csv("example_sample_sheet.csv", skip = 22) %>%
  select(Index_Plate_Well, I7_Index_ID, index, I5_Index_ID, index2)

plate_layout <- read_excel("TruSeq UDI plate index.xlsx", skip = 1) %>%
  column_to_rownames("...1") %>%
  as.matrix()

# Indices in rows F-H where supposed to be used, but instead the plate
# was rotate 180° and A-C where used instead. Mapping the intended to the
# actual indices by rotating the index matrix

# This rotates matrix 180°
rot_layout <- rev(plate_layout) %>%
  `dim<-`(dim(plate_layout))

intended_actual_mapping <- set_names(c(rot_layout), c(plate_layout))

indended_ss <- read_csv("sample_sheet_FC_05924.csv", skip = 22)

actual_ss <- indended_ss %>%
  select(Sample_ID, Sample_Well, I7_Index_ID) %>%
  mutate(I7_Index_ID = intended_actual_mapping[I7_Index_ID]) %>%
  left_join(index_barcodes, by = "I7_Index_ID")

clipr::write_clip(actual_ss)
