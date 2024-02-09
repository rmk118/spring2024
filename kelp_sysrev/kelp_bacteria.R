#Kelp microbiome review

library(PRISMAstatement)
library(tidyverse)

prisma(found = 1448,
       found_other = 106,
       no_dupes = 1483, 
       screened = 776, 
       screen_exclusions = 13, 
       full_text = 763,
       full_text_exclusions = 17, 
       qualitative = 746, 
       quantitative = 319,
       width = 800, height = 800)


screened_refs <- read_csv()

