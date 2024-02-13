#Kelp microbiome review

library(PRISMAstatement)
library(tidyverse)

prisma(found = 1448,
       found_other = 106,
       no_dupes = 1483, 
       screened = 1483, 
       screen_exclusions = 1308, 
       full_text = 175,
       full_text_exclusions = 125, 
       qualitative = 50, 
       quantitative = 50,
       width = 800, height = 800)

#files for full-text screening
inital_screen <- read_csv("./kelp_sysrev/initial_screen.csv")

