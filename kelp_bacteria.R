#Kelp microbiome review

library(metagear)
#library(RefManageR)
library(PRISMAstatement)
library(bib2df)

prisma(found = 5933,
       found_other = 106,
       no_dupes = 776, 
       screened = 776, 
       screen_exclusions = 13, 
       full_text = 763,
       full_text_exclusions = 17, 
       qualitative = 746, 
       quantitative = 319,
       width = 800, height = 800)

all_refs <- bib2df("retrieved_biome.bib")
