#Kelp microbiome review

library(metagear)
library(RefManageR)
library(PRISMAstatement)

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



ReadZotero(user = "8204423", .params=list(q = "kelp", key = "QJbKC91HmGH0hwnxNXOZHOFH",
                                          collection = "LMELYQTG", limit=5))
