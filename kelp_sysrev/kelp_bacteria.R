#Kelp microbiome review
#Ruby Krasnow
#Last modified: March 9, 2024

library(PRISMAstatement)
library(tidyverse)
library(patchwork)
library(clipr)
library(PNWColors)

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


# Core
phyla <- read.csv("./kelp_sysrev/microbiome_phyla.csv")
class_raw <- read.csv("./kelp_sysrev/microbiome_classes.csv")
fam_raw <- read.csv("./kelp_sysrev/microbiome_families.csv") %>% select(-First)
genera_raw <- read.csv("./kelp_sysrev/microbiome_genera.csv") %>% select(-First)

author_list <- c("Staufenberger", "Tourneroche", "Liu", "BD", "King", "Davis", "Schenk", "Park")


# Genera ------------------------------------------------------------------
genera_tidy <- genera_raw %>% 
  select(-Wiese) %>% 
  pivot_longer(cols=all_of(author_list), values_to = "present", names_to = "study")

genera <- genera_tidy %>% 
  group_by(genus) %>% 
  summarise(sum = sum(present, na.rm = TRUE)) %>% filter(sum >0) %>% 
  arrange(-sum)

genera2 <- left_join(genera_raw, genera) %>% 
  filter(!is.na(sum))%>% 
  select(-Wiese) %>% 
  mutate(across(all_of(author_list), ~replace_na(.x, FALSE))) %>% 
  filter(sum > 2) %>% 
  arrange(-sum)

consistent_gen <- genera2$genus

genera_heat <- genera_tidy %>%filter(genus %in% consistent_gen)

study_obs_gen <- genera_heat %>% 
  group_by(study) %>% 
  summarise(sum = sum(present, na.rm = TRUE))

study_gen <- study_obs_gen %>% arrange(-sum) %>% pull(study)


gen_plot <- ggplot(data=genera_heat, aes(x=factor(study, level=study_gen), y=factor(genus, level=consistent_gen), fill=present)) +
  geom_tile()+
  guides(fill="none")+
  theme_classic()+
  labs(x="Study", y="Genus", fill=NULL)+
  scale_fill_manual(values=pnw_palette(name="Sailboat",n=1,type="discrete"), na.value = "#dcdcdc")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
                           text = element_text(size = 13),
        axis.line = element_line(colour = "white"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
        axis.text.x = element_text(angle = -45, hjust=0.1),
        plot.margin = margin(0.2,1,0.1,0.2, "cm"))


# Families ----------------------------------------------------------------
fam_tidy <- fam_raw %>% 
  select(-Wiese) %>% 
  pivot_longer(cols=all_of(author_list), values_to = "present", names_to = "study")

fam <- fam_tidy %>% 
  group_by(family) %>% 
  summarise(sum = sum(present, na.rm = TRUE)) %>% filter(sum >0) %>% 
  arrange(-sum)

fam2 <- left_join(fam_raw, fam) %>% 
  filter(!is.na(sum))%>% 
  select(-Wiese) %>% 
  mutate(across(all_of(author_list), ~replace_na(.x, FALSE))) %>% 
  filter(sum > 2) %>% 
  arrange(-sum)

consistent_fam <- fam2$family

fam_heat <- fam_tidy %>%filter(family %in% consistent_fam)

study_obs_fam <- fam_heat %>% 
  group_by(study) %>% 
  summarise(sum = sum(present, na.rm = TRUE))

study_gen <- study_obs_fam %>% arrange(-sum) %>% pull(study)

fam_plot <- ggplot(data=fam_heat, aes(x=factor(study, level=study_gen), y=factor(family, level=consistent_fam), fill=present)) +
  geom_tile()+
  guides(fill="none")+
  theme_classic()+
  labs(x="Study", y="Family", fill=NULL)+
  scale_fill_manual(values=c("#015b58"), na.value = "#dcdcdc")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        text = element_text(size = 13),
        axis.line = element_line(colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = -45, hjust=0.1),
        plot.margin = margin(0.2,1,0.1,0.2, "cm"))

gen_plot+fam_plot+
  plot_annotation(tag_levels = 'A')


# Wiese -----------------------------------------------------------------

genera_raw %>% filter(genus %in% consistent_gen) %>% view()
fam_raw %>% filter(family %in% consistent_fam & Wiese==TRUE)
