

gbif_orig <- df_gbif[["gbif"]][["data"]][["Botrylloides_violaceus"]]


years_na <- occ2df(df_gbif) %>% 
  mutate(lng=as.double(longitude), lat=as.double(latitude), .keep="unused") %>% 
  mutate(year=year(date))


with_old <- years_na %>% 
  left_join(gbif_orig %>% select(key, year), by="key") %>% 
  rename(year=year.y) %>% 
  select(-year.x) %>% mutate(year = as.numeric(year))

tunicate_old <- st_as_sf(with_old %>% 
                       select(year, lng, lat),coords = c("lng","lat")) %>% filter(year>1980)

years_old <- tunicate_old %>% group_by(year) %>% summarise(year=mean(year)) # Group by year of obs

years_vec_old <- c(1988:2024) # Years with observations in the GBIF dataset

get_yrs_old <- function(yr) { 
  out <- years_old %>% filter(year <= yr) %>% st_union()
}
years2_old <- years_old %>% mutate(years_agg = map_vec(year, get_yrs_old))

init_old <- tunicate_old %>% 
  st_set_crs(4326) %>% # Define CRS
  filter(year==1988) %>% # Select only the year of introduction
  pull(geometry) # Pull the coordinates

dist_old <- tunicate_old %>% 
  filter(year>1988) %>% # Remove year of introduction
  st_set_crs(4326) %>% # Define CRS
  ungroup() %>% 
  mutate(dist = st_distance(geometry, init_pt)) %>% # Find the distance b/w the initial point and a given observation
  group_by(year) %>% 
  st_drop_geometry() %>% 
  select(year, dist) %>%  # Select only year and distance columns
  summarize(max_dist = max(dist)) %>% # Find the maximum distance b/w observations and initial point in each year
  mutate(max_dist=cummax(max_dist)) # Maximum distance the species has been observed from the initial point before and during a given year


speed_old <- dist_old %>% mutate(km=as.double(max_dist/1000), # Convert from m to km
                         speed=km/(year-1988)) %>% # Speed = distance (km)/time (yr)
  na.omit() %>% 
  mutate(outlier = is_outlier(speed), extreme=is_extreme(speed))


# additional data about timing of pre-2010s observations of B. violaceus:
# [1] Carman et al., 2010: An initial assessment of native and invasive tunicates in shellfish aquaculture of the North American east coast,  https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1439-0426.2010.01495.x

# [2] Carlton, 1989: Man's Role in Changing the Face of the Ocean: Biological Invasions and Implications for Conservation of Near-Shore Environments https://www.jstor.org/stable/2386170 - misidentified

# [3] Berman, 1992, Recent Invasions of the Gulf of Maine: Three Contrasting Ecological Histories, https://doi.org/10.1046/j.1523-1739.1992.06030435.x - misidentified

# [4] Bock 2010, Looking at both sides of the invasion: patterns of colonization in the violet tunicate Botrylloides violaceus, https://doi.org/10.1111/j.1365-294X.2010.04971.x 


# 44°56'47.6"N 66°54'19.7"W in 2009, from Martin et al. (2011): doi: 10.3391/ai.2011.6.4.05  

# http://www.aquaticinvasions.net/2011/AI_2011_6_4_Sephton_etal_Supplement.pdf 
# 2006: 43.71370 65.96947
# 2006: 43.72360 65.84057
# 2006: 43.56273 65.36245
# 2006: 44.36113 64.33362
# 2006: 44.37525 64.33069
# 2006: 45.51055 61.01768
# 2006: 45.50252 60.96400
# 2006: 46.62680 61.01597
# 2006: 46.28117 60.42500 
# 2006: 43.54773 65.43217
# 2006: 43.56273 65.36245 

# 2007
# 44.19365 66.16688
# 43.71370 65.96947
# 43.72360 65.84057
# 43.44470 65.63510
# 43.56273 65.36245
# 43.75780 65.32200
# 43.69950 65.10740
# 44.37525 64.33069
# 44.45598 64.30683 
# 44.53780 64.23840
# 45.33450 60.98610
# 46.20680 60.24900
# 46.62680 61.01597 

# 2008
# 44.19365 66.16688
# 43.81648 66.14785
# 43.71370 65.96947
# 43.72360 65.84057
# 43.44470 65.63510
# 43.56273 65.36245
# 43.69950 65.10740 
# 44.37525 64.33069
# 44.45598 64.30683
# 44.44782 64.37437
# 44.53780 64.23840
# 45.50710 60.96050
# 45.58870 60.96180
# 45.58333 60.74067 
# 45.66115 60.87440
# 46.82020 60.15450 
# 46.30850 60.28400
# 46.90320 60.46040
# 47.00010 60.46470 
# 46.62680 61.01597 
# 46.23000 61.31690

# 2009

