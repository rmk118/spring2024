

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

canada <- read_csv("./PDEs/sephton2011.csv",
                   col_names=TRUE,
                   col_types="cddllll")
canada <- canada %>% mutate(lng = -1*lon) %>% mutate(across(where(is.logical),~replace_na(.x, FALSE)))

canada <- canada %>% 
  pivot_longer(cols=where(is.logical), 
               values_to = "present", names_to = "year") %>%  filter(present==TRUE)%>%
  mutate(year=as.numeric(year)) %>% group_by(site) %>% slice_min(year) %>% ungroup()

# Convert df to sf
ca_tunicate <- st_as_sf(canada %>%
                       select(year, lng, lat), # Select relevant columns
                     coords = c("lng","lat")) # Specify which cols are coordinates

pal_ca <- colorNumeric(
  palette = c(pnw_palette("Sailboat", 4, type = "discrete")),
  domain = c(2006, 2007, 2008, 2009))


leaflet(data=ca_tunicate) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addCircles(color = ~pal_ca(year)) %>%
  addLegend("bottomright", pal = pal_ca, values = ~year,
            title = "Year",
            labFormat = labelFormat(big.mark = ""))
