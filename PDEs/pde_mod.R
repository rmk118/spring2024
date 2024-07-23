# Ruby Krasnow
# Code for PDE final project
# Last modified: April 29, 2024

#load packages
library(tidyverse)
library(gt)
library(sf)
library(spocc)
library(PNWColors)
library(leaflet)
library(webshot2)
library(ggmap)
library(rstatix)


# Predictions -------------------------------------------------------------
colony_day <- 0.225 #0.225 km per colony/day
raft_prob <- (2/3) # probability that a given colony will raft
hab_prob <- 0.5 # probability that a given rafting colony will reach a suitable habitat
km_year_raft <- colony_day*raft_prob*hab_prob*365
D_raft <- 2*(km_year_raft^2)/(pi)

r <- c(1.887,1.842,1.665,1.512,0.972,0.984, 0.867,0.778)
test_v <- function(x) {
  2*sqrt(x*D_raft)
}
preds <- map_vec(r, test_v)


# Observations ------------------------------------------------------------

# Define spatial boundary for observations
polygon_list = list(rbind(c(-72, 41), c(-59, 41), c(-59, 50), c (-72, 50), c(-72, 41)))
poly<- st_polygon(polygon_list)
poly<-st_sfc(poly, crs=4326)

# Query GBIF database
df_gbif <- occ(query = 'Botrylloides violaceus',
               from = 'gbif',
               has_coords = TRUE,
               geometry = c(-72, 41, -59, 50),
               limit = 1000,
               gbifopts = list(basisOfRecord ="HUMAN_OBSERVATION"))

# Convert lon and lat to numbers, remove observations with no date
df_gbif_cleaned <- occ2df(df_gbif) %>% 
  mutate(lng=as.double(longitude), lat=as.double(latitude), .keep="unused") %>% 
  filter(!is.na(date))

# Convert df to sf
tunicate <- st_as_sf(df_gbif_cleaned %>% 
                    select(date, lng, lat) %>% # Select relevant columns
                    mutate(year = year(date)),  # Add year column
                    coords = c("lng","lat")) # Specify which cols are coordinates

years <- tunicate %>% group_by(year) %>% summarise(year=mean(year)) # Group by year of obs
 
# Define function to get and merge all obs from previous years
get_yrs <- function(yr) { 
   out <- years %>% filter(year <= yr) %>% st_union()
   }
 
years_vec <- c(2005:2024) # Years with observations in the GBIF dataset

# Apply function to all years in the years_vec vector
years2 <- years %>% mutate(years_agg = map_vec(year, get_yrs))

init_pt <- tunicate %>% 
  st_set_crs(4326) %>% # Define CRS
  filter(year==2005) %>% # Select only the year of introduction
  pull(geometry) # Pull the coordinates

dist <- tunicate %>% 
  filter(year>2005) %>% # Remove year of introduction
  st_set_crs(4326) %>% # Define CRS
  ungroup() %>% 
  mutate(dist = st_distance(geometry, init_pt)) %>% # Find the distance b/w the initial point and a given observation
  group_by(year) %>% 
  st_drop_geometry() %>% 
  select(year, dist) %>%  # Select only year and distance columns
  summarize(max_dist = max(dist)) %>% # Find the maximum distance b/w observations and initial point in each year
  mutate(max_dist=cummax(max_dist)) # Maximum distance the species has been observed from the initial point before and during a given year

speed <- dist %>% mutate(km=as.double(max_dist/1000), # Convert from m to km
                         speed=km/(year-2005)) %>% # Speed = distance (km)/time (yr)
                         na.omit() %>% 
  mutate(outlier = is_outlier(speed), extreme=is_extreme(speed))

speed2 <- speed %>% 
  filter(year>2006)

# Tables -----------------------------------------------------------------------

## Table 1 ----------------------------------------------------

params <- tribble(
  ~Parameter, ~Value,
  "Distance traveled per colony per day (km)", 0.225,
  "Probability that a given colony will raft", 2/3,
  "Probability that a given rafting colony will reach a suitable habitat", 0.5,
  "Mean displacement (km/year)", 27.38,
  "D (km\U00B2/year)", 477.1
)

table1 <- gt(params)  %>% 
  fmt_number(
    columns = Value,
    decimals = 3,
    use_seps = FALSE,
    drop_trailing_zeros = TRUE
  ) %>% 
  cols_width(Parameter ~ pct(78), Value ~pct(22)) %>% 
  tab_options(table.font.size =12)

## Table 2 -------------------------------------------------------

obs <- c(summary(speed %>% filter(year>2006) %>% pull(speed)))
preds <- c(summary(preds))

summaries_df <-data.frame("Observed"=obs,
                 "Predicted"=preds) %>% 
  rownames_to_column() %>% 
  gt() %>% 
  fmt_number(
    columns = c(Observed, Predicted),
    decimals = 1,
    use_seps = FALSE)


# Figures -----------------------------------------------------------------

## Fig. 1 -------------------------------------------------------

# Create a continuous palette function
pal <- colorNumeric(
  palette = "viridis",
  domain = years_vec)

radii <- speed %>% select(year, max_dist, km) %>% arrange(-year)

map <- leaflet(options = leafletOptions(zoomControl = FALSE)) %>% 
  fitBounds(-72,41,-59.52,50) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addCircles(data=tunicate, color="black") %>% 
  addCircles(data=radii, 
             lng=-70.03979, 
             lat=43.75033, 
             radius=~max_dist, color = ~pal(year), fillOpacity = 0.1) %>%
  addLegend(data=radii, "bottomright", pal = pal, values = ~year,
            title = "Year",
            labFormat = labelFormat(big.mark = ""),
            opacity = 1)
map


## Fig. 2 -------------------------------------------------------

fig2 <- ggplot()+
  geom_histogram(aes(x=speed2 %>% pull(speed)), bins=23,
                 color=pnw_palette("Starfish")[2], fill=pnw_palette("Starfish")[2], alpha=0.5)+
  geom_vline(aes(xintercept = min(preds)),lty="dashed", color="black")+
  geom_vline(aes(xintercept = max(preds)),lty="dashed", color="black")+
  labs(x="Observed invasion speed (km/year)", y="Count", lty=NULL)+
  lims(x=c(0,NA))+
  theme_classic()+
  scale_linetype(guide=FALSE)+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        text = element_text(size = 13))
