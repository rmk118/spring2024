# Ruby Krasnow
# Code for PDE final project
# Last modified: April 15, 2024

#load packages
library(tidyverse)
library(spocc)
library(sf)
library(lwgeom)
library(leaflet)

#define spatial boundary for observations
polygon_list = list(rbind(c(-72, 41), c(-59.52, 41), c(-59.52, 50), c (-72, 50), c(-72, 41)))
poly<- st_polygon(polygon_list)
poly<-st_sfc(poly, crs=4326)


# GBIF map
df_gbif <- occ(query = 'Botrylloides violaceus',
               from = 'gbif', has_coords = TRUE,
               geometry = c(-72, 41, -59.52, 50), limit = 1000)

df_gbif_cleaned <- occ2df(df_gbif) %>% 
  mutate(lng=as.numeric(longitude), lat=as.numeric(latitude), .keep="unused") %>% 
  filter(!is.na(date))

tunicate <- st_as_sf(df_gbif_cleaned %>% 
                       select(date, lng, lat) %>% 
                       mutate(year = year(date)),  
                     coords = c("lng","lat"))

 years <- tunicate %>% group_by(year) %>% summarise(year=mean(year))
 

 get_yrs <- function(yr) {
   out <- years %>% filter(year <= yr) %>% st_union()
   }
 
 
years_vec <- c(2005:2024)

years2<- years %>% mutate(years_agg = map_vec(year, get_yrs))

circles <- years2 %>%
  rowwise() %>%
   mutate(circle = st_minimum_bounding_circle(years_agg), .keep="unused") %>% ungroup()

circles2 <- circles %>% st_drop_geometry() %>% 
  select(year, circle)

circles3 <- st_as_sf(circles2) %>% arrange(-year)

ggplot()+geom_sf(data=circles3)


# Create a continuous palette function
pal <- colorNumeric(
  palette = "viridis",
  domain = circles3$year)

m <- leaflet(circles3) %>% fitBounds(-72,41,-59.52,50)

m %>% 
  addProviderTiles(providers$Stadia.StamenTonerLite) %>% 
  addCircles(data=tunicate, color="black") %>% 
  addPolygons(fillOpacity = 0.2,
              color = ~pal(year))


