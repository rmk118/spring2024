# Ruby Krasnow
# Code for PDE final project

#load packages
library(spocc)
library(sf)
library(ggmap)

#define spatial boundary for observations
polygon_list = list(rbind(c(-72, 41), c(-59.52, 41), c(-59.52, 50), c (-72, 50), c(-72, 41)))
poly<- st_polygon(polygon_list)
poly<-st_sfc(poly, crs=4326)


# GBIF map
df_gbif <- occ(query = 'Botrylloides violaceus', from = 'gbif', has_coords = TRUE, geometry = poly, limit = 1000)
df_gbif_cleaned <- occ2df(df_gbif) %>% mutate(lon=as.numeric(longitude), lat=as.numeric(latitude)) %>% filter(!is.na(date))

tunicate_gbif <- df_gbif_cleaned %>% 
  mutate(lon=as.numeric(longitude), lat=as.numeric(latitude), date=ymd(date), .keep="unused")

bbox_gbif <- make_bbox(lon, lat, data = df_gbif_cleaned)
map_gbif <- get_stadiamap(bbox = bbox_gbif, maptype = "stamen_toner_lite", zoom=6)

ggmap(map_gbif) +
  geom_point(data = df_gbif_cleaned, aes(color=as.factor(year(date))))

# create circles with expanding radii by year
# mean inc. in radius/year ~ spreading speed