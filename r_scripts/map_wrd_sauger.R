setwd("/Users/samjohnson/Documents/NHD_H_Wyoming_State_Shape/Shape/")


all_wy_rivers <- readOGR("D:/Masters/misc/rivers/flowline/NHD__Flowline.shp")

install.packages("sf")
install.packages("ggspatial")
install.packages("rgeos") # no go
install.packages("SDraw") # no go
install.packages("rgdal") # no go
install.packages("ggmap")
install.packages("ggsn") # no go
install.packages("prettymapr")
install.packages(c("ggspatial", "prettymapr", "rosm"))
install.packages("maptiles")

library(sf)
library(ggplot2)
library(ggspatial)
library(rgeos)
library(SDraw)
library(rgdal)
library(ggmap)
library(ggsn)
library(prettymapr)
library(rosm)
library(maptiles)

citation("ggmap")
#' @Article{,
#'   author = {Kahle, David and Wickham, Hadley},
#'   title = {ggmap: Spatial Visualization with ggplot2},
#'   journal = {The R Journal},
#'   year = {2013},
#'   volume = {5},
#'   number = {1},
#'   pages = {144--161},
#'   url = {https://journal.r-project.org/archive/2013-1/kahle-wickham.pdf},
#' }

all_wy_rivers_0 <- st_read("NHDFlowline_0.shp")
all_wy_rivers_1 <- st_read("NHDFlowline_1.shp")
all_wy_rivers_2 <- st_read("NHDFlowline_2.shp")
library(dplyr)
all_wy_rivers <- bind_rows(all_wy_rivers, all_wy_rivers_1, all_wy_rivers_2)
roi <- all_wy_rivers[all_wy_rivers$gnis_name %in% c("Wind River","Little Wind River","Popo Agie River", "Little Popo Agie River", "Boysen Reservoir"),]

# test plot for rivers
ggplot() +
  geom_sf(data = roi, color = "blue") +
  theme_minimal()

res <- st_read("NHDWaterbody.shp")
head(table(res$gnis_name))
boysen <- res %>%
  filter(gnis_name == "Boysen Reservoir")
wrd <- bind_rows(roi, boysen)

roi <- st_transform(roi, 4326)
boysen <- st_transform(boysen, 4326)
roi_dissolved <- st_union(roi)


# get topo tiles
tiles <- get_tiles(
  roi,
  provider = "Esri.WorldTopoMap",
  zoom = 10
)

points_df <- data.frame(
  name = c("PA", "LW1", "LW2", "LW3", "LW4", "PCB", "BCB", "WRI"),
  lon = c(-108.50065, -108.43045, -108.43433, -108.41774, -108.41001, -108.161716, -108.165097, -108.181606),
  lat = c(42.94959, 42.97219, 42.97878, 42.98156, 42.98445,  43.248907,  43.284528, 43.178031),
  color = c(viridis::viridis(1), viridis::viridis(1), viridis::viridis(1), viridis::viridis(1), viridis::viridis(1), viridis::plasma(256)[180], viridis::plasma(256)[180], viridis::plasma(256)[180]),
  label = c("Spawning Aggregation (Sampled for F0 and F1 Sauger)", "Spawning Aggregation (Sampled for F0 and F1 Sauger)", "Spawning Aggregation (Sampled for F0 and F1 Sauger)", "Spawning Aggregation (Sampled for F0 and F1 Sauger)", "Spawning Aggregation (Sampled for F0 and F1 Sauger)", "Juvenile Nursery Habitat (Sampled for F2 Sauger)", "Juvenile Nursery Habitat (Sampled for F2 Sauger)", "Juvenile Nursery Habitat (Sampled for F2 Sauger)")
  )
points_sf <- st_as_sf(
  points_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lat/long
)
points_sf$label <- factor(
  points_sf$label,
  levels = sort(unique(points_sf$label), decreasing = TRUE)
)

map <- ggplot() +
       layer_spatial(tiles) +
       geom_sf(data = roi, color = "blue", linewidth = 1, lineend = "round") +
       geom_sf(data = boysen, fill = "blue", color = "blue") +
       geom_sf(data = points_sf, aes(color = label), size = 4) +
       scale_color_manual(values = c("Spawning Aggregation (Sampled for F0 and F1 Sauger)" = "#440154FF",
                                     "Juvenile Nursery Habitat (Sampled for F2 Sauger)" = "#F2844BFF"))+
      theme(legend.position = c(0.02, 0.98),
            legend.justification = c(0, 1),
            legend.position.inside = NULL)+
      coord_sf(xlim = c(-109.25, -108.00),
               ylim = c(42.55, 43.6),
               expand = FALSE,
               datum = sf::st_crs(4326)) +
      labs(x = "Longitude",
           y = "Latitude",
           color = NULL)+
      annotation_scale(location = "br")+
      annotation_north_arrow(location = "br",
                             which_north = "true",
                             style = north_arrow_fancy_orienteering(),
                             pad_y = unit(0.075, "npc"))
map
setwd("/Users/samjohnson/Desktop/")
ggsave("map.png", dpi = 600)
save.image(file = "WRDmap.RData")
# color key


