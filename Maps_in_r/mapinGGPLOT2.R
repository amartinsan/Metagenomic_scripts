if (!require(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))) 
  install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("rnaturalearth")
library("rnaturalearthdata")
library("ggplot2")
library("sf")
#Define the sites or points to plot in latitude and longitude

sites <- data.frame(x=c(-95.8369994,-95.24361111,-93.6105555,-94.3475000),y=c(25.88527778,25.892222,19.627222,18.94194444))
row.names(sites) = c("A1","A3","C11","D18")


"Plot world for example"
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf()


"plot the coordentaes of interest add, legend and text"

ggplot(data = world) +
  geom_sf(show.legend = T, fill="antiquewhite") + 
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 5,  shape = 24, fill = c("blue","darkred","darkgreen","orange")) + 
  geom_sf(data = states, fill = NA) +  
  annotate(geom = "text",x = -94, y = 23, label = "Gulf of Mexico",fontface="italic",color="grey22",size=9)+
  geom_text(data = sites, label=row.names(sites), aes(x = longitude, y = latitude), nudge_y = 0.3, color=c("blue","darkred","darkgreen","orange"),fontface="bold",size=5) +
  annotation_scale(location = "tr", width_hint = 0.5, tick_height = 5) + 
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-98.2, -89.5), ylim = c(18.08, 27.0))+ ggtitle("Sampling sites")+
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "aliceblue"), 
    axis.title = element_blank(),
    axis.text = element_text(size = 22),
    title = element_text(size = 18))
