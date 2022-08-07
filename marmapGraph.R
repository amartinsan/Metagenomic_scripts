
#Code to graph Gulf of Mexico using the NOAA and TOPO database with MARMAP
library(marmap)


#Get coordinates
GOM <- getNOAA.bathy(,lon1 = -98.1, lon2=-90, lat1 = 27 , lat2 = 18.08 ,resolution = 3 , keep = T)
#Color vector
blues<-c("steelblue4","deepskyblue4","deepskyblue3","deepskyblue2","deepskyblue","skyblue","lightblue")
#Plotgraph
plot(GOM,image=TRUE,land=TRUE,lwd=0.4,bpal=list(c(0,max(GOM),c("tan")),c(min(GOM),0,blues)))#blues es el vector con 
# highlight coastline
plot(GOM, lwd = 2, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline

#add sample coordinates to plot in decimal 
sampling <- data.frame(x=c(-95.8369994,-95.24361111,-93.6105555,-94.3475000),y=c(25.88527778,25.892222,19.627222,18.94194444))
#When the graph its to large (higher latitude), longitude coordinates make a space in the borders of the X axis
coord <- par("usr")
atl2 <- getNOAA.bathy(coord[1], coord[2], coord[3], coord[4], res = 3)
plot(atl2,image=TRUE,land=TRUE,lwd=0.4,bpal=list(c(0,max(GOM),c("tan","tan1")),c(min(GOM),0,blues)),drawlabels=T,xlab = )
#Graph the sampling coordinates
points(sampling$x, sampling$y, pch = 25, col = "black",bg = c("red","green","yellow","black"),cex=2.6)
#Define legend
legend("topright", fill=c("red","green","yellow","black"), c("A1-877 m","A3-1926 m","C11-857 m","D18-1281 m"),bty=0,bg="white")

#plot(GOM, lwd = 1, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline


#lot_ordination(D18_BAC_PHYSEQ,ordination = ordiD18,type = "samples",color = "Colorz", shape = "Seccion") +geom_point(size=4.5) + geom_text(label = D18_BAC_PHYSEQ@sam_data$Muestra,hjust = 0.5,vjust=1.5) +theme_bw() + theme(legend.title = element_blank(), axis.title = element_text(size = 18), axis.text = element_text(size = 15),panel.grid = element_blank()) + scale_color_brewer(palette = "Paired") 


ordiA1=ordinate(A1_BAC_physez,"bray",method = "NMDS")
plot_ordination(A1_BAC_physez,ordination = ordiA1,type = "samples", shape = "Seccion",color = "Muestra") +geom_point(size=4.5) + geom_text(label = D18_BAC_PHYSEQ@sam_data$Muestra,hjust = 0.5,vjust=1.5) +theme_bw() + theme(legend.title = element_blank(), axis.title = element_text(size = 18), axis.text = element_text(size = 15),panel.grid = element_blank()) + scale_color_brewer(palette = "Paired") 
