
#Microbiome Ordinations wiht Phyloseq and more
#!R

#load libraries
library("phyloseq")
library("microbiome")
library("tidyverse")
library("yingtools2")
library("umap")
library("dada2")
library("biomformat")
library("vegan")


##############
muestra=$PHYLOSEQfile
#basic ordinations (NMDS,PCoA,PCA,etc.), normal bray distance for example

ordi_NMDS <- ordinate(physeq = muestra,method = "NMDS","bray") 
ordi_PCoA <- ordinate(physeq = muestra,method = "PCoA","bray")
ordi_PCA <- ordinate(physeq = muestra,method = "RDA","bray")

####Plot with ggplot2##### you can sapply o for cycle to optimize

NMDSPlot = plot_ordination(muestra, ordi_NMDS) + geom_point(size=3) + theme_bw() 
    + theme(legend.title = element_blank())

PCoAplot = plot_ordination(muestra, ordi_PCoA) + geom_point(size=3) + theme_bw() 
    + theme(legend.title = element_blank())

PCAplot = plot_ordination(muestra, ordi_NMDS) + geom_point(size=3) + theme_bw() 
    + theme(legend.title = element_blank())

###UMAP , a little trickier ####
#get samples and OTUs well ASVs 
s <- get.samp(muestra)
otu <- get.otu(muestra,as.matrix=TRUE) %>% t()
umap <- umap(otu)
layou=data.frame(umap$layout) %>% rownames_to_column("sample") %>% rename(UMAP1=X1,UMAP2=X2) %>% left_join(s,by="sample")
UMAPplot = ggplot(layou,aes(x=UMAP1,y=UMAP2,color=ordennum)) + 
    geom_point(size=2,aes(shape=layou$batch)) + 
    theme_bw() + 
    theme(legend.title = element_blank(),)

#####Make a graph "interactive with plotly" #######
#Example

Muestras=rCLRdata
ordi_NMDS <- ordinate(Muestras,"NMDS","bray",)
p=plot_ordination(Muestras, ordi_NMDS, color = "ordennum",shape = "batch")  +  
    geom_point(size=2) + scale_color_viridis(option = "turbo") + 
    ggtitle("NMDS") + theme_bw() + 
    theme(legend.title = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size = 14))

#save plot as gglotly and to save as an .html file, maybe aplly this in shiny       
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "NMDS.html")

#Example on how to cycle a core microbiome graph of a a bunch of phyloseq objects in a vector


#The assing and paste is to save the graphs in a variable with the name of the vector variables

clientes <- c("URIEL","CLAUDIA","YESSICA","REGINA","PABLO","MONTSE","microbiofit")

for (i in clientes) {
  Muestra <- get(i)
  mtype <- as.character(Muestra@sam_data$orden[order(paste(Muestra@sam_data$nombre,Muestra@sam_data$orden))]) 
  div_alpha <- estimate_richness(Muestra, measures = c("Shannon", "Simpson"))
  # extract Shannon and Simpson indices from the alpha diversity table
  div_alpha <- as.data.frame(div_alpha[c("Shannon", "Simpson")])
  # add a column for the sample IDs
  div_alpha$Sample <- rownames(div_alpha)
  # convert data to long format
  div_alpha_long <- tidyr::gather(div_alpha, "measure", "value", -Sample)
  # create barplot using ggplot2
  assign(paste("alpha_", i, sep = ""), ggplot(div_alpha_long, aes(x = Sample, y = value, fill = measure)) +
           geom_bar(stat = "identity", position = "dodge") +
           scale_fill_manual(values = c("skyblue", "lightgreen")) +
           labs(x = "Sample", y = "Alpha diversity") + theme_bw() + scale_x_discrete(limits=mtype))
}



