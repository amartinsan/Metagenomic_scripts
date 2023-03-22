
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





