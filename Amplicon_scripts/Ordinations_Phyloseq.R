library(yingtools2)
library(umap)
library("phyloseq")
library("microbiome")
Muestras=rCLRdata
ordi_NMDS <- ordinate(Muestras,"NMDS","bray",)
p=plot_ordination(Muestras, ordi_NMDS, color = "ordennum",shape = "batch")  +  geom_point(size=2) + scale_color_viridis(option = "turbo") + ggtitle("NMDS") + theme_bw()+theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14))
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "NMDS.html")
#PCoA en vez de NMDS
ordi_PCoA <- ordinate(Muestras,"PCoA","bray")
p=plot_ordination(Muestras, ordi_PCoA, color = "ordennum",shape = "batch")  +  geom_point(size=2.5) + scale_color_viridis() + ggtitle("PCoA") + theme_bw() + theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14))
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "PCoA.html")
#PCA
ordi_PCA <- ordinate(Muestras,method = "RDA")
p=plot_ordination(Muestras, ordi_PCA, color = "ordennum",shape = "batch")  +  geom_point(size=2.5) + scale_color_viridis() + ggtitle("PCA") + theme_bw() + theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14))
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "PCA.html")

#UMAP
s <- get.samp(Muestras)
otu <- get.otu(Muestras,as.matrix=TRUE) %>% t()
umap <- umap(otu)
layou=data.frame(umap$layout) %>% rownames_to_column("sample") %>% rename(UMAP1=X1,UMAP2=X2) %>% left_join(s,by="sample")

p=ggplot(layou,aes(x=UMAP1,y=UMAP2,color=ordennum)) + geom_point(size=2,aes(shape=layou$batch)) + scale_color_viridis(option = "cividis") + theme_bw() + ggtitle("UMAP") + theme_bw() + theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14))

ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "UMAP.html")
