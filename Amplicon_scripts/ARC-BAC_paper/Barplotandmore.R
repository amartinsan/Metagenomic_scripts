#REDO
library("psadd") #plot_krona()
library("phyloseq")
library("dada2")
library("microbiome")
library("biomformat")
library("tidyverse")
#library(microbiomutilies)
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("") & Kingdom != "Unknown")http://127.0.0.1:30243/graphics/plot_zoom_png?width=1920&height=1057
abrel=transform(ps,"compositional")
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(34)
color_b=c("navyblue","firebrick","purple","orange","darkgreen","red","gold2","pink","brown",sample(color, 10))
y1 <- tax_glom(abrel, taxrank = 'Phylum') # agglomerate taxa
y4 <- psmelt(y1)
#y4$Phylum[y4$Abundance < 0.01] <- "Phylum < 1% abund." #rename genera with < 1% abundance
ggplot(data = y4, aes(x = Sample, y = Abundance, fill = reorder(Phylum, -Abundance))) + geom_bar(stat = "identity", position = "stack") +  scale_fill_manual(values = color_b) + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90,size = 19),axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 22),legend.text = element_text(size = 15))
