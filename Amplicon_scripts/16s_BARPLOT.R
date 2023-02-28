#Collapse low abudnace taxa into others category

#library("psadd") #plot_krona()
library("phyloseq")
#library("dada2")
library("microbiome")
library("microbiomeutilities")
library("biomformat")
library("tidyverse")

#using a phyloseq object COLLPASE LOW ABUNDANCE TAXA INTO ANOTHER CATEGORY

level=tax_glom(phyloseq,taxrank="Phylum") 


y1 <- tax_glom(microbiofit, taxrank = 'Phylum') # agglomerate taxa
#(y2 = merge_samples(y1, "Partition")) # merge samples on sample variable of interest
#y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
#y4 <- psmelt(y3) # create dataframe from phyloseq object
y4 <- psmelt(y1) # create dataframe from phyloseq object
y4$Phylum[y4$Abundance < 0.05] <- "Genera < 5% abund." #rename genera with < 1% abundance
