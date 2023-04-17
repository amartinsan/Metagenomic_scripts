#Normalize counts to CSS, way better for ordinations and not depending in rarefaction
library(metagenomeSeq)
MGS <- phyloseq_to_metagenomeSeq(physeq)
p <- cumNormStatFast(MGS)
MGS <- cumNorm(MGS, p= p)
normmalizado <- MRcounts(MGS, norm = T)
tablaOTUS <- as.data.frame(otu_table(physeq))
#copy phyloseq object
PhyloseqnormCSS <- physeq
#Substitute the otu_table with the normalized one
PhyloseqnormCSS@otu_table <- otu_table(normmalizado, taxa_are_rows = T)
#CHECK STEP
#tablAchek <- as.data.frame(PhyloseqnormCSS@otu_table)
#View(tablAchek)
