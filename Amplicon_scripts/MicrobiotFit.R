
library("phyloseq")
#library("dada2")
library("microbiome")
#library("mia")
library("biomformat")
library("tidyverse")
library("xlsx")

#load data
load("C:/Users/AMS/Downloads/2023-02-27-allsamples.rda")

#subset samples of interest microbiofit

#Uriel Alavarado
URI1=subset_samples(physeq = physeq,orden == "807")
URI2=subset_samples(physeq = physeq,orden == "910")
URIEL=merge_phyloseq(URI1,URI2)
#Claudia Godinez
CLAU1=subset_samples(physeq = physeq,orden == "802")
CLAU2=sua Pacheco
YES1=subset_samples(physeq = physeq,orden == "909")
CLAUDIA=merge_phyloseq(CLAU1,CLAU2)
#Yessicbset_samples(physeq = physeq,orden == "806")
YES2=subset_samples(physeq = physeq,orden == "907")
YESSICA=merge_phyloseq(YES1,YES2)
#Regina Garcia
REG1=subset_samples(physeq = physeq,orden == "823")
REG2=subset_samples(physeq = physeq,orden == "906")
REGINA=merge_phyloseq(REG1,REG2)

#Pablo Betrab
PAB1=subset_samples(physeq = physeq,orden == "805")
PAB2=subset_samples(physeq = physeq,orden == "898")
PABLO=merge_phyloseq(PAB1,PAB2)
#Montserrat Plasencia
MONT1=subset_samples(physeq = physeq,orden == "803")
MONT2=subset_samples(physeq = physeq,orden == "908")
MONTSE=merge_phyloseq(MONT1,MONT2)
#No second samples 
#UNA=subset_samples(physeq = physeq,orden == "833")
#DOS=subset_samples(physeq = physeq,orden == "816")
#TRES=subset_samples(physeq = physeq,orden == "809")
#CUATRO=subset_samples(physeq = physeq,orden == "818")

###ALL SAMPLES####
#microbiofit=merge_phyloseq(URI1,URI2,CLAU1,CLAU2,YES1,YES2,REG1,REG2,PAB1,PAB2,MONT1,MONT2,UNA,DOS,TRES,CUATRO) 
#Only samples with second analysis
#Merge samples 
microbiofit=merge_phyloseq(URI1,URI2,CLAU1,CLAU2,YES1,YES2,REG1,REG2,PAB1,PAB2,MONT1,MONT2) 
#Change sample names into order 
sample_names(microbiofit)=microbiofit@sam_data$orden
#Prueba it´s just for ordering de samples in the graph add it to the metadata 
prueba=c("primera","segunda","primera","segunda","primera","segunda","primera","segunda","primera","segunda","primera","segunda")
microbiofit@sam_data$prueba=prueba
metadatMicrofit=data.frame(sample_data(microbiofit))
#Vector for ordering samples along the X axis
mtype=as.character(metadatMicrofit$orden[order(paste(metadatMicrofit$nombre,metadatMicrofit$orden))]) 

#Uriel Alavarado
URI1=subset_samples(physeq = microbiofit,orden == "807")
URI2=subset_samples(physeq = microbiofit,orden == "910")
URIEL=merge_phyloseq(URI1,URI2)
#Claudia Godinez
CLAU1=subset_samples(physeq = microbiofit,orden == "802")
CLAU2=subset_samples(physeq = microbiofit,orden == "909")
CLAUDIA=merge_phyloseq(CLAU1,CLAU2)
#Yessica Pacheco
YES1=subset_samples(physeq = microbiofit,orden == "806")
YES2=subset_samples(physeq = microbiofit,orden == "907")
YESSICA=merge_phyloseq(YES1,YES2)
#Regina Garcia
REG1=subset_samples(physeq = microbiofit,orden == "823")
REG2=subset_samples(physeq = microbiofit,orden == "906")
REGINA=merge_phyloseq(REG1,REG2)

#Pablo Betrab
PAB1=subset_samples(physeq = microbiofit,orden == "805")
PAB2=subset_samples(physeq = microbiofit,orden == "898")
PABLO=merge_phyloseq(PAB1,PAB2)
#Montserrat Plasencia
MONT1=subset_samples(physeq = microbiofit,orden == "803")
MONT2=subset_samples(physeq = microbiofit,orden == "908")
MONTSE=merge_phyloseq(MONT1,MONT2)


#Barplot chunk
#The idea es to aggregat or glom taxa in a taxonomic level, melt the table and define low abundance taxa
y1 <- tax_glom(microbiofit, taxrank = 'Class') # agglomerate taxa
#(y2 = merge_samples(y1, "Partition")) # merge samples on sample variable of interest
#y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
#y4 <- psmelt(y3) # create dataframe from phyloseq object
y4 <- psmelt(y1) # create dataframe from phyloseq object
y4$Class[y4$Abundance < 0.04] <- "<4% abund." 

#Vector of random colors (NOT SO RANDOM)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(34)
color_b=c("navyblue","brown","purple","orange","darkgreen","red","gold2","pink","firebrick",sample(color, 10))

#ggplot(data=y4, aes(x=Sample, y=Abundance, fill=Class)) +  geom_bar(aes(), stat="identity", position="stack") 

#####Since Family level, there are many "unclassified" or "unknown" taxa, do not graph that
#niveles=c("Phylum","Class","Order","Family","Genus")

Muestra=CLAUDIA
niveles=c("Phylum","Class","Order", "Family")


#for cycle to plot every level
for (i in niveles) {
  y1 <- microbiome::aggregate_taxa(Muestra,level = i)
  y4 <- psmelt(y1)
  y4 <- y4 %>%
    mutate(Taxa = if_else(Abundance < 0.025, "<2.5% abund.", !!sym(i)))
  mtype=as.character(Muestra@sam_data$orden[order(paste(Muestra@sam_data$nombre,Muestra@sam_data$orden))]) 
  
  assign(paste("p", i, sep = ""), ggplot(y4, aes(x=Sample, y=Abundance, fill=Taxa)) + 
           geom_bar(stat="identity", position="stack") + scale_x_discrete(limits=mtype) +
           scale_fill_manual(values=color_b)+ theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15)))
}

###Taxonomy table###### 
###Save as .xlsx for clients

mab2OtuT <- data.frame(otu_table(microbiofit))
mab2Taxa <- data.frame(tax_table(microbiofit))
mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
#EXTRAS retirar columna de OTUid y Consensus
mab2join$Row.names = NULL
mab2join$Consensus = NULL
write.xlsx(mab2join,sheetName = "Sheet1",file = "Taxonomytable.txt")
#Remember to edit names and check file
##individual clientes########



clientes= c("URIEL","CLAUDIA","YESSICA","REGINA","PABLO","MONTSE","microbiofit")

for (i in clientes) {
  Muestra=get(i)
  mab2OtuT <- data.frame(otu_table(Muestra))
  mab2Taxa <- data.frame(tax_table(Muestra))
  mab2join <- merge(mab2Taxa,mab2OtuT,by = 0,all = T)
  #EXTRAS retirar columna de OTUid y Consensus
  mab2join$Row.names = NULL
  mab2join$Consensus = NULL
  write.xlsx(mab2join,sheetName = "Sheet1",file = paste("TAXO_TABLE",i,".xlsx",sep = ""))
} 

##############LAPPLY FOR OPTIMIZATION
clientes <- c("URIEL", "CLAUDIA", "YESSICA", "REGINA", "PABLO", "MONTSE")
physeq_list <- lapply(clientes, get)
physeq_combined <- merge_phyloseq(physeq_list)

mab2OtuT <- data.frame(otu_table(physeq_combined))
mab2Taxa <- data.frame(tax_table(physeq_combined))
mab2join <- merge(mab2Taxa, mab2OtuT, by = 0, all = TRUE)
#EXTRAS retirar columna de OTUid y Consensus
mab2join$Row.names <- NULL
mab2join$Consensus <- NULL
write.xlsx(mab2join, sheetName = "Sheet1", file = "TAXO_TABLE.xlsx")



#Krona per sample WINDOWS WONT RUN IT only linux with krona installed (check conda)

for (i in sample_names(microbiofit)) {
  KronaPerSAmple = subset_samples(microbiofit,Nombres==i)
  psadd::plot_krona(KronaPerSAmple,output = paste("/opt/Alchemycode/WORK",i,sep = "_"),variable = "Nombres")
}

##### Diversity measurments ###########

#First transform back relativa abundance to raw reads, for any problem with the metrics of alpha and beta diversity
microbiofit_counts = microbiofit

reads <- microbiofit@sam_data$lecturaseq
reads <- reads[match(rownames(tax_table(microbiofit)), names(reads))]
# extract OTU table from phyloseq object and convert to data frame
otu_df <- as.data.frame(otu_table(microbiofit))
# extract read counts from sample data
reads <- microbiofit@sam_data$lecturaseq
# Multiply each row of the matrix with the corresponding number of reads
otu_mat_counts <- round(t(t(otu_mat) * reads))
microbiofit_counts@otu_table = otu_table(otu_mat_counts, taxa_are_rows=T)
####
microbiofit_counts@sam_data$prueba=prueba


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



########################CORE MICROBIOME######################

Primeras=subset_samples(microbiofit_counts,prueba == "primera")
Segundas=subset_samples(microbiofit_counts,prueba == "segunda")


sets=c("Primeras","Segundas") 

for (i in sets) {
  Muestra = get(i)
  ps1.stool.rel <- microbiome::transform(Muestra, "compositional")
  ps1.stool.rel2 <- prune_taxa(taxa_sums(ps1.stool.rel) > 0, ps1.stool.rel)
  taxonomy <- as.data.frame(tax_table(ps1.stool.rel2))
  core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)
  core.abundance <- sample_sums(core(ps1.stool.rel2, detection = 0.001, prevalence = 50/100))
  prevalences <- seq(.05, 1, .05)
  detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
  ps1.stool.rel2.f <- microbiomeutilities::format_to_besthit(ps1.stool.rel2)
  prevalences <- seq(.05, 1, .05)
  detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
  
  assign(paste("Core",i,sep = "_"), p.core <- plot_core(ps1.stool.rel2.f, 
                                                        plot.type = "heatmap", 
                                                        colours = rev(brewer.pal(5, "Spectral")),
                                                        prevalences = prevalences,
                                                        detections = detections,
                                                        min.prevalence = .5) + 
           xlab("Detection Threshold (Relative Abundance (%))")  +
           ggtitle(paste("Core",i," ")) +
           theme_bw() + 
           theme(axis.text.y = element_text(face="italic"),axis.text.x = element_text(size = 12,angle = 45))
  )
}



p=pPhylum + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15), legend.text = element_text(size = 12), axis.title.y = element_text(size = 15)) + ggtitle("Phylum")
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "PHYLUM.html")
p=pClass + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15), legend.text = element_text(size = 12), axis.title.y = element_text(size = 15)) + ggtitle("Class")
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "CLASS.html")
p=pOrder + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15), legend.text = element_text(size = 12), axis.title.y = element_text(size = 15)) + ggtitle("Order")
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "ORDER.html")
p=pFamily + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15), legend.text = element_text(size = 12), axis.title.y = element_text(size = 15)) + ggtitle("Family")
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "FAMILY.html")
p=pGenus + theme_bw() + theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text = element_text(size = 15), legend.text = element_text(size = 12), axis.title.y = element_text(size = 15)) + ggtitle("Genus")
ax=ggplotly(p)
htmlwidgets::saveWidget(as_widget(ax), "GENUS.html")


ordi_NMDS <- ordinate(Phyloseqnorm,"NMDS","bray",)
ordi_PCoA <- ordinate(Phyloseqnorm,"PCoA","bray")
ordi_PCA <- ordinate(Phyloseqnorm,method = "RDA")

ordi_NMDS <- ordinate(physeq,"NMDS","bray",)


plot_ordination(Phyloseqnorm, ordi_NMDS,color = "orden",shape = "prueba")  +  geom_point(size=2)  + ggtitle("NMDS") + theme_bw()+theme(legend.title = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size = 14))
