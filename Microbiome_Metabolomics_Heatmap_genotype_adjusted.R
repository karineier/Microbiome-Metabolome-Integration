##### heatmap with genotype-adjusted ASV-module associations #####
setwd("/Users/karineier/Documents/Mecp2/Microbiome/Metabolomics_integration")

library(WGCNA)
library(phyloseq)

load("/Users/karineier/Documents/Mecp2/Data_integration/LMER_models_traits_v_WGCNA_females.RData")

FDRs.ASVs.genotype = read.csv("ASV_Module_adjPVal_gen_adj.csv")
FC.ASVs.genotype = read.csv("ASV_Module_FoldChange_gen_adj.csv")

load("/Users/karineier/Documents/Mecp2/Microbiome/phyloseq_nochim_silva.RData") # loading phyloseq data
ps = ps.silva.nochim

# remove Zymo/mock 
ps = subset_samples(ps, grepl("^JL.*Week", sample_names(ps)))
ps = prune_taxa(taxa_sums(ps)>0, ps)

head(tax_table(ps)) # checking taxa table

m = as(otu_table(ps), "matrix") # getting otu table and making it into a matrix - note: these are really ASVs
m = m + 1 # adding one to prevent log(0) issues
m = t(m) # transforming data so that columns = samples

m = m[,-c(50,53,65,223,242,682)] ## removing samples with very low ASV/read counts

pdata.0 = read.csv("sample_info_WGCNA.csv")
met.samples = read.csv("/Users/karineier/Documents/Mecp2/Metabolomics_Microbiome_Integration/JL787_0819_sample_key_mouse_fecal_16s_matched_with_metabolomic_samples.csv")

pdata.1 = pdata.0[-c(50,53,65,223,242,682),] # removing samples with very low ASV/read counts
rownames(pdata.1) = pdata.1$Sample_ID

pdata.2 = merge(pdata.1, met.samples, by=("SampleID"))
pdata.met.f = subset(pdata.2, Sex=="F")
pdata.met.m = subset(pdata.2, Sex=="M")

met.colnames = as.character(pdata.met.f$SampleID) ### names of samples that we need to select and use from microbiome data

m.met.f = m[,met.colnames]

FDRs.ASVs.genotype.1 = t(FDRs.ASVs.genotype[,c(9:13)])
FC.ASVs.genotype.1 = t(FC.ASVs.genotype[,c(9:13)])

modules = c("Salmon", "Green", "Blue", "Turquoise", "Grey")
rownames(FDRs.ASVs.genotype.1) = modules
rownames(FC.ASVs.genotype.1) = modules

textMatrix.genotype = paste(ifelse((signif(FDRs.ASVs.genotype.1, 1))<0.05, "*", ""), sep="")
dim(textMatrix.genotype) = dim(FDRs.ASVs.genotype.1)

tiff("ASV_Heatmap_Females_Modules_sigforGenotype_gen_adj_Family.tiff", res=400, height=4, width=10, units="in")
map1 = labeledHeatmap(Matrix = FC.ASVs.genotype.1,
                      xLabels = FC.ASVs.genotype$Family,
                      yLabels = rownames(FC.ASVs.genotype.1),
                      ySymbols = rownames(FC.ASVs.genotype.1),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.35,
                      zlim = c(-8,8),
                      cex.lab.y = 1,
                      main=paste("ASVs Associated with Metabolomic Modules: Females"))
dev.off()

tiff("ASV_Heatmap_Females_Modules_sigforGenotype_gen_adj_Genus.tiff", res=400, height=4, width=10, units="in")
map1 = labeledHeatmap(Matrix = FC.ASVs.genotype.1,
                      xLabels = FC.ASVs.genotype$Genus,
                      yLabels = rownames(FC.ASVs.genotype.1),
                      ySymbols = rownames(FC.ASVs.genotype.1),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.35,
                      zlim = c(-8,8),
                      cex.lab.y = 1,
                      main=paste("ASVs Associated with Metabolomic Modules: Females"))
dev.off()

write.table(FC.ASVs.genotype, file="FC_ASVs_associated_w_genotype-assoc_metabolomic_modules.txt", col.names=T, sep="\t")

### Subsetting to only ASVs that are associated with main effect of genotype ###

overlap = read.csv("Genotype_MetabolomeModules_OverlappingASVs_Females_venn_result27564.csv")

colnames(FDRs.ASVs.genotype.1) = FDRs.ASVs.genotype$ASV
colnames(FC.ASVs.genotype.1) = FC.ASVs.genotype$ASV

FDRs.ASVs.genotype.sigonly = FDRs.ASVs.genotype.1[,overlap$ASVs]
FC.ASVs.genotype.sigonly = FC.ASVs.genotype.1[,overlap$ASVs]

overlap$ASV = overlap$ASVs
genotype.sigonly.annot = merge(overlap, FC.ASVs.genotype, by="ASV")

genotype.sigonly.annot = genotype.sigonly.annot[
  with(genotype.sigonly.annot, order(Phylum, Class, Order, Family)),
  ]

textMatrix.genotype = paste(ifelse((signif(FDRs.ASVs.genotype.sigonly, 1))<0.05, "*", ""), sep="")
dim(textMatrix.genotype) = dim(FDRs.ASVs.genotype.sigonly)

tiff("ASV_Heatmap_Females_Modules_sigforGenotype_gen_adj_GensigASVsOnly - Family.tiff", res=400, height=3, width=7, units="in")
map2 = labeledHeatmap(Matrix = FC.ASVs.genotype.sigonly,
                      xLabels = genotype.sigonly.annot$Family,
                      yLabels = rownames(FC.ASVs.genotype.sigonly),
                      ySymbols = rownames(FC.ASVs.genotype.sigonly),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.5,
                      zlim = c(-8,8),
                      cex.lab.y = 1,
                      main=paste("ASVs Associated with Metabolomic Modules: Females"))
dev.off()

tiff("ASV_Heatmap_Females_Modules_sigforGenotype_gen_adj_GensigASVsOnly - Genus.tiff", res=400, height=3, width=7, units="in")
map2 = labeledHeatmap(Matrix = FC.ASVs.genotype.sigonly,
                      xLabels = genotype.sigonly.annot$Genus,
                      yLabels = rownames(FC.ASVs.genotype.sigonly),
                      ySymbols = rownames(FC.ASVs.genotype.sigonly),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.5,
                      zlim = c(-8,8),
                      cex.lab.y = 1,
                      main=paste("ASVs Associated with Metabolomic Modules: Females"))
dev.off()

tiff("ASV_Heatmap_Females_Modules_sigforGenotype_gen_adj_GensigASVsOnly - ASV.tiff", res=400, height=3, width=7, units="in")
map2 = labeledHeatmap(Matrix = FC.ASVs.genotype.sigonly,
                      xLabels = genotype.sigonly.annot$ASV,
                      yLabels = rownames(FC.ASVs.genotype.sigonly),
                      ySymbols = rownames(FC.ASVs.genotype.sigonly),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.5,
                      zlim = c(-8,8),
                      cex.lab.y = 1,
                      main=paste("ASVs Associated with Metabolomic Modules: Females"))
dev.off()

