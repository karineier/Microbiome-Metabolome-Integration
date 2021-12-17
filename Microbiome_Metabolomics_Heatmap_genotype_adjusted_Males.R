##### heatmap with genotype-adjusted ASV-module associations #####
setwd("/Users/karineier/Documents/Mecp2/Microbiome/Metabolomics_integration")

library(WGCNA)

FDRs.ASVs.genotype = read.csv("ASV_Module_adjPVal_gen_adj_Males.csv")
FC.ASVs.genotype = read.csv("ASV_Module_FoldChange_gen_adj_Males.csv")

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

met.colnames = as.character(pdata.met.m$SampleID) ### names of samples that we need to select and use from microbiome data

m.met.m = m[,met.colnames]

taxonomy = tax_table(ps, errorIfNULL=FALSE)
if( !is.null(taxonomy)) {
  taxonomy = data.frame(as(taxonomy, "matrix"))
}

taxonomy$ASV = rownames(taxonomy)

FC.ASVs.genotype.annot = t(FC.ASVs.genotype)
FC.ASVs.genotype.annot = FC.ASVs.genotype.annot[-1,]
FC.ASVs.genotype.annot = as.data.frame(FC.ASVs.genotype.annot)
FC.ASVs.genotype.annot$ASV = rownames(FC.ASVs.genotype.annot)
FC.ASVs.genotype.annot = merge(FC.ASVs.genotype.annot, taxonomy, by="ASV")

FC.ASVs.genotype.annot = FC.ASVs.genotype.annot[
  with(FC.ASVs.genotype.annot, order(Phylum, Class, Order, Family)),
  ]

rownames(FC.ASVs.genotype) = FC.ASVs.genotype$Module
FC.ASVs.genotype.1 = FC.ASVs.genotype[,-1]

rownames(FDRs.ASVs.genotype) = FDRs.ASVs.genotype$Module
FDRs.ASVs.genotype.1 = FDRs.ASVs.genotype[,-1]

textMatrix.genotype = paste(ifelse((signif(FDRs.ASVs.genotype.1, 1))<0.05, "*", ""), sep="")
dim(textMatrix.genotype) = dim(FDRs.ASVs.genotype.1)

tiff("ASV_Heatmap_Males_Modules_sigforGenotype_gen_adj.tiff", res=400, height=4, width=10, units="in")
map1 = labeledHeatmap(Matrix = FC.ASVs.genotype.1,
                      xLabels = FC.ASVs.genotype.annot$Family,
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
                      main=paste("ASVs Associated with Metabolomic Modules: Males"))
dev.off()

write.table(FC.ASVs.genotype.1, file="FC_ASVs_associated_w_genotype-assoc_metabolomic_modules_Males.txt", col.names=T, sep="\t")

### Subsetting to only ASVs that are associated with main effect of genotype ###

overlap = read.csv("Genotype_MetabolomeModules_OverlappingASVs_Males_venn_result.csv")

overlap$ASV = overlap$ASVs

genotype.sigonly.annot = merge(overlap, FC.ASVs.genotype.annot, by="ASV")

FDRs.ASVs.genotype.2 = as.data.frame(t(FDRs.ASVs.genotype[,-1]))
FDRs.ASVs.genotype.2$ASV = rownames(FDRs.ASVs.genotype.2)

genotype.sigonly.annot.1 = merge(genotype.sigonly.annot, FDRs.ASVs.genotype.2, by="ASV")

genotype.sigonly.annot.2 = genotype.sigonly.annot.1[
  with(genotype.sigonly.annot.1, order(Phylum, Class, Order, Family)),
]

FC.ASVs.genotype.sigonly.heatmapMatrix = matrix(as.numeric(t(genotype.sigonly.annot.2[,c(3:8)])), ncol = ncol(t(genotype.sigonly.annot.2[,c(3:8)])))
FDRs.ASVs.genotype.sigonly = t(genotype.sigonly.annot.2[,c(16:21)])

textMatrix.genotype = paste(ifelse((signif(FDRs.ASVs.genotype.sigonly, 1))<0.05, "*", ""), sep="")
dim(textMatrix.genotype) = dim(FDRs.ASVs.genotype.sigonly)

hubs = read.csv(glue::glue("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/Males/Hub_metabolites.csv"))

moduleorder = data.frame(module = c("tan", "red", "green", "lightgreen", "salmon", "darkorange"))

hubs.0 = plyr::join(moduleorder, hubs)
hubs.0$annotation = replace_na(hubs.0$annotation, "Unassigned Module")


tiff("ASV_Heatmap_Males_Modules_sigforGenotype_gen_adj_GensigASVsOnly - Family.tiff", res=400, height=3, width=9, units="in")
map2 = labeledHeatmap(Matrix = FC.ASVs.genotype.sigonly.heatmapMatrix,
                      xLabels = genotype.sigonly.annot.2$Family,
                      yLabels = hubs.0$annotation,
                      ySymbols = rownames(FC.ASVs.genotype.sigonly.heatmapMatrix),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.5,
                      zlim = c(-8,8),
                      cex.lab.y = 0.49,
                      main=paste("ASVs Associated with Metabolomic Modules: Males"))
dev.off()

tiff("ASV_Heatmap_Males_Modules_sigforGenotype_gen_adj_GensigASVsOnly - ASV.tiff", res=400, height=3, width=9, units="in")
map2 = labeledHeatmap(Matrix = FC.ASVs.genotype.sigonly.heatmapMatrix,
                      xLabels = genotype.sigonly.annot.2$ASV,
                      yLabels = hubs.0$annotation,
                      ySymbols = rownames(FC.ASVs.genotype.sigonly.heatmapMatrix),
                      colorLabels = FALSE,
                      colors = blueWhiteRed(50),
                      textMatrix = textMatrix.genotype,
                      setStdMargins = FALSE,
                      cex.text = 0.65,
                      cex.lab.x = 0.5,
                      zlim = c(-8,8),
                      cex.lab.y = 0.8,
                      main=paste("ASVs Associated with Metabolomic Modules: Males"))
dev.off()



