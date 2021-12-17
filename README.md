# Microbiome-Metabolome-Integration
Integrative analysis of longitudinal microbiome and metabolome data.

## Background and Approach
Untargeted metabolomics data is challenging to interpret since it is typically comprised of thousands of known and unknown metabolites. An additional challenge is in integrating this high-dimensional data with other high-dimensional data (e.g., microbiome, transcriptome). Thus, I applied a data reduction technique, weighted co-network analysis, to group metabolites that are highly correlated with one another into modules. I applied the [WGCNA package in R](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) to metabolomics data (after log-transformation and filtering). After grouping metabolites into modules, an Eigennode is calculated for each module for each biological replicate, which allows for modeling as a continuous variable in statistical models. 

