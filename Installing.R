if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MetaboAnalystR")

library(gage)
library(pathview)
library(gageData)

# Vector of logFC values with KEGG compound IDs
kegg_vector <- setNames(test_extract$logFC, test_extract$kegg_map)

# Run pathway enrichment
data(kegg.sets.hs)  # human KEGG pathways
gage_res <- gage(kegg_vector, gsets = kegg.sets.hs, ref = NULL, samp = NULL)

# View enriched pathways
head(gage_res$greater)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gage")
BiocManager::install("gageData")



# Install
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MetaboAnalystR")

library(MetaboAnalystR)

# Initialize
mSet <- InitDataObjects("conc", "pathora", FALSE)
mSet <- Setup.MapData(mSet, dataSet = test_extract$HMDB)
mSet <- CrossReferencing(mSet, "name")
mSet <- CreateMappingResultTable(mSet)

# Do pathway analysis
mSet <- SetMetabolomeFilter(mSet, "hsa")  # hsa = human
mSet <- PerformPathwayAnalysis(mSet)

# View results
head(mSet$analSet$path.mat)


metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

# Create a personal library folder
dir.create("~/Rlibs", showWarnings = FALSE)

# Tell R to use it
.libPaths(c("~/Rlibs", .libPaths()))

cat('.libPaths(c("~/Rlibs", .libPaths()))\n', file = "~/.Rprofile", append = TRUE)


source("https://bioconductor.org/biocLite.R")
biocLite("SSPA")

install.packages("~/Downloads/SSPA_1.14.2.tar.gz", repos = NULL, type = "source")

install.packages("pacman")

library(pacman)
pacman::p_load(
  "impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore",
  "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes", "BiocParallel",
  "MSnbase", "multtest", "RBGL", "edgeR", "fgsea"
)

# This is what I'm doing to pass through gfortran problem
# sudo mkdir -p /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0
# 
# sudo ln -s /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0/libgfortran.dylib \
# /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0/libgfortran.dylib
# 
# sudo ln -s /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0/libquadmath.dylib \
# /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0/libquadmath.dylib
# 
# sudo ln -s /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0/libemutls_w.a \
# /opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0/libemutls_w.a



devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)

library(MetaboAnalystR)

