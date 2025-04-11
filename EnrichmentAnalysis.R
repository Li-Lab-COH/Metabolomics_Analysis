
library(MetaboAnalystR)
library(dplyr)

# Step 1: Create an empty MetaboAnalyst object
mSet <- InitDataObjects("conc", "msetora", FALSE)

# Step 2: Use your significant compound names (assume they are in `Compounds`)
sig_cmpds <- unique(test_extract$Compounds)

# Optional: remove NAs
sig_cmpds <- sig_cmpds[!is.na(sig_cmpds)]

# Step 3: Set compound names for enrichment
# Input type can be: "hmdb", "kegg", or "name"
mSet <- SetMetabolomeFilter(mSet, "hsa")  # "hsa" for human pathway
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway")  # Options: "smpdb_pathway", "kegg_pathway", etc.
mSet <- Setup.MapData(mSet, sig_cmpds)

# Step 4: Perform enrichment analysis
mSet <- CrossReferencing(mSet, "name")  # use "name" if you're passing compound names
mSet <- CreateMappingResultTable(mSet)
mSet <- CalculateOraScore(mSet)

# Step 5: View or extract results
enrich_results <- mSet$analSet$ora.mat
head(enrich_results)

# Optional: Plot results
PlotEnrichDotPlot(mSet, "pathway_enrich_plot.png", format = "png", dpi = 300)








# 1. Load the library
library(MetaboAnalystR)

# 2. Initialize data object for ORA (over-representation analysis)
mSet <- InitDataObjects("conc", "msetora", FALSE)

# 3. Set organism (use "hsa" for human)
mSet <- SetMetabolomeFilter(mSet, "hsa")

# 4. Set the metabolite set library
# Options include: "smpdb_pathway", "kegg_pathway", "hmdb_pathway", "metaboanalyst_pathway"
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway")  # you can also try "kegg_pathway"

# 5. Prepare your compound list
sig_cmpds <- unique(test_extract$Compounds)  # Use compound names
sig_cmpds <- sig_cmpds[!is.na(sig_cmpds)]

# 6. Map compounds to known IDs
mSet <- Setup.MapData(mSet, sig_cmpds)
mSet <- CrossReferencing(mSet, "name")  # You could try "hmdb" or "kegg" if more accurate

# 7. Create mapping table (this loads pathway info from MetaboAnalyst)
mSet <- CreateMappingResultTable(mSet)

# 8. This step is crucial and may be missing in some older pipelines:
# Manually assign the pathway lib type to avoid the error
mSet$dataSet$path.lib <- "smpdb_pathway"       # Must match SetCurrentMsetLib()
mSet$pathwaylibtype <- "smpdb"                 # Use "kegg" if using KEGG pathway

# 9. Now calculate enrichment
mSet <- CalculateOraScore(mSet)

# 10. View or save results
head(mSet$analSet$ora.mat)
write.csv(mSet$analSet$ora.mat, "metaboanalyst_enrichment_results.csv", row.names = FALSE)










# Pull the query vector and hit values
query_vec <- mSet$name.map$query.vec
hit_vals <- mSet$name.map$hit.values

# Identify indices where no match was found (i.e., NA in hit.values)
unmatched_idx <- which(is.na(hit_vals))

# Extract the corresponding compound names from the query vector
unmatched_names <- query_vec[unmatched_idx]

# View the result
print(unmatched_names)



all_names <- mSet$dataSet$cmpd
matched_names <- names(mSet$dataSet$map.table)
unmatched_names <- setdiff(all_names, matched_names)




639/9672

dim(test_extract)
length(test_extract$Compounds)
length(test_extract$HMDB[test_extract$HMDB != "-"])
length(test_extract$kegg_map[test_extract$kegg_map != "-"])
sum(test_extract$kegg_map == "-" & test_extract$HMDB == "-")
sum(test_extract$PubChem.CID == "-")


test_extract %>%
  filter(kegg_map == "-", HMDB == "-") %>%
  pull(Compounds)







# Load required libraries
library(dplyr)
library(tidyr)
library(webchem)   # To interact with PubChem (retrieving properties, synonyms)
library(KEGGREST)  # To query KEGG for pathway mappings
library(purrr)

# 1. Preprocess your test_extract data --------------------------------------
# Remove columns that contain individual measurements ("IF" columns)
metadata_cols <- grep("^IF", colnames(test_extract), invert = TRUE, value = TRUE)
metabolite_meta <- test_extract %>% select(all_of(metadata_cols))

# For demonstration, we will work with these key columns:
# - 'Index'            : Unique identifier for each metabolite
# - 'Compounds'        : The vendor compound name
# - 'PubChem.CID'      : PubChem compound ID (as provided)
# - 'HMDB'             : HMDB IDs, if available
# - 'Formula'          : Chemical formula, etc.
metabolites <- metabolite_meta %>%
  select(Index, Compounds, PubChem.CID, HMDB, CAS, Formula)

# Ensure PubChem.CID is a character string
metabolites <- metabolites %>% mutate(PubChem.CID = as.character(PubChem.CID))

# 2. Retrieve Chemical Properties from PubChem -----------------------------
# Here we query PubChem for common properties: Canonical SMILES, InChIKey, and IUPACName.
# This step uses the 'pc_prop' function from webchem.


# Remove metabolites with missing or invalid PubChem CIDs (e.g., "-")
metabolites <- metabolites %>%
  filter(PubChem.CID != "-", !is.na(PubChem.CID)) %>%
  mutate(PubChem.CID = as.character(PubChem.CID))

message("Retrieving PubChem properties...")
pc_props <- pc_prop(metabolites$PubChem.CID, 
                    properties = c("CanonicalSMILES", "InChIKey", "IUPACName"),
                    verbose = TRUE)

pc_props2 <- pc_prop(metabolites$PubChem.CID,
                    verbose = TRUE)

# Merge the retrieved properties with our metabolite table;
# note that pc_props$CID is the PubChem CID from PubChem, which should match our PubChem.CID column.

# Ensure both are character type before joining
pc_props <- pc_props %>% mutate(CID = as.character(CID))
metabolites <- left_join(metabolites, pc_props, by = c("PubChem.CID" = "CID"))



# 3. Extract KEGG IDs from PubChem Synonyms ----------------------------------
# Many PubChem entries include cross-reference synonyms like "KEGG:C00031".
# Use the pc_synonyms() function to fetch synonyms for each PubChem CID.
message("Retrieving PubChem synonyms to extract KEGG IDs...")
synonyms_list <- pc_synonyms(metabolites$PubChem.CID)
pc_synonyms(131833064)
# Define a helper function to extract a KEGG identifier from a vector of synonyms
extract_kegg <- function(syns) {
  if (is.null(syns) || length(syns) == 0) {
    return(NA_character_)
  }
  # Look for any synonym containing a pattern like "KEGG"
  # Some synonyms might appear as "KEGG:C00031" or similar.
  kegg_syn <- grep("KEGG[:_]", syns, value = TRUE)
  if (length(kegg_syn) == 0) {
    return(NA_character_)
  } else {
    # Extract the KEGG compound ID (e.g., "C00031") using a regex
    kegg_id <- sub(".*(C\\d+).*", "\\1", kegg_syn[1])
    return(kegg_id)
  }
}

# Apply the extraction function for each metabolite
metabolites$KEGG <- map_chr(synonyms_list, extract_kegg)

# 4. Map KEGG IDs to Pathways -----------------------------------------------
# For those metabolites where a KEGG ID was retrieved, use KEGGREST to fetch pathway associations.
# We query KEGG with the endpoint keggLink("pathway", "cpd:{KEGG}").

# Only work with rows that have a non-missing KEGG ID
kegg_metabs <- metabolites %>% filter(!is.na(KEGG))

# Define a function to retrieve pathways for a given KEGG compound ID.
get_kegg_map <- function(kegg_id) {
  tryCatch({
    # Prepend the KEGG compound prefix ("cpd:")
    query_id <- paste0("cpd:", kegg_id)
    # keggLink returns a named character vector (or possibly empty if nothing is found)
    paths <- keggLink("pathway", query_id)
    # If pathways are found, remove the "path:" prefix and collapse multiple IDs with a separator.
    if (length(paths) > 0) {
      # Remove the "path:" prefix
      clean_paths <- sub("path:", "", paths)
      # Collapse if there are multiple pathways
      return(paste(unique(clean_paths), collapse = "; "))
    } else {
      return(NA_character_)
    }
  }, error = function(e) {
    return(NA_character_)
  })
}

# Query KEGG for each compound in kegg_metabs
message("Mapping KEGG IDs to pathways...")
kegg_metabs <- kegg_metabs %>% 
  rowwise() %>% 
  mutate(Pathways = get_kegg_map(KEGG)) %>%
  ungroup()

# Merge the pathway information back into the complete metabolites table.
metabolites <- left_join(metabolites, kegg_metabs %>% select(Index, Pathways), by = "Index")

# 5. Save the final annotation table -----------------------------------------
# At this point, the `metabolites` data frame now has these annotation columns:
# - CanonicalSMILES, InChIKey, IUPACName (from PubChem)
# - KEGG (the extracted KEGG compound ID, e.g., "C00031")
# - Pathways (a semicolon-separated string of KEGG pathway IDs)
#
# This annotation table can now be used to build a TERM2GENE mapping for enrichment analysis.
message("Final annotated metabolite table:")
print(head(metabolites))

# (Optional) Write the annotation table to a file for further inspection.
write.csv(metabolites, file = "metabolite_annotations.csv", row.names = FALSE)


#------------------ Kegg ---------------------------


# Create a long-format TERM2GENE table
term2gene <- test_extract %>%
  select(Index, kegg_map) %>%
  filter(kegg_map != "-") %>%
  mutate(kegg_map = strsplit(kegg_map, ",")) %>%
  unnest(kegg_map) %>%
  rename(TERM = kegg_map, GENE = Index)  # "GENE" will be metabolite index

foreground_metabs <- metabolites$Index

# All metabolites that were measured (e.g., from the input to limma)
background_metabs <- unique(metlData_log2$Index)  # or from raw data

library(clusterProfiler)

# Run enrichment (ORA) using your custom TERM2GENE object
enrich_res <- enricher(
  gene          = foreground_metabs,
  universe      = background_metabs,
  TERM2GENE     = term2gene,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1  # adjust as needed
)


hmdb_ids <- test_extract$HMDB[test_extract$HMDB != "-"]

mSet <- InitDataObjects("mset", "pathora", FALSE)
mSet <- SetMetabolomeFilter(mSet, "hsa") # hsa = Homo sapiens
mSet <- Setup.MapData(mSet, hmdb_ids)
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- CreateMappingResultTable(mSet)
mSet <- SetKEGG.PathLib(mSet, "hsa", lib.version = "current")
mSet <- CalculateOraScore(mSet, nodeImp = "rbc", method = "fisher")


# View results
res <- mSet$analSet$ora.mat
head(res)









library(MetaboAnalystR)

# Step 1: Initialize the object
mSet <- InitDataObjects("mset", "pathora", FALSE)
mSet <- SetMetabolomeFilter(mSet, "hsa")  # "hsa" for Homo sapiens

# Step 2: Provide your HMDB IDs
hmdb_ids <- test_extract$HMDB
mSet <- Setup.MapData(mSet, hmdb_ids)

# Step 3: Map IDs
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- CreateMappingResultTable(mSet)

# Check mapped compounds
mapped <- mSet$dataSet$map.table
print(head(mapped))

# Step 4: Set KEGG library
mSet <- SetKEGG.PathLib(mSet, "hsa", lib.version = "latest")  # or "3.0" if needed

# Step 5: Perform enrichment
mSet <- CalculateOraScore(mSet, nodeImp = "rbc", method = "fisher")

# Step 6: View results
head(mSet$analSet$ora.mat)

ora_res <- mSet$analSet$ora.mat
print(head(ora_res))


ora_res_df <- as.data.frame(ora_res)
sig_pathways <- ora_res_df[ora_res_df$`Raw p` < 0.05, ]

sig_pathways <- ora_res[ora_res$`Raw p` < 0.05, ]

library(ggplot2)

ggplot(sig_pathways, aes(x = reorder(rownames(sig_pathways), -`Raw p`), y = -log10(`Raw p`))) +
  geom_col() +
  coord_flip() +
  labs(x = "Pathway", y = "-log10(p-value)", title = "Enriched Metabolic Pathways") +
  theme_minimal()

