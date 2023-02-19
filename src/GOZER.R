# Header Info ------------------------------------------------------------------
#
# GOZER
# GSEA and ORA Zero-Entry R script
#
# a //FeAR// R-script - 11-Jan-2023
#
# Based on "Biomedical Knowledge Mining using GOSemSim and clusterProfiler"
# by Guangchuang Yu et al.
# References:
# https://github.com/YuLab-SMU/biomedical-knowledge-mining-book
# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
# http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
#
# Pipeline for the hypergeometric-test-based enrichment analysis of GO-terms and
# KEGG pathways.
# GO and KEGG are the most frequently used DBs for functional analysis. They are
# typically the first choice because of their long-standing curation and
# availability for a wide range of species.
#
# GENE ONTOLOGY
#   The clusterProfiler package implements the enrichGO() function for gene
#   ontology over-representation analysis (ORA).
#   GO analyses (groupGO(), enrichGO() and gseGO()) support any organisms that
#   have an OrgDb object available (e.g., "org.Hs.eg.db",  "org.Mm.eg.db", ...).
#
# KEGG PATHWAY
#   The clusterProfiler package supports all organisms that have KEGG annotation
#   data available in the KEGG database.
#
# Warning:
#   The KEGG FTP service is not freely available for academic use since 2012,
#   and there are many software packages using out-dated KEGG annotation data.
#   The clusterProfiler package supports downloading the latest online version
#   of KEGG data using the KEGG website, which is freely available for academic
#   users. Both the KEGG pathway and module are supported in clusterProfiler.





# Org selection and Pkg loading ------------------------------------------------

# The famous package from CMA-laboratories
library(cmatools)

# The main package for ORA and GSEA in R (by Guangchuang Yu)
# Provides bitr(), enrichGO(), simplify(), enrichKEGG(), browseKEGG(), ...
library(clusterProfiler)

# NOTE:
# If getOption("clusterProfiler.download.method") defaults to "libcurl",
# this may prevent downloading KEGG data...
# You can try setting the "auto" option to use enrichKEGG() function.
getOption("clusterProfiler.download.method")
R.utils::setOption("clusterProfiler.download.method", "auto")
getOption("clusterProfiler.download.method")

# clusterProfiler companion package to plot ORA and GSEA results (by Yu G.)
# Provides pairwise_termsim(), emapplot(), goplot(), cnetplot(), dotplot(), ...
# All the visualization methods are developed based on 'ggplot2' graphics.
library(enrichplot)
#library(ggplot2)

# To download KEGG enrichment plots locally (by Luo and Brouwer 2013)
# Currently only KEGG pathways are implemented, but pathways from Reactome,
# NCI and other databases will (hopefully) be supported in the future.
library(pathview)
library(ggnewscale)

# Select one organism among the following:
  # human
  # mouse
  # rat
  # fly
  # arabidopsis
  # yeast
  # zebrafish
  # worm
  # e.coli (strain K12)
org <- "human"

if (org == "human") {
  org_db <- "org.Hs.eg.db"  # for GO
  org_code <- "hsa"         # for KEGG
} else if (org == "mouse") {
  org_db <- "org.Mm.eg.db"
  org_code <- "mmu"
} else if (org == "rat") {
  org_db <- "org.Rn.eg.db"
  org_code <- "rno"
} else if (org == "fly") {
  org_db <- "org.Dm.eg.db"
  org_code <- "dme"
} else if (org == "arabidopsis") {
  org_db <- "org.At.tair.db"
  org_code <- "ath"
} else if (org == "yeast") {
  org_db <- "org.Sc.sgd.db"
  org_code <- "sce"
} else if (org == "zebrafish") {
  org_db <- "org.Dr.eg.db"
  org_code <- "dre"
} else if (org == "worm") {
  org_db <- "org.Ce.eg.db"
  org_code <- "cel"
} else if (org == "e.coli") {
  org_db <- "org.EcK12.eg.db"
  org_code <- "eco"
}

# Load the organism-specific annotation package
library(org_db, character.only = TRUE)

# Source custom functions
# NOTE: This way of sourcing only works within RStudio !!
script.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(file.path(script.dir, "set_specificity.R", fsep = .Platform$file.sep))

# Desktop local path
local.Desktop <- paste(Sys.getenv("USERPROFILE"), "Desktop",
                       sep = .Platform$file.sep)
# User's R local path
local.R <- paste(Sys.getenv("R_USER"), "R", sep = .Platform$file.sep)

# [OR]
local.R <- paste(Sys.getenv("R_LIBS_USER"), "R", sep = .Platform$file.sep)





# Print dates and versions -----------------------------------------------------
{
  cat("\nclusterProfiler ver.:", toString(packageVersion("clusterProfiler")))
  cat("\nclusterProfiler date:", toString(packageDate("clusterProfiler")))
  cat("\n")
  
  # These packages are usually updated biannually...
  cat("\n", org_db, " ver.: ", toString(packageVersion(org_db)), sep = "")
  cat("\n", org_db, " date: ", toString(packageDate(org_db)), sep = "")
  
  # Search for date format of the following types: 2022-09-26 OR 2011-Mar15
  regex <- "\\d{4}-\\w{2,3}-?\\d{2}"
  # NOTE: the double \ used to escape regex special symbols are because strings
  #       in R also use \ as an escape character (\n, \t ...), so to search for
  #       a digit in a string (\d) you must also escape the \ character!
  
  # Gene Ontology (GO)
  cat("\n  |__ Gene Ontology DB: ")
  # NOTE: stringr package allows to easily extract just the regex match rather
  #       than the whole string as grep() function does in R
  GOann_obj <- sub(".db", "GO", org_db, fixed = TRUE)
  stringr::str_extract(paste(help_as_text(GOann_obj), collapse = " "),
                       regex) |> cat()
  
  # KEGG Pathways
  cat("\n  |__ KEGG Pathways DB: ")
  PATHann_obj <- sub(".db", "PATH", org_db, fixed = TRUE)
  stringr::str_extract(paste(help_as_text(PATHann_obj), collapse = " "),
                       regex) |> cat()
  
  cat("\n\n")
}

# Here is an alternative way to retrieve DB time stamps without need to parse
# objects' help pages, that can be also used as a double check:
#grep("GOSOURCEDATE", capture.output(eval(parse(text = org_db))), value = T)
#grep("KEGGSOURCEDATE", capture.output(eval(parse(text = org_db))), value = T)




  
# Load data (DEA results) ------------------------------------------------------

# Choose one of the following:

# Load the example 'geneList' provided by DOSE package (uses Entrez gene IDs)
# DOSE == Disease Ontology Semantic and Enrichment analysis
data(geneList, package = "DOSE")
DEA_results <- data.frame(Entrez_ID = names(geneList), logFC = geneList)
lms(DEA_results, cols = Inf)
cols_of_interest <- c(Ensembl_ID = NULL,
                      Entrez_ID = 1,
                      GeneSymbol = NULL,
                      logFC = 2,
                      adj_pval = NULL)

# [[[ OR ]]]

# Use the DEG list from my Rheumatoid Arthritis microarray study.
# It should have been already loaded along with `cmatools` package.
# Etanercept (anti-TNFa) vs Methotrexate (MTX) contrast, significant DEGs only
DEA_results <- cmatools::DEGs_stat
lms(DEA_results, cols = Inf)
cols_of_interest <- c(Ensembl_ID = NULL,
                      Entrez_ID = NULL,
                      GeneSymbol = 1,
                      logFC = 3,
                      adj_pval = 5)

# [[[ OR ]]]

# Browse local folders for a DEG list in Tab-Separated-Values format.
# You can even start from a "full DEA" result table rather than using the subset
# of significant DEGs, since the subsetting is implemented later in this script.
myFile <- rstudioapi::selectFile(caption = "Select DEG list",
                                 label = "Select",
                                 path = local.Desktop,
                                 filter = "All Files (*)",
                                 existing = TRUE)
# Substitute ":" with "\t" (docker4seq standard), then read the table
DEA_results <- read.table(text = gsub(":", "\t", readLines(myFile)),
                          sep = "\t", header = TRUE, row.names = NULL)
lms(DEA_results, cols = Inf)
# For each of the following features, enter the proper index or NULL for missing
cols_of_interest <- c(Ensembl_ID = 2,
                      Entrez_ID = NULL,
                      GeneSymbol = 1,
                      logFC = 4,
                      adj_pval = 8)





# Reshape data -----------------------------------------------------------------

# Reorganize data and get rid of the unused columns
DEA_results <- DEA_results[,cols_of_interest]
colnames(DEA_results) <- names(cols_of_interest)
lms(DEA_results, cols = Inf)

# Remove rows with one or more NAs in their numeric columns
# (NAs are only allowed in annotation columns)
DEA_results <- DEA_results[complete.cases(
  DEA_results[,sapply(DEA_results, class) == "numeric"]),]
lms(DEA_results, cols = Inf)

# Extract significant DEGs (possibly separating UPs from DOWNs)
thr <- 0.5 # Threshold on absolute log2 Fold Change
DEG_sign <- "Down" # Choose among: "Up", "Down", "All"

DEGs <- DEA_results[DEA_results$adj_pval < 0.05 & abs(DEA_results$logFC) > thr,]
if (DEG_sign == "Up") {
  DEGs <- DEGs[DEGs$logFC > thr, ]
} else if (DEG_sign == "Down") {
  DEGs <- DEGs[DEGs$logFC < -thr, ]
}
lms(DEGs)

# Add more gene IDs
if (isEmpty(grep("ID", names(cols_of_interest)))) {
  # When only GeneSymbols are available use them to retrieve ENSG and ENTREZ IDs.
  # Possible duplicated symbols among DEGs are not an issue for bitr() since
  # they are automatically removed, however duplicates DEGs will be removed by
  # GOZER before merging the data frames in order to disambiguate accessory
  # information (i.e., logFC, adj_pval)
  
  cat("\nAttaching both Ensembl_IDs and Entrez_IDs...\n\n")
  
  # Biological ID Translator to add new IDs of interest (i.e., Entrez Gene ID)
  # see keytypes(org.Hs.eg.db) for a complete list of the available ID types
  more_ids <- bitr(DEGs$GeneSymbol,
                   fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org_db,
                   drop = FALSE)
  
  # WARNING: the following stats could be slightly inaccurate if one or more NAs
  # were present among gene IDs in the DEG list...
  one2many <- dnues(more_ids$SYMBOL)
  ENSmatch <- sum(is.na(more_ids$ENSEMBL))  # this...
  ENTmatch <- sum(is.na(more_ids$ENTREZID)) # ...and this could be better :-/
  ENSmany2one <- sum(duplicated(na.omit(more_ids$ENSEMBL)))
  ENTmany2one <- sum(duplicated(na.omit(more_ids$ENTREZID)))
  cat("\nTranslation Report",
      "\n------------------",
      "\n  ", one2many[1], " out of ", length(unique(DEGs$GeneSymbol)),
      "\tstarting unique GeneSymbols have been mapped onto ",
      "\n  ", one2many[2],
      " distinct ID pairs (ENSEMBL, ENTREZID), producing a ",
      dim(more_ids)[1], " x ", dim(more_ids)[2],
      "\n  translation matrix having:",
      "\n    ", ENSmatch, "\tNAs among Ensembl IDs",
      "\n    ", ENTmatch, "\tNAs among Entrez Gene IDs",
      "\n    ", ENSmany2one, "\tDuplicated Ensembl IDs",
      "\n    ", ENTmany2one, "\tDuplicated Entrez Gene IDs",
      "\n\n", sep = "")
  
  # Remove duplicates from DEG list, always keeping the most significant one
  DEGs <- DEGs[order(DEGs$adj_pval), ]
  DEGs <- DEGs[!duplicated(DEGs$GeneSymbol), ]
  
  # Left outer join: return all rows from the left table, and any rows with
  # matching keys from the right table
  DEGs2 <- merge(more_ids, DEGs,
                 by.x = "SYMBOL", by.y = "GeneSymbol",
                 all.x = TRUE)
  
} else if (!is.na(cols_of_interest["Ensembl_ID"]) &
           is.na(cols_of_interest["Entrez_ID"])) {
  # When Ensembl_IDs are present use them to retrieve ENTREZ IDs.
  
  cat("\nAttaching Entrez_IDs...\n\n")
  more_ids <- bitr(DEGs$Ensembl_ID,
                   fromType = "ENSEMBL",
                   toType = c("ENTREZID"),
                   OrgDb = org_db,
                   drop = FALSE)
  
  one2many <- dnues(more_ids$ENSEMBL)
  ENTmatch <- sum(is.na(more_ids$ENTREZID))
  ENTmany2one <- sum(duplicated(na.omit(more_ids$ENTREZID)))
  cat("\nTranslation Report",
      "\n------------------",
      "\n  ", one2many[1], " out of ", length(unique(DEGs$GeneSymbol)),
      "\tstarting unique GeneSymbols have been mapped onto ",
      "\n  ", one2many[2],
      " distinct ID pairs (ENSEMBL, ENTREZID), producing a ",
      dim(more_ids)[1], " x ", dim(more_ids)[2],
      "\n  translation matrix having:",
      "\n    ", ENTmatch, "\tEnsembl IDs have no Entrez Gene ID matches",
      "\n    ", ENTmany2one, "\tDuplicated Entrez Gene IDs",
      "\n\n", sep = "")
  
  # Remove duplicates from DEG list, always keeping the most significant one
  DEGs <- DEGs[order(DEGs$adj_pval), ]
  DEGs <- DEGs[!duplicated(DEGs$Ensembl_ID), ]
  
  # Left outer join: return all rows from the left table, and any rows with
  # matching keys from the right table
  DEGs2 <- merge(more_ids, DEGs,
                 by.x = "ENSEMBL", by.y = "Ensembl_ID",
                 all.x = TRUE)
  
} else if (is.na(cols_of_interest["Ensembl_ID"]) &
           !is.na(cols_of_interest["Entrez_ID"])) {
  # When Entrez_IDs are present use them to retrieve Ensembl_IDs.
  
  cat("\nAttaching Ensembl_IDs... TO BE IMPLEMENTED\n\n")
}

lms(DEGs2, cols = Inf)

# Prepare the slim lists (as 'named atomic vectors') to feed to enrich_functions
# (NOTE: data frames wouldn't have worked with enrichplot::cnetplot())
# for GO
ENS_indx <- is.na(DEGs2$ENSEMBL)
DEGs2GO <- as.vector(DEGs2$logFC[!ENS_indx])
names(DEGs2GO) <- DEGs2$ENSEMBL[!ENS_indx]
lms(DEGs2GO, cols = Inf)
# and for KEGG
ENT_indx <- is.na(DEGs2$ENTREZID)
DEGs2KEGG <- as.vector(DEGs2$logFC[!ENT_indx])
names(DEGs2KEGG) <- DEGs2$ENTREZID[!ENT_indx]
lms(DEGs2KEGG, cols = Inf)




  
# GO classification ------------------------------------------------------------
#
#   In groupGO() function, 'gene' input parameter wants a character vector of
#   gene IDs (usually significant DEGs). Any gene ID type supported by the
#   corresponding 'OrgDb' can be directly used in GO analyses. Users just need
#   to specify the 'keyType' parameter. If 'readable' is set to TRUE, the input
#   gene IDs will be converted to gene symbols in the output matrix (see slots
#   @result$geneID and @gene2Symbol).
#
#   groupGO() output 'ggo' is a list of 3 (one per ontology) groupGOResult S4
#   objects, each of them featuring 15 slots. Use View(ego) to explore their
#   structure.

GO <- c("CC", "MF", "BP")

ggo <- list()   # Distribution of DEGs at a specific GO level x
ngo <- vector() # Size of the ontology at level x
neg <- vector() # Non-empty GO-terms at level x
for (ontlg in GO) {
  ggo[[ontlg]] <- groupGO(gene     = names(DEGs2GO),
                          OrgDb    = org_db,
                          keyType  = "ENSEMBL",
                          ont      = ontlg,
                          level    = 2,
                          readable = TRUE)
  lms(ggo[[ontlg]], cols = 4, name = ontlg)
  ngo[ontlg] <- dim(ggo[[ontlg]])[1]
  neg[ontlg] <- sum(ggo[[ontlg]][,3] != 0)
}





# ORA using GO -----------------------------------------------------------------
# 
#   In enrichGO() function, 'gene' input parameter wants a character vector of
#   gene IDs (usually significant DEGs). Any gene ID type supported by the
#   corresponding 'OrgDb' can be directly used in GO analyses. Users just need
#   to specify the 'keyType' parameter. If 'readable' is set to TRUE, the input
#   gene IDs will be converted to gene symbols in the output matrix (see slots
#   @result$geneID and @gene2Symbol).
#
#   enrichGO() output 'ego' is a list of 3 (one per ontology) enrichResult S4
#   objects, each of them featuring 14 slots + a 'root' slot containing only the
#   significant GO terms (kind of living subset of @result slot). Use View(ego)
#   to explore the structure.

GO <- c("CC", "MF", "BP")

ego <- list()   # Enrichment by Gene Ontology
dgo <- vector() # Number of significant GO terms
for (ontlg in GO) {
  ego[[ontlg]] <- enrichGO(gene          = names(DEGs2GO),
                           OrgDb         = org_db,
                           keyType       = 'ENSEMBL',
                           ont           = ontlg,
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
  lms(ego[[ontlg]], cols = 7, name = ontlg)
  dgo[ontlg] <- dim(ego[[ontlg]])[1]
}

# NOTE
# 'pvalueCutoff' is the ADJUSTED p-value cutoff on enrichment tests to report as
# as significant. 'qvalueCutoff' is an additional cutoff for which test must
# pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted
# pvalues and iii) qvalueCutoff on qvalues to be reported as significant.

# Visualize enriched GO terms as directed acyclic graphs (original GO structure).
# showCategory parameter sets the number of GO terms to show, starting from the
# most significant one in the list. NOTE: since some of the most significant
# terms could have in turn one or more significant parents, the actual number of
# colored nodes shown in the plot is generally larger than 'showCategory'. 
goplot(ego$CC, showCategory = 4)
goplot(ego$MF, showCategory = 5, geom = "label")
goplot(ego$BP, showCategory = 10)





# Redundancy reduction --------------------------------------------------------- 
#
#   Get simplified (non-redundant) Enrichment Maps (clusters of similar terms).
#   pairwise_termsim() computes the similarity matrix among GO terms and sets it
#   to the @termsim (and @method) slot of an enrichResult object. Then,
#   simplify() uses this information to prune redundant terms. The function
#   internally calls the GOSemSim package (Yu et al. 2010) to calculate semantic
#   similarity among GO terms and remove those highly similar terms by keeping
#   one representative term.
#   The criteria for the selection of the representative term can be:
#    - 'pval': the most significant term as in REVIGO (built-in function)
#    - the most general/specific term (extended feature)
#   The so obtained "Simplified ontologies" are appended to the 'ego' list.
#   Finally, emapplot() is used to draw both (before- and after-pruning)
#   Enrichment Maps for comparison.

# Pruning criteria; choose among: "pval", "general", "specific"
criterium <- "general"
cutoff <- 0.7 # Prune terms with a similarity above this cutoff

# Average Level Report (ALR) matrix
ALR <- matrix(data = NA, nrow = 3, ncol = 6)
rownames(ALR) <- GO
colnames(ALR) <- c("Pre_AverageLevel",
                   "Post_AverageLevel",
                   "Pre_Mean(Generality)",
                   "Post_Mean(Generality)",
                   "Pre_Mean(Specificity)",
                   "Post_Mean(Specificity)")

for (ontlg in GO) {
  ego[[ontlg]] <- pairwise_termsim(ego[[ontlg]])
  ego[[ontlg]] <- set_specificity(ego[[ontlg]], ont_list.dir = script.dir)
  ontlg_s <- paste0(ontlg, "_s")
  
  if (criterium == "pval") {
    # If GO terms have high similarity select the more significant one (REVIGO)
    ego[[ontlg_s]] <- simplify(ego[[ontlg]],
                               cutoff     = cutoff,
                               by         = "p.adjust",
                               select_fun = min)
  } else if (criterium == "general") {
    # If GO terms have high similarity select the more general one
    ego[[ontlg_s]] <- simplify(ego[[ontlg]],
                               cutoff     = cutoff,
                               by         = "first_GO_level",
                               select_fun = min)
  } else if (criterium == "specific") {
    # If GO terms have a high similarity select the more specific one
    ego[[ontlg_s]] <- simplify(ego[[ontlg]],
                               cutoff     = cutoff,
                               by         = "last_GO_level",
                               select_fun = max)
  }
  dgo[ontlg_s] <- dim(ego[[ontlg_s]])[1] # Update size vector
  
  # Average level of the significant terms (as an index of specificity)
  ALR[ontlg, 1] <- mean(sapply(ego[[ontlg]]$GO_levels, mean))
  ALR[ontlg, 2] <- mean(sapply(ego[[ontlg_s]]$GO_levels, mean))
  ALR[ontlg, 3] <- mean(ego[[ontlg]]$first_GO_level)
  ALR[ontlg, 4] <- mean(ego[[ontlg_s]]$first_GO_level)
  ALR[ontlg, 5] <- mean(ego[[ontlg]]$last_GO_level)
  ALR[ontlg, 6] <- mean(ego[[ontlg_s]]$last_GO_level)
  
  p1 <- emapplot(ego[[ontlg]], cex_label_category = 0.8)
  p2 <- emapplot(ego[[ontlg_s]], cex_label_category = 0.8)
  print(cowplot::plot_grid(p1, p2, ncol = 2,
                     labels = c("Full Map",
                                paste0("Simplified by \'",
                                       criterium, "\' criterium")),
                     label_size = 12, vjust = 5) +
          cowplot::draw_label(paste(ontlg, "Enrichment Map"),
                              fontface = 'bold', size = 20,
                              x = 0.5, y = 1, hjust = 0.5, vjust = 1))
}

{
cat("\nAverage Level Report\nBefore and After simplification\nby \"",
    criterium, "\" criterium\n\n", sep = "")
print(ALR)
cat("\n")
}





# ORA using KEGG ---------------------------------------------------------------
#
#   The 'gene' parameter is a character vector of gene IDs. Input ID type can be
#   only 'kegg', 'ncbi-geneid', 'ncbi-proteinid' or 'uniprot'. Unlike enrichGO()
#   there is no readable 'parameter' for enrichKEGG(). However, users can use
#   the setReadable() function if there is an OrgDb available for the species.
#
#   NOTE
#   The 'kegg' is the primary ID used in KEGG database. The data source of KEGG
#   was from NCBI. A rule of thumb for the 'kegg' ID is Entrez Gene ID for
#   eukaryote species and Locus ID for prokaryotes.

kk <- enrichKEGG(gene         = names(DEGs2KEGG),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 use_internal_data = FALSE) # To use latest online KEGG data !!
lms(kk, cols = 8, name = "KEGG")
dkk <- dim(kk)[1]
kkx <- setReadable(kk, org_db, "ENTREZID") # Convert gene IDs to Symbols
lms(kkx, cols = 8, name = "KEGG")
dkkx <- dim(kkx)[1]

# KEGG Module is a collection of manually defined function units. In some
# situation, KEGG Modules have a more straightforward interpretation.
mkk <- enrichMKEGG(gene         = names(DEGs2KEGG),
                   organism     = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
lms(mkk, cols = 7, name = "KEGG Module")
dmkk <- dim(mkk)[1]
mkkx <- setReadable(mkk, org_db, "ENTREZID") # Convert gene IDs to Symbols
lms(mkkx, cols = 7, name = "KEGG Module")
dmkkx <- dim(mkkx)[1]

# Visualize KEGG pathways with DEGs on KEGG website through browseKEGG()...
{
inkegg <- menu(paste0(kk[,1], " - ", kk[,2]),
               title = "Choose a pathway!", graphics = TRUE)
browseKEGG(kk, kk[inkegg,1])
}

# ...or locally by pathview()
# Here, 'gene.data' parameter can be, as usual, a character vector of gene IDs,
# but if logFCs are also supplied (i.e., a numeric vector with gene IDs as
# names), DEGs in the pathway will be highlighted according to the colorscale
# specified through 'limit', 'bins', 'low', 'mid', 'high', 'na.col' parameters. 
# In case of multiple-gene nodes the logFC average value will be used.
# NOTE: KEGG pathways, in general, includes genes, enzymes, and compounds.
# pathview() accepts compound data through the 'cpd.data' parameter and, if both
# 'cpd.data' and 'gene.data' are provided, genes and compounds are highlighted
# together in the same pathway, using different colorscales.
{
#	Directory to store KEGG pathway data files (.xml) and image files (.png) of
# the 'clean' (i.e., original) KEGG pathways
kegg.folder <- paste(local.R, "KEGG", sep = .Platform$file.sep)
if (!file.exists(kegg.folder)) {
  dir.create(kegg.folder)
  cat("\nNew folder '", basename(kegg.folder), "' has been created inside:\n",
      dirname(kegg.folder), "\n\n", sep = "")
}
# KEGG pathways with gene contributions will be saved in the current WD
setwd(local.Desktop)
# Let's plot now...
inkegg <- menu(c("ALL", paste0(kk[,1], " - ", kk[,2])),
               title = "Choose a pathway!", graphics = TRUE) - 1
if (inkegg > 0) {
  pathway.id <- kk[inkegg,1]
} else if (inkegg == 0) {
  pathway.id <- kk[,1]
} else {
  stop("Aborted")
}
path.desc <- pathview(gene.data   = DEGs2KEGG,
                      pathway.id  = pathway.id,
                      species     = org_code,
                      kegg.native = TRUE, # Try FALSE to get a vectorial output
                      kegg.dir    = kegg.folder,
                      limit       = list(gene=max(abs(DEGs2KEGG)), cpd=1),
                      bins        = list(gene = 20, cpd=10),
                      low         = list(gene = "green", cpd = "blue"),
                      mid         = list(gene = "gray", cpd = "gray"),
                      high        = list(gene = "red", cpd = "yellow"),
                      na.col      = "transparent")
}





# Visualization through enrichplot ---------------------------------------------

# The enrichplot package implements several methods to visualize enriched terms.
# Most of them are general methods that can be used on GO, KEGG, MSigDb, and
# other gene set annotations. Specifically, the enrichplot package supports
# visualizing enrichment results obtained from DOSE, clusterProfiler, ReactomePA,
# and meshes. Both over representation analysis (ORA) and gene set enrichment
# analysis (GSEA) are supported.

# Add title:      + ggtitle("tile_tile_tile_")
# Remove legend:  + theme(legend.position = 'none')

# Maximum number of top-categories (terms) to show
cat_num <- 30
# If users are interested to show some specific pathways (e.g., excluding some
# unimportant pathways among the top categories), users can pass a vector of
# selected pathways to the 'showCategory' parameter (using ID or Description).
# Another solution is using the filter verb to extract a subset of the result as
# described in Chapter 16.
# cat_num <- ego$Description[1:10]



# 1. Bar Plot
#   Bar plot is the most widely used method to visualize enriched terms. It
#   depicts the enrichment scores (e.g., p-values) and gene count or ratio as
#   bar height and color
  barplot(ego$CC, showCategory = cat_num, x = "Count")
  barplot(ego$MF, showCategory = cat_num)
  barplot(ego$BP, showCategory = cat_num)
  barplot(ego$CC, showCategory = cat_num)
  barplot(ego$MF, showCategory = cat_num)
  barplot(ego$BP, showCategory = cat_num)
  barplot(kk, showCategory = cat_num)
  barplot(mkk, showCategory = cat_num)
  
  # You can use other variables derived through mutate() as bar height or color
  # (color redefinition currently not working...)
  mutate(ego$CC, qscore = -log(p.adjust, base = 10)) |>
    barplot(x = "qscore", showCategory = cat_num)
  
  # The Fold Enrichment for a given term (or gene set) is defined as the ratio
  # between the relative amount of DEGs in the given gene set (k/n) and its
  # expected value under the null hypothesis, i.e., the relative size of the
  # gene set in the universe (or background) (K/N).
  # These two ratios are already available in the slot @results$GneRatio and
  # @results$BgRatio, however they are in string form... In order to get the
  # numerical values we need to parse and evaluate expressions
  
  ego$CC@result$FE <-
    sapply(
      parse(
        text = paste0(ego$CC@result$GeneRatio,
                      "/(", ego$CC@result$BgRatio, ")")),
      eval)
  barplot(ego$CC, showCategory = cat_num, x = "FE")
  
  # non vengono stampate le xlabels...
  # Controllare le dimensioni degli slot ...Perch??
  # @geneSets sono > delle righe di result ?



# 2. Dot Plot
  dotplot(ego$CC, showCategory = cat_num)
  dotplot(ego$MF, showCategory = cat_num)
  dotplot(ego$BP, showCategory = cat_num)
  dotplot(kk, showCategory = cat_num)
  dotplot(mkk, showCategory = cat_num)



# 3. Gene-Concept Network
  # NOTE: 'categorySize = pvalue' parameter is  (still?) currently not
  # implemented... so circle size has no meaning
  cnetplot(ego$CC, foldChange = DEGs2GO)
  cnetplot(ego$MF, foldChange = DEGs2GO, circular = TRUE, colorEdge = TRUE)
  cnetplot(ego$BP, foldChange = DEGs2GO)
  cnetplot(kkx, foldChange = DEGs2KEGG)
  cnetplot(mkkx, foldChange = DEGs2KEGG)
  
  
  
# 4. Heatplot
  heatplot(ego$CC, foldChange = DEGs2GO, showCategory = cat_num)
  heatplot(ego$CC_s, foldChange = DEGs2GO, showCategory = cat_num)
  heatplot(ego$MF, showCategory = cat_num)
  heatplot(kkx, showCategory = cat_num)
  heatplot(mkkx, showCategory = cat_num)


    
# 5. Tree plot
  treeplot(ego$CC, nCluster = 5)
  treeplot(ego$MF, nCluster = 5)
  treeplot(ego$BP, nCluster = 7)
  treeplot(ego$BP, nCluster = 7, hclust_method = "average")


  
# 6. Enrichment Map
#     Enrichment map organizes enriched terms into a network with edges
#     connecting overlapping gene sets. Mutually overlapping gene sets tend to
#     cluster together, making it easier to identify functional modules. When
#     the similarity between terms meets a certain threshold (default is 0.2,
#     adjusted by parameter 'min_edge'), there will be edges between terms. The
#     stronger the similarity, the shorter and thicker the edges. The similarity
#     between terms is obtained by function 'pairwise_termsim'.
  
  emapplot(ego$CC)
  emapplot(ego$MF)
  emapplot(ego$BP, layout = "kk")

  

  
  