# These functions uses getGOLevel() to quantify the generality/specificity of GO
# terms. Such an information will be used by simplify() function as criteria to
# reduce redundancy. The main idea is that if two terms have an high similarity,
# the more general/specific term will be retained. This is likely to select the
# parent/child GO terms, respectively.
#
# Modified from the original script suggested by TylerSagendorf (19 Sep 2021)
# https://github.com/YuLab-SMU/clusterProfiler/issues/372
# See RedundancyReduction_OriginalPost.tnotes

# For each ontology, this function gets a list of all GO terms at each level
# from 1:14 (keep in mind that each ID may be present in multiple levels,
# depending on its relationship to other terms in the DAG).
# This creates a list of three lists named 'ont_list' that is saved locally as a
# binary .RData file (namely in the 'user-defined target.dir') to allow a quick
# access by set_specificity() function (see below).
#
# NOTE_1: 'getGOLevel' (used to query GO IDs at a specific level) is not an
#         exported function from 'namespace:clusterProfiler' which is why the
#         triple colon are used to access it.
# NOTE_2: It fails to get terms from level 15 and above for some reason (i.e.,
#         ontology ending -> try to catch the exception...)

ont_list_generate <- function(target.dir = getwd(), depth = 14) {
  
  GO <- c("CC", "MF", "BP")
  ont_list <- lapply(GO, function(onts) {
    lapply(1:depth, function(level) {
      clusterProfiler:::getGOLevel(onts, level)
    })
  })
  names(ont_list) <- GO
  
  save(ont_list, file = paste0(target.dir, "/ont_list.RData"))
}

# Starting from the ontology-specific list 'ont_list' of all the GO terms
# present at each ontology level (terms as a function of levels), this function
# get the "inverse" information (levels as a function of terms) and adds three
# new variables to the @result slot of a given object x of class 'enrichResult':
# - GO_levels: for all GO_terms, the vector of levels to which they take part
# - first_GO_level: for all GO_terms, the lower level they belong to (generality)
# - last_GO_level: for all GO_terms, the higher level they belong to (specificity)
#
# In order to speed up computation it is much better to load the starting data
# structure 'ont_list' (a list of three lists) from somewhere rather than
# computing it from scratch each time the function is invoked! For this reason,
# the file 'ont_list.RData' will be searched in the the 'ont_list.dir' directory
# (usually 'script.dir', that is the same directory containing these functions
# and the main script HORAS).
# Use ont_list_generate() function to recreate or update the ont_list file.

set_specificity <- function(x, ont_list.dir = getwd()) {
  
  # Package for ORA and GSEA in R (by Guangchuang Yu)
  # here used for magrittr pipe operator and dplyr::mutate() function 
  library(clusterProfiler) 
  
  ont_list.path <- paste0(ont_list.dir, "/ont_list.RData")
  
  if (!file.exists(ont_list.path)) {
    stop(paste0("ont_list.RData object not found.\n",
        "Check path or create a new file using ont_list_generate() function"))
  }
  
  # else...
  load(ont_list.path)
  
  # x is an object of class 'enrichResult' (see structure)
  ont <- x@ontology
  
  # Look for the occurrence of each GO term (ID) in each of the 14 levels of the
  # ontology considered here,
  x %<>% # see magrittr::`%<>%` operator to pipe and then assign the result
    mutate(
      # For all GO terms (ID)...
      GO_levels = sapply(ID, function(ID) { # "x@result$ID" implied
        # ...list the levels in which they can be found.
        which(unlist(lapply(ont_list[[ont]], function(GOlevel) ID %in% GOlevel)))
      }),
      # ...then create two other new variables holding just the most general...
      first_GO_level = unlist(lapply(GO_levels, min)),
      # ...and the most specific levels 
      last_GO_level = unlist(lapply(GO_levels, max))
    )
  
  return(x)
}






