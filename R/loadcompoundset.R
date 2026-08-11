#' Function to generate compound set for pathway analysis
#' 
#' This function list the pathways along with the list of compounds all the
#'  pathway contains
#' @param entitycode three char code of organism for kegg
#' @import KEGGREST
#' @return compound_sets is returned


# loadcsets <- function() {
#     #list of kegg pathway linked from compound databases
#     allpathwayncompound<-keggLink("pathway", "compound") #ref pathways for cpd
#     names(allpathwayncompound) <- str_remove(pattern = "cpd:", 
#                                             names(allpathwayncompound)  )
#     compound_sets <- split(names(allpathwayncompound), 
#                         unname(allpathwayncompound))
#     names(compound_sets)<- str_remove(pattern = "path:",names(compound_sets) )
# 
#     pathway_names <- unname(keggList("pathway"))[match(names(compound_sets),
#                                                 names(keggList("pathway")))]
#                             
#     names(compound_sets) <- str_replace_all(names(compound_sets), "$", 
#                                             paste0(" ",pathway_names))
#     return(compound_sets)
# }


#maps every KEGG pathway available in mouse to the list of KEGG compounds (metabolites) that participate in that pathway.

loadcsets <- function(entitycode) {
    mmu_pathways <- keggList("pathway", entitycode) # a master list of all available pathway IDs. 
    mmu_ids <- str_extract(names(mmu_pathways), "[0-9]{5}")
    allpathwayncompound <- keggLink("pathway", "compound") # which pathway contain which compound
    names(allpathwayncompound) <- str_remove(names(allpathwayncompound), "cpd:")
    
    compound_sets <- split(names(allpathwayncompound),
                           unname(allpathwayncompound))
    names(compound_sets) <- str_remove(names(compound_sets), "path:")
    
    # --- Filter to mouse-available pathways only ---
    map_ids <- str_extract(names(compound_sets), "[0-9]{5}")
    compound_sets <- compound_sets[map_ids %in% mmu_ids]
    
    # Re-attach names
    pathway_names <- unname(keggList("pathway"))[
        match(str_extract(names(compound_sets), "[0-9]{5}"),
              str_extract(names(keggList("pathway")), "[0-9]{5}"))
    ]
    names(compound_sets) <- paste0(
        str_extract(names(compound_sets), "[0-9]{5}"),
        " ", pathway_names
    )
    names(compound_sets) <- paste0(entitycode, names(compound_sets))
    
    return(compound_sets)
}
