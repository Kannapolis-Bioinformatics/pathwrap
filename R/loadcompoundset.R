#' Function to generate compound set for pathway analysis
#' 
#' This function list the pathways along with the list of compounds all the
#'  pathway contains
#' @import KEGGREST
#' @return compound_sets is returned


loadcsets <- function() {
    #list of kegg pathway linked from compound databases
    allpathwayncompound<-keggLink("pathway", "compound") #ref pathways for cpd
    names(allpathwayncompound) <- str_remove(pattern = "cpd:", 
                                            names(allpathwayncompound)  )
    compound_sets <- split(names(allpathwayncompound), 
                        unname(allpathwayncompound))
    names(compound_sets)<- str_remove(pattern = "path:",names(compound_sets) )

    pathway_names <- unname(keggList("pathway"))[match(names(compound_sets),
                                                names(keggList("pathway")))]
                            
    names(compound_sets) <- str_replace_all(names(compound_sets), "$", 
                                            paste0(" ",pathway_names))
    return(compound_sets)
}
