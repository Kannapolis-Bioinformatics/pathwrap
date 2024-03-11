#' Function to generate compound set for pathway analysis
#' 
#' This function list the pathways along with the list of compounds the pathway contains
#' @import KEGGREST
#' @return compound_sets is returned


loadcsets <- function() {
    allpathwayncompound<-keggLink("pathway", "compound") #reference pathways for compound
    names(allpathwayncompound) <- str_remove(pattern = "cpd:", names(allpathwayncompound)  )
    #allpathwayncompound <- str_remove(pattern = "path:", unname(allpathwayncompound))
    
    compound_sets <- split(names(allpathwayncompound), unname(allpathwayncompound))
    names(compound_sets)<- str_remove(pattern = "path:",names(compound_sets) )
    return(compound_sets)
}

