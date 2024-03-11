#' Check and plot pathways with the data
#' 
#' This function checks if pathways can be downloaded from the kegg and 
#' plots by rendering the data
#' @import pathview
#' @param gage.dir directory in which the pathways are saved
#' @param entity organism whose pathway is plotted
#' @param path.ids list of pathway ids to plot
#' @param fc_matrix genedata to plot
#' @param cpd_data compound data to plot
#' @return invisiblenull

plotpathways <- function(gage.dir,entity,path.ids, fc_matrix,cpd_data = NULL){
    #gage.dir <- file.path(gage.dir , "KEGG",fsep = .Platform$file.sep)
    message(paste0("STEP 7: visualizing the pathway", " in ", entity,
                collapse=""))
    for (pid in na.omit(path.ids[seq_len(6)])){
        tryCatch({
            message(paste0("Plotting pathview for ", pid, collapse=""))
            pathview(gene.data = fc_matrix,pathway.id = pid,
                            species = entity,out.suffix = "pathview")
            Files <- list.files(path = getwd(),  full.names = TRUE,pattern =pid)
            if (length(Files)!=0){
                newName <- gsub(dirname(Files)[1], gage.dir, Files)
		print(paste0("renaming file", Files, " to ", newName ))

                file.rename(Files, newName)
		}
            
        }, error = function(e) {
             message(c("ERROR: Pathview failed on", pid, collapse=""))
        }) }
    return(invisible(x=NULL))
}
