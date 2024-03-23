#' Check and plot pathways with the data
#' 
#' This function checks if pathways can be downloaded from the kegg and 
#' plots by rendering the data
#' @import pathview
#' @param gage.dir directory in which the pathways are saved
#' @param entity organism whose pathway is plotted
#' @param path.ids list of pathway ids to plot
#' @param cpd_data compound data to plot
#' @param gene_data gene information to be plotted in pathview
#' @param gene_id_type what type of gene id is used in gene_data
#' @param cpd.idtype what type of cpd.id is used in cpd.data
#' @return invisiblenull

plotpathways <- function(gage.dir, entity, path.ids, gene_data, cpd_data, 
                        cpd.idtype = "kegg", gene_id_type = "entrez") {
    check_and_warn <- function(condition, message) {
        if (condition) {
            warning(message)
        }
    }
    species_code <- kegg.species.code(entity)
    for (pid in na.omit(path.ids[seq_len(6)])) {
        message(c("Plotting pathview for ", pid, collapse=""))
        tryCatch({
            download.kegg(kegg.dir = gage.dir, pathway.id = pid, 
                        species = species_code)
            pathview(gene.data = gene_data, pathway.id = pid, 
                    species = species_code,
                    out.suffix = "pathview", kegg.dir = gage.dir, 
                    cpd.data = cpd_data,
                    cpd.idtype = cpd.idtype, gene_id_type = gene_id_type)
            
        }, error = function(w) {
            check_and_warn(TRUE, paste("ERROR: Pathview failed on", pid))
        })
        Files <- list.files(path = getwd(), full.names = TRUE, 
                            pattern = pid)
        
        if (length(Files) != 0) {
            newName <- gsub(dirname(Files)[1], gage.dir, Files)
            file.rename(Files, newName)   
        }
    }
    return(invisible(NULL))
}
