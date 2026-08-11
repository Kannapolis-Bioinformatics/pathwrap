#' Check and plot pathways with the data
#' 
#' This function checks if pathways can be downloaded from the kegg and 
#' plots by rendering the data
#' @import pathview
#' @param kegg.dir directory in which the pathways are saved
#' @param entity organism whose pathway is plotted
#' @param path.ids list of pathway ids to plot c("00056","00098")
#' @param cpd_data compound data to plot
#' @param gene_data gene information to be plotted in pathview
#' @param gene_id_type what type of gene id is used in gene_data
#' @param cpd_id_type what type of cpd.id is used in cpd.data
#' @return invisiblenull

plotpathways <- function(kegg.dir, entity, path.ids, gene_data, cpd_data, 
                         cpd_id_type = "kegg", gene_id_type = "entrez") {
    check_and_warn <- function(condition, message) {
        if (condition) {
            warning(message)
        }
    }
    species_code <- kegg.species.code(entity)
    message("Working in ", kegg.dir , collapse="")
    for (pid in na.omit(path.ids) ){
            if (!file.exists(file.path(kegg.dir, paste0(unname(species_code),pid, ".pathview.png")))){
            message(c("Plotting pathview for ", pid, collapse=""))
            tryCatch({
                download.kegg(kegg.dir = kegg.dir, pathway.id = pid, #files are downloaded in kegg.dir
                              species = species_code)
                pathview(gene.data = gene_data, pathway.id = pid, #.pathview is written in getwd()
                         species = species_code,
                         out.suffix = "pathview", kegg.dir = kegg.dir, 
                         cpd.data = cpd_data)
                #,keys.align="y",match.data=F,multi.state= T, 
                 #        same.layer=T,  cpd.idtype = cpd_id_type, gene_id_type = gene_id_type,kegg.native= TRUE, split.group=TRUE
                 #       )
                
            }, error = function(w) {
                check_and_warn(TRUE, paste("Pathview failed on", pid))
            })
            Files <- list.files(path = getwd() ,pattern=".pathview.", full.names = TRUE, 
                                )
            print(Files)
            
            if (length(Files) != 0) {
                newName <- gsub(dirname(Files)[1], kegg.dir, Files)
                file.rename(Files, newName)   
                message(newName , " is created")
            }
        }
    }
    return(invisible(NULL))
}
