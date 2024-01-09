#' Run GAGE and Pathview
#'
#' `run_pathway` runs GAGE for pathway analysis
#' GAGE is based upon the user supplied comparision
#' method for the species specified.
#' The biological process, cellular component and
#' molecular function analysis for GO terms are done separately.
#'
#' KEGG disease, KEGG signalling and metabolism pathways are analysed separately
#' Top enriched pathways with  "q.val" < 0.1 are visualized using pathview.
#'
#' @param gsets: gene sets to analyse
#' @param work.dir : directory where results will be stored
#' @param same.dir : if the direction is same for GAGE analysis, GAGE parameter
#' @param compare : GAGE parameter
#'
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'
run_gage <- function(gsets, work.dir, same.dir,compare = compare,
                     fc_matrix = exp.fc, entity)) {
    fc.kegg.p <- gage( exp.fc, gsets = gsets, ref = NULL,samp = NULL,
                       same.dir = same.dir, compare = compare)
    sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(
        fc.kegg.p$greater[, "q.val"])
    path.ids <- rownames(fc.kegg.p$greater)[sel]
    anla_type <- "KEGG"
    if (same.dir == TRUE) {
        anla_type <- "GO"
        gage.dir <- file.path(gage.dir , "GO",fsep = .Platform$file.sep)
        sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(
            fc.kegg.p$less[, "q.val"])
        path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
        write.table(fc.kegg.p$less,
                    sep = "\t", file = file.path(work.dir, paste0(
    "fc.", anla_type, ".p.less.txt"), fsep = .Platform$file.sep))
        path.ids <- c(path.ids[seq_len(3)], path.ids.l[seq_len(3)])
    }
    path.ids <- substr(path.ids, 1, 8)
    write.table(fc.kegg.p$greater, sep = "\t",
                file = file.path(work.dir, paste0(
    "fc.", anla_type,"p.greater.txt"), fsep = .Platform$file.sep))
    # visualize top 3 pathways
    # run pathview only for KEGG pathways
    if (same.dir == FALSE) {
        gage.dir <- file.path(gage.dir , "KEGG",fsep = .Platform$file.sep)
        message(paste0("STEP 7: visualizing the pathway", " in ", entity,
                       collapse=""))
        for (pid in na.omit(path.ids.2[1:6])) {
            tryCatch({
                message(paste0("Plotting pathview for ", pid))
                    pathview::pathview(gene.data = exp.fc,pathway.id = pid,
                        species = entity,out.suffix = "pathview")
                }, error = function(e) {
                    message(c("ERROR: Pathview failed on", pid))
                })
        Files <- list.files(path = getwd(),  full.names = TRUE,
                            pattern =paste0(pid, ".pathview.png$"))
        newName <- gsub(dirname(Files), work.dir, Files)
        file.rename(Files, newName)
        }
    }
    kegg.sig <- sigGeneSet(fc.kegg.p,
    outname = paste0(entity, anla_type, ".sig", basename(work.dir)
                           ), pdf.size = c(17, 17), heatmap = FALSE)
    # wont give heatmap for fold change used in gage
    write.table(kegg.sig$greater,
                file = file.path(gage.dir, paste0(anla_type, ".sig.txt"),
                                 fsep = .Platform$file.sep),sep = "\t")
}
