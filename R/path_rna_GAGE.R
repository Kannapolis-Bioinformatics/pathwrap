#' Run GAGE and Pathview
#'
#' `run_pathway` runs GAGE for pathway analysis
#' GAGE is based upon the user supplied comparision
#' method for the species specified.
#' The biological process, cellular component and
#' molecular function analysis for GO terms are done separately.
#'
#' KEGG disease, KEGG signalling and metabolism pathways are analysed separately
#' Top enriched pathways with    "q.val" < 0.1 are visualized using pathview.
#'
#' @param gsets : gene sets to analyse
#' @param work.dir : directory where results will be stored
#' @param same.dir : if the direction is same for GAGE analysis, GAGE parameter
#' @param compare : GAGE parameter
#' @param    fc_matrix fold change values
#' @param entity organism of study, scientific name
#' @import gage
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'
run_gage <- function(gsets, work.dir, same.dir, fc_matrix,
        compare = compare, entity) {
        exp.fc<- fc_matrix
        #set workdir
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
                        fc.kegg.p$less[, "q.val"]
                )
                path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
                write.table(fc.kegg.p$less,
                                        sep = "\t",
                                        file = file.path(work.dir, paste0(
        "fc.", anla_type, ".p.less.txt"), fsep = .Platform$file.sep))
                path.ids <- c(path.ids[seq_len(3)], path.ids.l[seq_len(3)])
        }
        path.ids <- substr(path.ids, 1, 8)
        write.table(fc.kegg.p$greater,
                                sep = "\t",
                                file = file.path(work.dir, paste0(
        "fc.", anla_type,"p.greater.txt"), fsep = .Platform$file.sep))
        # visualize top 3 pathways
        # run pathview only for KEGG pathways
        if (same.dir == FALSE) {
            gage.dir <- file.path(gage.dir , "KEGG",fsep = .Platform$file.sep)
            message(paste0("STEP 7: visualizing the pathway", " in ", entity,
            collapse=""))
                pv.out.list <- vapply(
                        na.omit(path.ids[seq_len(6)]),
                        function(pid) {
                                pathview(kegg.dir = work.dir,gene.data = exp.fc,
                                pathway.id = pid, species = entity, 
                                out.suffix = paste0(entity, pid))[1]
                        }, data.frame(length(na.omit(path.ids[seq_len(6)]))))}
        kegg.sig <- sigGeneSet(fc.kegg.p,
        outname = paste0(entity, anla_type, ".sig", basename(work.dir)
            ), pdf.size = c(17, 17), heatmap = FALSE)
        # wont give heatmap for fold change used in gage
        write.table(kegg.sig$greater,
            file = file.path(gage.dir, paste0(anla_type, ".sig.txt"),
            fsep = .Platform$file.sep),sep = "\t")
}
