#' Run GAGE and Pathview
#'
#' This function runs generally applicable gene set enrichment for pathway analysis :
#' GAGE based upon the user supplied comparision method for the species specified.
#' The biological process, cellular component and molecular function analysis for GO terms are done seperately.
#' Also, KEGG disease and KEGG signalling and metabolism pathways are analysed seperately.
#' Finally the top enriched pathways with  "q.val" < 0.1 are visualized using pathview.
#'
#' @param entity : scientific name of the species
#' @param exp.fc : log fold expression values
#' @param compare : how the comparision is done for GAGE, see gage for details
#' @param gage.dir : directory in which GAGE will be run
#' @param cnts : counts of genes to use in pathview
#' @param grp.idx : index of the reference and sample for differential analysis
#'
#' @import stats
#' @import utils
#' @import gage
#' @import pathview
#'
#'

run_pathway <- function(entity, exp.fc, compare, gage.dir, cnts, grp.idx) {
    ref <- which(grp.idx == "reference")
    samp <- which(grp.idx == "sample")

    kegg.gs <- kegg.gsets(entity)
    run_gset_analysis <- function(gsets, work.dir, same.dir, compare = compare) {
        # setwd(work.dir)
        fc.kegg.p <- gage(exp.fc, gsets = gsets, ref = NULL, samp = NULL, same.dir = same.dir, compare = compare)
        sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
        path.ids <- rownames(fc.kegg.p$greater)[sel]
        anla_type <- "KEGG"
        if (same.dir == TRUE) {
            anla_type <- "GO"
            sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[, "q.val"])
            path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
            write.table(fc.kegg.p$less, file = file.path(work.dir, paste0("fc.", anla_type, ".p.less.txt")), sep = "\t")
            path.ids <- c(path.ids[seq_len(3)], path.ids.l[seq_len(3)])
        }

        path.ids <- substr(path.ids, 1, 8)
        # path.ids <- gsub("[^0-9.-]", "", sapply(stringr::str_split(path.ids, " ", 2),"[[",1))
        write.table(fc.kegg.p$greater, file = file.path(work.dir, paste0("fc.", anla_type, "p.greater.txt")), sep = "\t")
        # visualize top 3 pathways
        # run pathview only for KEGG pathways
        if (same.dir == FALSE) {
            message(paste0("STEP 7: visualizing the pathway", " in ", entity))
            pv.out.list <- sapply(na.omit(path.ids[seq_len(6)]), function(pid) pathview(gene.data = exp.fc, pathway.id = pid, species = entity, out.suffix = paste0(entity, pid)))
        }
        kegg.sig <- sigGeneSet(fc.kegg.p, outname = paste0(entity, anla_type, ".sig", basename(work.dir)), pdf.size = c(17, 17), heatmap = FALSE) # wont give heatmap for fold change used in gage
        write.table(kegg.sig$greater, file = file.path(gage.dir, paste0(anla_type, ".sig.txt")), sep = "\t")
    }

    # TO DO try using lapply for all function call of kegg pathways
    run_gset_analysis(kegg.gs$kg.sets[kegg.gs$sigmet.idx], sig_n_met, same.dir = FALSE, compare = compare)
    run_gset_analysis(kegg.gs$kg.sets[kegg.gs$dise.idx], disease, same.dir = FALSE, compare = compare)
    run_gset_analysis(kegg.gs$kg.sets[kegg.gs$sig.idx], signalling, same.dir = FALSE, compare = compare)
    run_gset_analysis(kegg.gs$kg.sets[kegg.gs$met.idx], metabolism, same.dir = FALSE, compare = compare)
    ############################################################################################

    keggcode_sel <- unname(korg[which(korg[, 4] == entity), 3])
    data(bods, package = "gage", envir = environment())
    common_name_species <- bods[, 2][which(bods[, 3] == keggcode_sel)]

    go.gs <- go.gsets(common_name_species)
    go.bp <- go.gs$go.sets[go.gs$go.subs$BP]
    go.mf <- go.gs$go.sets[go.gs$go.subs$MF]
    go.cc <- go.gs$go.sets[go.gs$go.subs$CC]

    run_gset_analysis(go.bp, biological_process, same.dir = TRUE, compare = compare)
    run_gset_analysis(go.mf, molecular_function, same.dir = TRUE, compare = compare)
    run_gset_analysis(go.cc, cellular_component, same.dir = TRUE, compare = compare)
}
