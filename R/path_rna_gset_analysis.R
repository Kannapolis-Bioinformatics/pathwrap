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
#' @param entity : scientific name of the species
#' @param exp.fc : log fold expression values
#' @param compare : how the comparison is done for GAGE, see gage for details
#' @param gage.dir : directory in which GAGE will be run
#' @param cnts : counts of genes to use in pathview
#' @param grp.idx : index of the reference and sample for differential analysis
#'
#' @import stats
#' @import utils
#' @import gage
#' @import pathview
#' @return nothing returned
#'

run_pathway <- function(entity, exp.fc, compare, gage.dir, cnts, grp.idx) {
    ref <- which(grp.idx == "reference")
    samp <- which(grp.idx == "sample")

    kegg.gs <- kegg.gsets(entity)
    
    signmetinkegg <- kegg.gs$kg.sets[kegg.gs$sigmet.idx]
    diseaseinkegg <- kegg.gs$kg.sets[kegg.gs$dise.idx]
    siginkegg <- kegg.gs$kg.sets[kegg.gs$sig.idx]
    metainkegg <- kegg.gs$kg.sets[kegg.gs$met.idx]

    # TO DO try using lapply for all function call of kegg pathways
    run_gage(signmetinkegg, sig_n_met, fc_matrix = exp.fc,entity,
                same.dir = FALSE, compare = compare)
    run_gage(diseaseinkegg, disease,fc_matrix = exp.fc,entity,
        same.dir = FALSE,compare = compare )
    run_gage(siginkegg, signalling,fc_matrix = exp.fc,entity,
        same.dir = FALSE,compare = compare   )
    run_gage(metainkegg, metabolism,fc_matrix = exp.fc,entity,
        same.dir = FALSE,compare = compare)
    #GO ANALYSIS
    keggcode_sel <- unname(korg[which(korg[, 4] == entity), 3])
    data(bods, package = "gage", envir = environment())
    common_name_species <- bods[, 2][which(bods[, 3] == keggcode_sel)]

    go.gs <- go.gsets(common_name_species)
    go.bp <- go.gs$go.sets[go.gs$go.subs$BP]
    go.mf <- go.gs$go.sets[go.gs$go.subs$MF]
    go.cc <- go.gs$go.sets[go.gs$go.subs$CC]

    run_gage_analysis(go.bp, biological_process,fc_matrix = exp.fc,
        same.dir = TRUE,compare = compare,entity)
    run_gage(go.mf, molecular_function,same.dir = TRUE,
        compare = compare,fc_matrix = exp.fc, entity)
    run_gage(go.cc, cellular_component,same.dir = TRUE,
        compare = compare, fc_matrix = exp.fc,entity)
    return(invisible(NULL))
}
