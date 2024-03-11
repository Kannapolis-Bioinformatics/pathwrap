#' script to run the combine analysis of compound and gene
#' This is internal function of pathviewwrap
#' 
#' @param gpath_ids top pathways from gene based analysis
#' @param cpath_ids significance compound sets ids
#' @param pgs.gene p-value of analysis of gene set 
#' @param pgs.cpd p-value of analysis of compound set
#' @param somdir directory where data is stored
#' @param gsets pathways on which analysis is based
#' @param gage.out retsult from gage analysis of gene
#' @param gage.out.cpd retsult from gage analysis of compound
#' @param qcut the threshold value to determine significance of gene/cpd sets
#' @return path.ids significant pathids from combined analysis

run_combinedpath_analysis<- function(gpath_ids, cpath_ids,gsets, pgs.gene,
                            pgs.cpd, somdir, gage.out, gage.out.cpd,qcut=0.2) {
    nmax <- 6
    path.ids=c(gpath_ids,cpath_ids)
    nsig.c=0
    if(!is.null(gpath_ids) & !is.null(cpath_ids)){
        pnames=names(gsets)
        pmat=cbind(pgs.gene[pnames], pgs.cpd[pnames])
        log.pmat=-log(pmat)
        nc <- apply(log.pmat,1, function(x) sum(!is.na(x)))
        sg.glob <- apply(log.pmat, 1, sum, na.rm=T)
        pvals <- pgamma(sg.glob, shape = nc, rate = 1, lower.tail = FALSE)
        qvals=p.adjust(pvals, method = "BH")
        colnames(pmat)=paste0("p.", c("gene","cpd"))
        
        gmat=cmat=cbind(gage.out[pnames,1:2], gage.out.cpd[,1][pnames], 
                        gage.out.cpd[,2][pnames])
        colnames(gmat)=paste0(rep(c("stat.mean", "set.size"),2), 
                              rep(c(".gene",".cpd"),each=2))
        
        combo.out=cbind(gmat, pmat, p.global=pvals, q.global=qvals)
        combo.out=combo.out[order(pvals),]
        combo.out[is.nan(combo.out)]=NA
        write.table(combo.out,file=file.path(somdir, "combo.res.tsv"), sep="\t", 
                    col.names=NA,  quote=T)
        
        sig.c=combo.out[,"q.global"]<qcut & !is.na(combo.out[,"q.global"])
        nsig.c=sum(sig.c, na.rm=T)
        if(nsig.c>0){
            combo.out.sig=data.frame(combo.out)[sig.c,]
            write.table(combo.out.sig, file = file.path(somdir, "combo.res.sig.tsv"), sep="\t",
                        col.names=NA, quote = FALSE)
            path.ids=rownames(combo.out.sig)
        }
    }
 

    path.ids=gsub(".+([0-9]{5}).+", "\\1",path.ids)
    globs=grep("^01[1-2]", path.ids)
    #rm.glob=length(globs)>0 & (!is.null(gene.d) | !args2$kegg)
    #if(rm.glob) path.ids=path.ids[- globs]
    if(length(globs)>0) path.ids=path.ids[- globs]
    path.ids=unique(path.ids)
    if(length(path.ids)>nmax & nsig.c<1){
        lgp=length(gpath_ids)
        lcp=length(cpath_ids)
        mn=min(lgp,lcp,round(nmax/2))
        path.ids=c(gpath_ids[1:mn],cpath_ids[1:mn])
        if(lgp>mn) path.ids=c(path.ids,gpath_ids[(mn+1):lgp])
        if(lcp>mn) path.ids=c(path.ids,cpath_ids[(mn+1):lcp])
        path.ids=gsub(".+([0-9]{5}).+", "\\1",path.ids)
        globs=grep("^01[1-2]", path.ids)
        #    rm.glob=length(globs)>0 & (!is.null(gene.d) | !args2$kegg)
        #    if(rm.glob) path.ids=path.ids[- globs]
        if(length(globs)>0) path.ids=path.ids[- globs]
        path.ids=unique(path.ids)
    }
    path.ids=path.ids[1:min(nmax, length(path.ids))]
    return(path.ids)
}