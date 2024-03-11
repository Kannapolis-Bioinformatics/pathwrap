run_cpathway<- function( cdatapath,cpd_id_type= "KEGG COMPOUND accession",csamp,
                        cref, ccompare="paired" , cset_dir){
    if (grepl(".csv$", x = cdatapath)){
        cdata <- read.csv(cdatapath, header = T, row.names = 1)
    } else{
        cdata<- read.table(file = cdatapath,sep="\t",header= TRUE,row.names = 1)
    }
    
    if (cpd_id_type!="KEGG COMPOUND accession"){
        cpd_idmap<-cpd2kegg(in.ids=rownames(cdata),in.type=toupper(cpd_id_type))
        didx<-duplicated(cpd_idmap[,1])
        cpd_idmap<-cpd_idmap[!didx,]
        cdata <-mol.sum(cdata, cpd_idmap)
    }
    
    # same.dir if to test for changes in a gene set toward a single direction
    if(!is.null(cref)){
        ncsamp<-length(csamp)
        ncref<-length(cref)
        if(ccompare=="paired" & ncsamp==ncref) 
            {cdata<-cdata[,csamp]- cdata[,cref]}
        else if (ncref==1) cdata<- cdata[,csamp]- cdata[,cref]
        else cdata <-cdata[,csamp]- rowMeans(cdata[,cref])
    }
    csets <- loadcsets()
    gage_return_obj_c <- run_gage(csets, cset_dir , same.dir = FALSE,
            compare = ccompare, fc_matrix = cdata,entity= NA)
    return(gage_return_obj_c)
    
    
}
