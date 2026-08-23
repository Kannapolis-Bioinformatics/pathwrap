#' process phenofile 
#' this function extracts information from user supplied phenofile ;
#' information include samplename, location of raw files, if the reads are paired, 
#' which samples belong to which class and which sample is pair of other sample
#'
#' @param phenofile path to the phenofile where raw data files path is stored 
#' @return a list of things
#'
process_phenofile <- function(phenofile){
    
    if (!file.exists(phenofile)) { ### TO DO make sure reference is first ANum
        message("Please provide phenofile with Class information")}
    coldata <- read.table(phenofile, sep = "\t", header = TRUE)
    if (colnames(coldata)[ncol(coldata)] != "Class") {
        message("Please make sure class information is in last column with
                colname 'Class' if you want to run differential analysis")
    
    }
    if (ncol(coldata)==3){
        if (!all(colnames(coldata) %in% c( "SampleName", "FileName","Class"))){
            message("the phenofile column names must be ")
            message("SampleName ", "FileName ",   "Class " )}
        
    } else{
        if (!all(colnames(coldata)%in%c("SampleName ","FileName1 ","FileName2 ",
                                        "Class"  ))){
            message("the phenofile column names must be ") 
            message("SampleName ", "FileName1 ", "FileName2 ", "Class")}
            
    }
    
    coldata$Class <- as.factor(coldata$Class)
    paired_info <- coldata$PairedInfo
    SampleName <- coldata$SampleName
    if ("PairedInfo"%in%colnames(coldata)){
        coldata<- subset(coldata, select = -c(PairedInfo) )
    }
    filenames <- as.data.frame(coldata[, -c(1, ncol(coldata))])
    if (dim(filenames)[2] == 1) {
        endness <- "SE"
        fq.dir <- dirname(filenames[1, 1])
        colnames(filenames) <-  "FileName"  
    } else if (dim(filenames)[2] == 2) {
        endness <- "PE"
        fq.dir <- dirname(filenames$FileName1[1]) 
    } else if (dim(filenames)[2] > 5) {
        message('Please make sure there are max 5 columns with colnames')
        message('"SampleName" "FileName1" "FileName2" "PairedInfo" "Class"')
        message("if the sample are paired")
    } else{
        endness<- NA
        fq.dir <- NA
    }
    return(list("endness"  =endness , "fq.dir" = fq.dir, "SampleName"= SampleName,
                "FileName" = filenames, "paired_info"= paired_info , "coldata" = coldata ))
    
}