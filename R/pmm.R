pmm <- function(df.data, response, weight = "None", ignore = 3, simplify = TRUE,
                gene.col = "GeneID", condition.col = "condition")
{

    ## CHECKING INPUT ##
    if(!is.data.frame(df.data))
        stop("The data must be given as a data.frame.")
    if(!is.character(gene.col))
        stop("The argument for gene.col needs to be a character.")
    if(!gene.col %in% colnames(df.data))
        stop(paste("The data needs to contain a column that indicates",
                   "the genes. Specify the column name in the argument",
                   "gene.col."))
    if(length(unique(df.data[,gene.col])) < 2)
        stop("The data needs to have at least two different genes.")
    if(!is.character(condition.col))
        stop("The argument for condition.col needs to be a character.")
    if(!condition.col %in% colnames(df.data))
        stop(paste("The data needs to contain a column that indicates",
                   "the condition. Specify the column name in the",
                   "argument condition.col."))
    if(length(unique(df.data[,condition.col])) < 2)
        stop("The data needs to have at least two different conditions.")
    if(!is.character(response))
        stop("The argument for response needs to be a character.")
    if(!response %in% colnames(df.data))
        stop(paste("The data frame does not contain a column named",
                   "as the given response argument."))



    ## PREPARATION ##
    ## Reading Data
    df.data$y <- df.data[,response]     # response
    df.data <- df.data[!is.na(df.data$y),] # removes rows with NAs.

    ## Removing genes with less then k replicates (k = ignore)
    rm.reps <- function(df.data,n){
        df.data$IDxxx <- paste(df.data[,condition.col],df.data[,gene.col],
                            sep = ":")
        t1 <- table(df.data$IDxxx)
        df.data <- df.data[df.data$IDxxx %in% names(t1)[t1 >= n],]
        return(droplevels(df.data))
    }
    df.data <- rm.reps(df.data,ignore)
    if(nrow(df.data) == 0)
        stop("The number of measurements per gene and conditions is less than ",
             ignore," repetitions. Ignore was set too high.")


    ## Defining factors
    df.data$condition <- as.factor(df.data[,condition.col])
    df.data$GeneID <- as.factor(df.data[,gene.col])  # GeneIDs
    condition <- levels(df.data$condition)          # conditions
    if(length(condition) <= 1 | length(levels(df.data$GeneID)) <= 1)
        stop(paste("There has too be more than one gene and one condition",
                   "to fit the pmm.",
                   "(Check your data and the argument ignore.)"))

    ## Checking if weights are required
    if(!is.character(weight))
        stop("The argument weight needs to be a character")
    if(!(weight %in% colnames(df.data)) & weight != "None")
        stop("There is no column with the weights.")
    if(!weight %in% colnames(df.data)){
        df.data$NOweight <- rep(1,nrow(df.data))
        weight <- "NOweight"
    }


    ## FIT ##
    ## Fitting the linear mixed model
    fit.lmm <- lmer(y ~ condition + (1|GeneID) + (1|condition:GeneID),
                    data = df.data,                   # data
                    weights = df.data[,weight],       # weights
                    verbose = FALSE)


    ## EXTRACTING EFFECTS ##
    ## Extracting a_g
    ag <- data.frame(ag = ranef(fit.lmm)[["GeneID"]][,1],
                     GeneID = rownames(ranef(fit.lmm)[[2]]))
    ## Extracting b_cg
    tp1 <- function(name) substr(name,gregexpr(":",name)[[1]]+1,nchar(name))
    bcg <- data.frame(bcg = ranef(fit.lmm)[["condition:GeneID"]][,1],
                      GenePathID = rownames(ranef(fit.lmm)[[1]]))
    bcg$GeneID <- sapply(as.character(bcg$GenePathID),tp1)
    ## Calculating c_cg = a_g + b_cg
    ccg <- merge(bcg,ag,by = "GeneID")
    ccg$ccg <- ccg$ag + ccg$bcg
    ## Generating a data.frame with the results
    tp2 <- function(name) substr(name,1,gregexpr(":",name)[[1]]-1)
    ccg$condition <- sapply(as.character(ccg$GenePathID),tp2)
    ccg <- ccg[,c("GeneID","condition","ag","bcg","ccg")]
    ccg.matrix <- reshape(ccg[,c("condition","GeneID","ccg")],
                          direction = "wide",
                          timevar = "condition", idvar = "GeneID")

    ## FALSE DISCOVERY RATE  ##
    genes <- ccg.matrix$GeneID
    ## for each condition we calculate seperately the fdr
    for(i in condition){
        ## Removing genes that are not measured for condition i
        sub.ccg.matrix <- na.omit(ccg.matrix[,c("GeneID",
                                    paste("ccg.",i,sep=""))])
        ## Computing false discovery rate       
        dat <- local.fdr(sub.ccg.matrix[,paste("ccg.",i,sep="")])
        ## Adding fdr to ccg values
        sub.ccg.matrix <- cbind(sub.ccg.matrix,fdr=dat) 
        ## Adding genes again that are not measured for condition i
        nmiss <- sum(is.na(ccg.matrix[,paste("ccg.",i,sep="")]))
        dmiss <- ccg.matrix[is.na(ccg.matrix[,paste("ccg.",i,sep="")]),
                            c("GeneID",paste("ccg.",i,sep=""))]
        sub.ccg.matrix <- rbind(sub.ccg.matrix,
                                cbind(dmiss,fdr = rep(NA,nmiss)))
        ## Setting a unique name for column containing the fdr
        colnames(sub.ccg.matrix)[colnames(sub.ccg.matrix)=="fdr"] <-
            paste("fdr",i,sep=".")
        ## Adding new column to ccg.matrix
        ccg.matrix <- merge(ccg.matrix,
                  sub.ccg.matrix[,c("GeneID",paste("fdr.",i,sep=""))],
                              by = "GeneID")
    }
    ## Sorting columns by conditions
    ccg.matrix <- ccg.matrix[,c("GeneID",
                                paste(rep(c("ccg","fdr"),length(condition)),
                                        rep(condition,each = 2),sep="."))]


    ## OUTPUT ##
    ## returns the c_pg matrix and false discovery rate
    attr(ccg.matrix,"doc") <- "PMM.Result"
    if(simplify)
        return(ccg.matrix)
    else
        return(list(ccg.matrix = ccg.matrix, lmm = fit.lmm,ccg = ccg))
}
