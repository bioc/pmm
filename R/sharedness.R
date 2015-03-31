sharedness <- function(fit, threshold = 0.2, na.action = "na.omit"){

    ## NA-Handling for genes that were not measured for all conditions
    if(na.action == "na.omit") fit <- na.omit(fit)

    ## Checking if input is a result of pmm
    if(is.null(attr(fit,"doc"))){
        stop("Your input matrix is not a result of pmm.")
    }else if(attr(fit,"doc")!= "PMM.Result"){
        stop("Your input matrix is not a result of pmm.")
    }

    ## Extracting the conditions from the fit-object.
    condition <- colnames(fit)[grep("ccg",colnames(fit))]
    condition <- sapply(condition,
                        function(name) strsplit(name,split = "ccg.")[[1]][2])

    ## Selecting all hits i.e. genes that have at least in one condition
    ## fdr below the given threshold.
    fdr <- fit[,grep("fdr",colnames(fit))]
    if(na.action == "use"){
        for(p in condition){
            nmiss <- sum(is.na(fdr[,paste("fdr",p,sep = ".")]))
            if(nmiss > 0)
                fdr[is.na(fdr[,paste("fdr",p,sep = ".")]),
                    paste("fdr",p,sep = ".")] <- rep(2,nmiss)
        }
    }

    hit.genes <- function(fdr, threshold){
        fit$GeneID[fdr < threshold]
    }
    hits <- unique(unlist(apply(fdr,2,hit.genes,threshold = threshold)))

    ## Constructing a data.frame with all genes that have at least in
    ## one condition fdr below the given threshold.
    d.hit <- NULL
    for(p in condition){
        add <- fit[fit[,"GeneID"] %in% hits,
                   c(paste("ccg",p,sep = "."),
                     paste("fdr",p,sep = "."))]
        colnames(add) <- c("ccg","fdr")
        d.hit <- rbind(d.hit,cbind(condition = rep(p,length(hits)),
                                   GeneID = sort(hits),add))
    }

    ## If there are no hits, then output an error.
    if(nrow(d.hit) == 0)
        stop("No genes have FDR below the given threshold.")

    ## Calculating the sharedness for all genes.
    hit.mean <- tapply(d.hit$fdr,d.hit$GeneID,mean, na.rm = TRUE)
    hit.prop <- tapply(d.hit$fdr,d.hit$GeneID,
                       function(x) sum(x<1,na.rm = TRUE)/length(x))
    sh <- (hit.prop + 1 - hit.mean)/2

    ## Generating one data.frame with all results
    sh <- data.frame(GeneID = names(sh), Sharedness = as.numeric(sh))

    ## Output: Returns sharedness
    return(sh)
}
