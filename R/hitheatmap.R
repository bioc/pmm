hitheatmap <- function(fit, threshold = 0.2, sharedness.score = FALSE,
                       main = "", na.action = "use", ...)
{

    ## NA-Handling for genes that were not measured for all conditions
    if(na.action == "na.omit") fit <- na.omit(fit)

    ## Check if input is a result of pmm
    if(is.null(attr(fit,"doc"))){
        stop("Your input matrix is not a result of pmm.")
    }else if(attr(fit,"doc")!= "PMM.Result"){
        stop("Your input matrix is not a result of pmm.")
    }

    ## Selects hits from fit, i.e. genes that have fdr values below the
    ## given threshold.
    fdr <- fit[,c(grep("fdr",colnames(fit)))]
    hits <- function(fdr, threshold){
        below.th <- apply(fdr,2,function(x) which(x < threshold))
        return(fit[unique(unlist(below.th)),"GeneID"])
    }
    hit.genes <- hits(fdr,threshold)

    ## Saving ccg values of hit genes in one data.frame
    d.hit <- fit[fit[,"GeneID"] %in% hit.genes,
                 c(1,grep("ccg",colnames(fit)))]
    d.hit <- d.hit[order(apply(d.hit[,grep("ccg",colnames(d.hit))],1,mean,
                               na.rm = TRUE)),]

    ## If there are no hits, then output an error.
    if(nrow(d.hit) == 0)
        stop("No genes have FDR below the given threshold.")

    ## Caculating Sharedness
    if(sharedness.score){
        sh <- sharedness(fit,threshold, na.action = na.action)
    }

    ## Plotting Options
    par(mar = c(5,10,4.1,6), ...)

    ## Define color for up and down hits and sharedness
    ramp <- colorRamp(c("white","indianred","red","red2","red4","darkred"))
    red <- rgb(ramp(seq(0, 1, length = 102)), maxColorValue = 255)
    ramp <- colorRamp(c("white","lightblue", "blue","blue2","blue4","darkblue"))
    blue <- rgb(ramp(seq(0, 1, length = 102)), maxColorValue = 255)
    ramp <- colorRamp(c("white","lightgreen","green","green2","green4",
                        "darkgreen"))
    green <- rgb(ramp(seq(0, 1, length = 102)), maxColorValue = 255)


    ## Number of genes and conditions
    condition <- unique(substring(colnames(d.hit)[2:ncol(d.hit)],5))
    N <- nrow(d.hit)
    cn <- length(condition)

    ## Prepare Plot Window
    sy <- 1
    if(sharedness.score) sy <- 0
    plot(y = c(sy-0.5,cn+0.5), x = c(1,N+3.1), type = "n",
         bty ="n", xaxt ="n", yaxt = "n", main = main,
         xlab = "", ylab = "")

    ## Plotting Results
    for(j in 1:cn){
        cc <- condition[cn - j + 1]
        for(g in 1:N){
            ind <- round(100/(max(abs(d.hit[,-1]),na.rm = TRUE))*
                             abs(fit[fit$GeneID == d.hit$GeneID[g],
                                     paste("ccg",cc,sep=".")]),0)
            if(!is.na(ind)){
                if(ind != 0){
                    color <- ifelse(fit[fit$GeneID == d.hit$GeneID[g],
                                        paste("ccg",cc,sep=".")] > 0,
                                    red[ind], blue[ind])
                }else{
                    color <- "white"
                }
                marker <- ifelse(fit[fit$GeneID == d.hit$GeneID[g],
                                     paste("fdr",cc,sep=".")] < threshold,
                                 "yellow",color[g])
                markersize <- ifelse(fit[fit$GeneID == d.hit$GeneID[g],
                                         paste("fdr",cc,sep=".")] < threshold,
                                     1.2,0.1)
                polygon(x = c(g - 0.5, g - 0.5, g + 0.5, g + 0.5),
                        y = c(j - 0.5, j + 0.5, j + 0.5, j - 0.5),
                        border = "white", col = color, lwd = 1.2)
                points(g,j,pch = 8, cex = markersize, col = marker)
            }else{
                color <- "white"
                polygon(x = c(g - 0.5, g - 0.5, g + 0.5, g + 0.5),
                        y = c(j - 0.5, j + 0.5, j + 0.5, j - 0.5),
                        border = "white", col = color, lwd = 1.2)
                text(g,j,"NA", cex = 0.4)
            }
        }
    }

    ## Adding Sharedness
    if(sharedness.score){
        for(g in 1:N){
            ind <- 100*abs(sh[sh$GeneID == d.hit$GeneID[g],"Sharedness"])
            polygon(x = c(g - 0.5, g - 0.5, g + 0.5, g + 0.5),
                    y = c(-0.3, 0.3, 0.3,  -0.3),
                    border = "white", col = green[ind],lwd = 2)
        }
    }

    ## Adding Legend
    col.legend <- colorRampPalette(c("darkblue","blue4","blue2","blue",
                                     "lightblue","white","indianred",
                                     "red","red2","red4","darkred"))(200)
    color.bar <- function(lut, min, max, max.scale = NULL) {
        scale <- (length(lut)-1)/(max-min)
        if(!is.null(max.scale)){
            axis(4, round(seq(min, max, len=3),2),
                 c(-round(max.scale,1),0,round(max.scale,1)),las=1,tck = -0.01)}
        for (i in 1:(length(lut)-1)) {
            y <- (i-1)/scale + min
            rect(N+2.5,y,N+3.5,y+1/scale, col=lut[i], border=NA)
        }
    }
    color.bar(col.legend,1/3*cn,cn,max(abs(d.hit[,-1]),na.rm = TRUE))
    axis(4,cn+0.5,"c_cg values",las=2,tck=0,mgp=c(0,-1.5,0))
    points(N+3, 1/4*cn, pch=8, col = "yellow", cex = 1.2)
    axis(4,1/4*cn,paste("FDR <",threshold),las=2,tck=0,mgp=c(0,-0.1,0))
    if(sharedness.score){
        col.sh.legend <- colorRampPalette(c("white","lightgreen","green",
                                            "green2","green4","darkgreen"))(200)
        color.bar(col.sh.legend,-0.5,1/6*cn)
        axis(4,c(-0.5,1/6*cn),c(0,1),las=1,tck = -0.01)
    }

    ## Axis labelling
    axis(1,1:N,d.hit$GeneID,las = 2, tick = FALSE, ...)
    axis(2,1:cn,condition[cn:1], las = 2, tick = FALSE, ...)
    if(sharedness.score) axis(2,0,"Sharedness", las = 2, tick = FALSE)
}
